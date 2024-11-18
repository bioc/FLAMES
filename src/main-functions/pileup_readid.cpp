#include "htslib/sam.h"
#include <Rcpp.h>
#include <array>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>
// [[Rcpp::plugins(cpp17)]]

// x <- FLAMES:::variant_count_tb(
//   bam_path =
//   '/vast/scratch/users/wang.ch/RaCHseq_sup_fastq/FLAMES_out/sample15_align2genome_filtered.bam',
//   seqname = 'chr21', pos = 34880579, indel = F, verbose = F)

const std::vector<char> BASES = {'A', 'T', 'C', 'G', '-'};
constexpr int BASES_SIZE = 5;
constexpr int MAX_EDIT_DISTANCE = 3;

// Code for fast edit distance calculation for short sequences modified from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
// s2 is always assumned to be the shorter string (barcode)
bool within_edit_dist(const std::string &s1, const std::string &s2) {

  const std::string_view s1_view(s1);
  const std::string_view s2_view(s2);

  const std::size_t len1 = s1_view.size() + 1;
  const std::size_t len2 = s2_view.size() + 1;

  std::vector<unsigned int> dist_holder(len1 * len2);
  // initialise the edit distance matrix.
  // penalise for gaps at the start and end of the shorter sequence (j)
  // but not for shifting the start/end of the longer sequence (i,0)
  dist_holder[0] = 0; //[0][0]
  for (std::size_t j = 1; j < len2; ++j)
    dist_holder[j] = j; //[0][j];
  for (std::size_t i = 1; i < len1; ++i)
    dist_holder[i * len2] = i; //[i][0];

  // loop over the distance matrix elements and calculate running distance
  for (std::size_t j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (std::size_t i = 1; i < len1; ++i) {
      unsigned int sub =
          (s1_view[i - 1] == s2_view[j - 1]) ? 0 : 1; // match / mismatch score

      const unsigned int &top_left = dist_holder[(i - 1) * len2 + (j - 1)];
      const unsigned int &left = dist_holder[i * len2 + (j - 1)];
      const unsigned int &top = dist_holder[(i - 1) * len2 + j];

      unsigned int min_value = std::min({top + 1, left + 1, top_left + sub});
      dist_holder[i * len2 + j] = min_value;

      if (min_value <= MAX_EDIT_DISTANCE)
        any_below_threshold = true;
    }
    if (!any_below_threshold) { // early exit to save time.
      return false;
    }
  }
  return dist_holder.back() <= MAX_EDIT_DISTANCE;
}

// barcode -> UMI -> array of allele counts (A, T, C, G, -)
using SNPTable = std::unordered_map<
    std::string,
    std::unordered_map<std::string, std::array<unsigned int, BASES_SIZE>>>;

// barcode -> UMI -> variant -> count
using IndelTable = std::unordered_map<
    std::string,
    std::unordered_map<std::string,
                       std::unordered_map<std::string, unsigned int>>>;

// for each barcode, group UMIs within a certain edit distance
// only compare each UMI to the first UMI in the group
template <typename TableType>
std::unordered_map<std::string, std::vector<std::vector<std::string>>>
group_umis(const TableType &table, int max_distance) {

  std::unordered_map<std::string, std::vector<std::vector<std::string>>>
      grouped_umis;
  unsigned int end;

  for (const auto &barcode :
       table) { // first: barcode, second: UMI unordered_map
    for (const auto &umi : barcode.second) { // first: UMI, second: SNP / indel
      // first time seeing this barcode
      if (grouped_umis.find(barcode.first) == grouped_umis.end()) {
        grouped_umis[barcode.first].push_back({umi.first});
        continue;
      }
      // if UMI matches with any eixsting group's first UMI, add to the group
      bool added = false;
      for (auto &group : grouped_umis[barcode.first]) {
        if (within_edit_dist(umi.first, group[0])) {
          group.push_back(umi.first);
          added = true;
          break;
        }
      }
      // start a new group if no match with any existing group
      if (!added) {
        grouped_umis[barcode.first].push_back({umi.first});
      }
    }
  }

  return grouped_umis;
}

unsigned int grouped_umi_voting_snp(
    const std::vector<std::string> &umis,
    const std::unordered_map<std::string, std::array<unsigned int, BASES_SIZE>>
        &subtable) {
  // sum up the array of 5 for all umis
  std::array<unsigned int, BASES_SIZE> umi_sum = {0};
  for (const auto &umi : umis) {
    for (int i = 0; i < BASES_SIZE; i++) {
      umi_sum[i] += subtable.at(umi)[i];
    }
  }

  // find the index of the max element and return
  // TODO: what if there are multiple max elements?
  return std::distance(umi_sum.begin(),
                       std::max_element(umi_sum.begin(), umi_sum.end()));
}

unsigned int grouped_umi_voting_indel(
    const std::vector<std::string> &umis,
    const std::unordered_map<std::string,
                             std::unordered_map<std::string, unsigned int>>
        subtable,
    const std::unordered_set<std::string> &indel_keys) {

  // sum up counts to an unordered_map
  std::unordered_map<std::string, unsigned int> umi_sum;
  for (const auto &umi : umis) {
    for (const auto &indel : subtable.at(umi)) {
      umi_sum[indel.first] += indel.second;
    }
  }

  // find the index of the max element and return
  std::string max_key = "."; // should not matter
  unsigned int max_val = 0;
  for (const auto &indel : umi_sum) {
    if (indel.second > max_val) {
      max_val = indel.second;
      max_key = indel.first;
    }
  }
  return std::distance(indel_keys.begin(), indel_keys.find(max_key));
}

typedef struct plpconf {
  const char *inname;
  samFile *infile;
  sam_hdr_t *in_samhdr;
  hts_idx_t *in_idx;
  hts_itr_t *iter;
} plpconf;

// copied from htslib/samples/pileup.c
// not sure what these are for
inline int plpconstructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  return 0;
}
inline int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  return 0;
}

int readdata(void *data, bam1_t *b) {
  plpconf *conf = (plpconf *)data;
  if (!conf || !conf->infile) {
    return -2; // cant read data
  }

  bool skip = true;
  int ret = 0;
  do {
    ret = conf->iter ? sam_itr_next(conf->infile, conf->iter, b)
                     : sam_read1(conf->infile, conf->in_samhdr, b);
    if (b->core.tid < 0 ||
        (b->core.flag & BAM_FUNMAP)) { // exclude unmapped reads
      skip = true;
      continue;
    } else {
      skip = false;
    }
  } while (skip);

  return ret;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix variant_count_matrix_cpp(Rcpp::String bam_path,
                                             Rcpp::String seqname, int pos,
                                             bool indel, bool verbose) {
  if (pos > 1) {
    // htslib is 0-based
    pos = pos - 1;
  }
  plpconf conf = {0};
  bam_plp_t plpiter = NULL;
  const bam_pileup1_t *plp = NULL;

  conf.inname = bam_path.get_cstring();
  const char *seqname_c = seqname.get_cstring();

  if (!(conf.infile = sam_open(conf.inname, "r"))) {
    Rcpp::stop("Failed to open file %s\n", conf.inname);
  }
  if (!(conf.in_samhdr = sam_hdr_read(conf.infile))) {
    Rcpp::stop("Failed to read header from file!\n");
  }
  if (!(plpiter = bam_plp_init(readdata, &conf))) {
    Rcpp::stop("Failed to initialize pileup data\n");
  }
  bam_plp_set_maxcnt(plpiter, INT_MAX); // caps at 2b
  if (!(conf.in_idx = sam_index_load(conf.infile, conf.inname))) {
    Rcpp::stop("Failed to load index for %s", conf.inname);
  }
  int tid = bam_name2id(conf.in_samhdr, seqname_c);
  if ((conf.iter = sam_itr_queryi(conf.in_idx, tid, pos, pos + 1)) == 0) {
    Rcpp::stop("Failed to parse region");
  }

  // set constructor destructor callbacks
  bam_plp_constructor(plpiter, plpconstructor);
  bam_plp_destructor(plpiter, plpdestructor);

  // char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
  const char seq_nt16_char[] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                                'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

  // fetch an alignment to determine the UMI length first
  bam1_t *bam_tmp = NULL;
  if (!(bam_tmp = bam_init1())) {
    Rcpp::stop("Failed to initialize bamdata");
  }
  if (sam_itr_next(conf.infile, sam_itr_queryi(conf.in_idx, tid, pos, pos + 1),
                   bam_tmp) < 0) {
    Rcpp::stop("Failed to fetch an alignment");
  }
  std::string read_id = bam_get_qname(bam_tmp);
  const std::size_t umi_idx = read_id.find("#");
  const std::size_t id_idx = read_id.find("_");
  if (umi_idx == std::string::npos || id_idx == std::string::npos) {
    Rcpp::stop("Unexpected read id format: %s", read_id);
  }
  if (verbose) {
    Rcpp::Rcout << "Checking read ID format:\n"
                << "ReadID: " << read_id << "\n"
                << "barcode: " << read_id.substr(0, id_idx) << "\n"
                << "UMI: " << read_id.substr(id_idx + 1, umi_idx - id_idx - 1)
                << "\n";
  }
  bam_destroy1(bam_tmp);

  // all values in BASES, repsctively
  std::unordered_map<char, int> allele_to_idx;
  std::unordered_map<int, std::string> idx_to_allele;
  for (int i = 0; i < BASES.size(); i++) {
    allele_to_idx[BASES[i]] = i;
    idx_to_allele[i] = BASES[i];
  }

  SNPTable snps;
  IndelTable indels;

  // all barcodes
  std::set<std::string> barcodes;

  // keep track of all possible indel variants
  std::unordered_set<std::string> indel_keys;
  indel_keys.insert(".");

  int actual_tid = -1, n = -1, actual_pos = -1;
  while ((plp = bam_plp_auto(plpiter, &actual_tid, &actual_pos, &n))) {
    if (actual_tid != tid || actual_pos != pos) {
      continue;
    }

    // iterate all reads and print the read id and the allele
    for (int j = 0; j < n; j++) {

      if (plp[j].is_refskip) {
        continue;
      }

      const std::string read_id = bam_get_qname(plp[j].b);
      const std::string barcode = read_id.substr(0, id_idx);
      const std::string umi = read_id.substr(id_idx + 1, umi_idx - id_idx - 1);
      barcodes.insert(barcode);

      if (!indel) { // count SNV
        if (plp[j].is_del) {
          snps[barcode][umi][allele_to_idx['-']]++;
        } else {
          char allele =
              seq_nt16_char[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)];
          snps[barcode][umi][allele_to_idx[allele]]++;
        }
      } else { // count indel
        if (plp[j].indel == 0) {
          indels[barcode][umi]["."]++;
        } else if (plp[j].indel > 0) {
          // key as '+[inserted_bases]' e.g. '+A' or '+AT'
          std::string inserted_bases = "";
          for (int k = 1; k <= plp[j].indel; k++) {
            inserted_bases +=
                seq_nt16_char[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos + k)];
          }
          indels[barcode][umi]["+" + inserted_bases]++;
          indel_keys.insert("+" + inserted_bases);
        } else if (plp[j].indel < 0) {
          // key as '-[deleted_number_of_bases]' e.g. '-2'
          indels[barcode][umi]["-" + std::to_string(-plp[j].indel)]++;
          indel_keys.insert("-" + std::to_string(-plp[j].indel));
        }
      }
    }
  }

  Rcpp::NumericMatrix ret(indel ? indel_keys.size() : BASES.size(),
                          barcodes.size());
  Rcpp::colnames(ret) = Rcpp::wrap(barcodes);
  if (indel) {
    Rcpp::rownames(ret) = Rcpp::wrap(indel_keys);
  } else {
    Rcpp::rownames(ret) = Rcpp::wrap(BASES);
  }

  std::unordered_map<std::string, std::vector<std::vector<std::string>>>
      grouped_umis;
  if (indel) {
    grouped_umis = group_umis(indels, MAX_EDIT_DISTANCE);
  } else {
    grouped_umis = group_umis(snps, MAX_EDIT_DISTANCE);
  }

  // perform majority voting for reads from the same UMI group and barcode
  if (!indel) {
    for (const auto &barcode : grouped_umis) {
      for (const auto &umi_group : barcode.second) {
        unsigned int max_idx =
            grouped_umi_voting_snp(umi_group, snps[barcode.first]);
        auto bc_it = std::find(barcodes.begin(), barcodes.end(), barcode.first);
        if (bc_it == barcodes.end()) {
          Rcpp::stop("Barcode not found in barcodes");
        } else {
          unsigned int bc_idx = std::distance(barcodes.begin(), bc_it);
          ret(max_idx, bc_idx)++;
        }
      }
    }
  } else {
    for (const auto &barcode : grouped_umis) {
      for (const auto &umi_group : barcode.second) {
        unsigned int max_idx = grouped_umi_voting_indel(
            umi_group, indels[barcode.first], indel_keys);
        auto bc_it = std::find(barcodes.begin(), barcodes.end(), barcode.first);
        if (bc_it == barcodes.end()) {
          Rcpp::stop("Barcode not found in barcodes");
        } else {
          unsigned int bc_idx = std::distance(barcodes.begin(), bc_it);
          ret(max_idx, bc_idx)++;
        }
      }
    }
  }

  if (conf.in_samhdr) {
    sam_hdr_destroy(conf.in_samhdr);
  }
  if (conf.infile) {
    sam_close(conf.infile);
  }
  if (conf.in_idx) {
    hts_idx_destroy(conf.in_idx);
  }
  if (conf.iter) {
    hts_itr_destroy(conf.iter);
  }
  if (plpiter) {
    bam_plp_destroy(plpiter);
  }

  return ret;
}
