#ifndef BAM2R_PILEUP_H
#define BAM2R_PILEUP_H

#include "htslib/hts.h"
#include "htslib/khash.h"
#include "htslib/sam.h"


KHASH_MAP_INIT_STR(strh, uint8_t)


struct NTParams{ // fixed params
  // removed s param as never used
	const int beg, end, bq_boundary, head_clip;

	int len() const noexcept {
	  return end - beg;
	}
};

constexpr int N_COUNTS_FIELD = 11;
// I'm dubious that including the htsFile in this struct is useful
// removed i, some sort of unused counter
struct NTTable {
  const NTParams params;  // separated fixed params into here
	int* const counts;  // data array start
	htsFile *in;
};


// NOTE: from htslib - const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";  (this is indexed to get the base char)
// struct to cleave dependency between pileup function and BAM for testing
// TODO: consider making these all const (frozen) via a constructor (but don't it's not necessary at this time)
struct PileupRead {
  int32_t qpos;
  const char* qname;
  int32_t qlen;
  uint32_t map_q;
  uint8_t base_nt16i;  // 0-15
  uint8_t base_q;
  int indel;
  bool rev, is_del, is_head, is_tail;
};


// exposed for testing
void score_pile(
  const PileupRead& pile,
  int* counts,
  const NTParams& params,
  khash_t(strh)* overlap_table
);


void bam2R_pileup_function(
  const bam_pileup1_t* pileups_ptr,
  int pos,
  int n_pileups,
  NTTable& nttable
);


#endif
