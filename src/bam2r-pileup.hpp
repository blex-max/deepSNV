// luca revisions

#ifndef BAM2R_PILEUP_H
#define BAM2R_PILEUP_H

#include "htslib/hts.h"
#include "htslib/sam.h"
#include <cstdint>
#include <stdexcept>
#include <string>

#define NT_A 1
#define NT_C 2
#define NT_G 4
#define NT_T 8
#define AMBIG_NT 255

#define COUNT_IS_DEL 4 // *
#define COUNT_N 5 // N
#define COUNT_INS 6 // +
#define COUNT_DEL 7 // -
#define COUNT_HEAD 8 // ^
#define COUNT_TAIL 9 // $
#define COUNT_QUALITY 10 // Q

// TODO: use the 4 additional bits in the base value for the SET flag
#define FLAG_UNSET 0
#define FLAG_POS_FAIL (1 << 0) // Position fail
#define FLAG_QUAL_FAIL (1 << 1) // Quality fail
#define FLAG_REV (1 << 2) // Reverse orientation
#define FLAG_FDEL (1 << 3) // Followed by a deletion
#define FLAG_FINS (1 << 4) // Followed by an insertion
#define FLAG_HEAD (1 << 5) // Head
#define FLAG_TAIL (1 << 6) // Tail
#define FLAG_IS_DEL (1 << 7) // Is a deleted base

#define UNDEFINED_VALUE UINT8_MAX


struct NTParams {
    // fixed params
    // removed s param as never used
    const int beg, end, bq_bound, head_clip_bound;

    NTParams (int beg,
              int end,
              int bq_bound,
              int head_clip_bound)
        : beg (beg),
          end (end),
          bq_bound (bq_bound),
          head_clip_bound (head_clip_bound) {
        if (!(end > beg)) {
            throw std::invalid_argument ("end must be greater than beginning");
        }
    }

    int len () const noexcept { return end - beg; }
};

constexpr int N_COUNTS_FIELD = 11;
// I'm dubious that including the htsFile in this struct is useful
// removed i, some sort of unused counter
struct NTTable {
    const NTParams params; // separated fixed params into here
    int *const counts; // data array start
    htsFile *in;
};

uint8_t get_pileup_base (const bam_pileup1_t &p);

uint8_t get_pileup_base_quality (const bam_pileup1_t &p);


// cleave tie to bam pointer
struct PileupReadInfo {
    int32_t qpos;
    std::string qname;
    int32_t qlen;
    uint8_t map_q;
    uint8_t base_nt16i; // 0-15
    uint8_t base_q;
    int indel;
    bool rev, is_del, is_head, is_tail;

    static PileupReadInfo from_pileup (const bam_pileup1_t &p) {
        return PileupReadInfo{p.qpos,
                              std::string (bam_get_qname (p.b)),
                              p.b->core.l_qseq,
                              p.b->core.qual,
                              get_pileup_base (p),
                              get_pileup_base_quality (p),
                              p.indel,
                              bam_is_rev (p.b),
                              p.is_del != 0,
                              p.is_head != 0,
                              p.is_tail != 0};
    }
};


uint8_t get_pileup_flag (const NTParams &params,
                         const PileupReadInfo &p);

struct BaseInfo {
    uint8_t base = UNDEFINED_VALUE;
    uint8_t base_quality;
    uint8_t flag;
    uint8_t map_quality;

    void from_pinfo (const PileupReadInfo &p,
                     const NTParams &params) {
        auto pri_flag = get_pileup_flag (params, p);
        uint8_t input_base;
        if ((pri_flag & (FLAG_QUAL_FAIL | FLAG_POS_FAIL)) != 0) { // squash to ambig if fail
            input_base = AMBIG_NT;
        } else {
            input_base = p.base_nt16i;
        }
        this->base = input_base;
        this->base_quality = p.base_q;
        this->flag = pri_flag;
        this->map_quality = p.map_q;
    }
};

struct BaseInfoPair {
    BaseInfo bases[2];
};

class BalancedPairCounter {
    /*
	If calls are equivalent (ish), alternate between counting first seen and second seen.
	Does not directly account for strand
	Does not consider equivalence beyond same base - though all qual/pos fails are converted
	into ambiguous bases upstream so the better call will always be taken.
	But no opinion on is_tail, is_head, etc.
	*/
    bool pair_toggle = true;

  public:
    void score_pair (const BaseInfoPair info,
                     const uint64_t param_length,
                     int *counts);
};

int bam2R_pileup_function (const bam_pileup1_t *pileups_ptr,
                           int pos,
                           int n_pileups,
                           NTTable &nttable);

#endif
