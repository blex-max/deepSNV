// luca revisions

// NOTE: The R-side interface relies on only the counts array
// all other structure is only internal representation

/*
if a base fails filters store base as ambiguous (N).
NOTE: initially tried different code for incoming ambiguous (N), and
ambiguous due to these filters (255), but that meant that overlapping
ambiguous bases from read pairs would be double counted since the different
reasons for ambiguity aren't reflected in the counts array. This was
discussed and it was decided that this should not occur. NOTE: this
implementation means an overlapping base encountered after a good base that
fails these filters, and is therefore flipped to N/15, is treated as a
different base and therefore increments the ambiguous counter in the
results array. This has been discussed and determined to be the correct
approach (subject to testing).
*/

/*
Changes:
- all non-canonical bases count as N (any valid 4-bit code)
*/

#include "bam2r-pileup.hpp"
#include "htslib/sam.h"
#include <cstdint>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <unordered_map>

// htslib 4-bit-encoding values

static const uint8_t base_to_count_field[16] = {
    COUNT_N, COUNT_A, COUNT_C, COUNT_N, COUNT_G, COUNT_N, COUNT_N, COUNT_N,
    COUNT_T, COUNT_N, COUNT_N, COUNT_N, COUNT_N, COUNT_N, COUNT_N, COUNT_N};


#define PILEUP_PAIR_STRIDE 2
#define PILEUP_PAIR_FLAG 0
#define PILEUP_PAIR_BASE 1

/* LOOKUP TABLES */

// Indexed as: map[is_del][is_head][is_tail]
static const uint8_t bam_pileup_to_flag[2][2][2] = {
    // is_del = 0
    {// is_head = 0
     {FLAG_UNSET, FLAG_TAIL},
     // is_head = 1
     {FLAG_HEAD, FLAG_HEAD | FLAG_TAIL}},
    // is_del = 1
    {// is_head = 0
     {FLAG_IS_DEL, FLAG_IS_DEL | FLAG_TAIL},
     // is_head = 1
     {FLAG_IS_DEL | FLAG_HEAD, FLAG_IS_DEL | FLAG_HEAD | FLAG_TAIL}}};

static const uint8_t indel_to_flag[3] = {
    [0] = FLAG_FDEL, // negative
    [1] = FLAG_UNSET, // zero
    [2] = FLAG_FINS // positive
};

static const uint8_t position_fail_to_flag[2] = {FLAG_UNSET, FLAG_POS_FAIL};
static const uint8_t quality_fail_to_flag[2] = {FLAG_UNSET, FLAG_QUAL_FAIL};
static const uint8_t reverse_to_flag[2] = {FLAG_UNSET, FLAG_REV};

/* HELPER FUNCTIONS */


static int nttable_get_offset (const NTTable *nttable,
                               const int pos) {
    // TODO: boundary check
    return pos - nttable->params.beg;
}

static int *nttable_get_counts (const NTTable *nttable,
                                const int pos) {
    return nttable->counts + nttable_get_offset (nttable, pos);
}

uint8_t get_pileup_flag (const NTParams &params,
                         const PileupReadInfo &p) {
    return bam_pileup_to_flag[p.is_del][p.is_head][p.is_tail] | reverse_to_flag[p.rev] |
        indel_to_flag[(p.indel > 0) + (p.indel >= 0)] |
        quality_fail_to_flag[p.base_q <= params.bq_bound] |
        // Position fail (should the test be separate for forward and reverse?)
        position_fail_to_flag[(p.qpos < params.head_clip_bound ||
                               (p.base_q && p.qlen - p.qpos < params.head_clip_bound))];
}

uint8_t get_pileup_base (const bam_pileup1_t &p) { return bam_seqi (bam_get_seq (p.b), p.qpos); }

uint8_t get_pileup_base_quality (const bam_pileup1_t &p) { return bam_get_qual (p.b)[p.qpos]; }

void base_set (BaseInfo &b,
               const NTParams &params,
               const PileupReadInfo &p) {
    b.base = p.base_nt16i;
    b.flag = get_pileup_flag (params, p);
    b.map_quality = p.map_q;
    b.base_quality = p.base_q;
}

void collate_alleles (const NTParams &params,
                      const PileupReadInfo &p,
                      std::unordered_map<std::string,
                                         BaseInfoPair> &m) {
    // forward member goes into bases[0], reverse into bases[1]
    const int to_set = (int)p.rev; // indexes into the base values
    const int other = 1 - to_set;

    // n.b. BaseInfoPair ctor inits .base to UNDEFINED_VALUE
    auto emp = m.emplace (p.qname, BaseInfoPair{}); // could be more efficient
    auto kv = emp.first;
    // if there was already a key, emplace fails and nothing inserted.
    bool qname_new_to_map = emp.second;
    BaseInfoPair &pi = kv->second;

    auto b_toset = pi.bases[to_set].base;
    auto b_oth = pi.bases[other].base;
    // std::string debug_toset = (b_toset == UNDEFINED_VALUE) ? "UNDEF" : std::string(&seq_nt16_str[b_toset]);
    // std::string debug_oth = (b_oth == UNDEFINED_VALUE) ? "UNDEF" : std::string(&seq_nt16_str[b_oth]);

    if (!qname_new_to_map) { // qname seen before
        if (b_toset != UNDEFINED_VALUE) {
            throw std::runtime_error ("duplicate qname on same strand! " + p.qname);
        }
        if (b_oth == UNDEFINED_VALUE) {
            throw std::runtime_error ("pair map malformed (other strand unset): " + p.qname);
        }
    }
    // fill new
    base_set (pi.bases[to_set], params, p);
}


void score_single (const BaseInfo b,
                   const uint64_t c_offset,
                   int *counts) {
    const uint64_t strand_offset = (b.flag & FLAG_REV) ? c_offset * N_COUNTS_FIELD : 0;

    counts[strand_offset + c_offset * COUNT_HEAD] += (b.flag & FLAG_HEAD) != FLAG_UNSET;
    counts[strand_offset + c_offset * COUNT_TAIL] += (b.flag & FLAG_TAIL) != FLAG_UNSET;

    if (b.flag & FLAG_POS_FAIL) {
        counts[strand_offset + c_offset * COUNT_N]++;
    } else {
        if (b.flag & FLAG_IS_DEL) {
            counts[strand_offset + c_offset * COUNT_IS_DEL]++;
        } else {
            if (b.flag & FLAG_QUAL_FAIL) {
                counts[strand_offset + c_offset * COUNT_N]++;
            } else {
                // ASSUMPTION: base is 4 bit (in [0, 15])
                counts[strand_offset + c_offset * (uint64_t)base_to_count_field[b.base]]++;
            }

            // NOTE: what about multi-base deletions (is_del follwed by negative indel?)?
            counts[strand_offset + c_offset * COUNT_DEL] += (b.flag & FLAG_FDEL) != 0;
            counts[strand_offset + c_offset * COUNT_INS] += (b.flag & FLAG_FINS) != 0;
        }
        counts[strand_offset + c_offset * COUNT_QUALITY] += (int)b.map_quality;
    }
}


void BalancedPairCounter::score_pair (const BaseInfoPair info,
                                      const uint64_t param_length,
                                      int *counts) {
    // NOTE: the first item is ALWAYS set, because they are set in order of appearence
    const BaseInfo a = info.bases[0];
    const BaseInfo b = info.bases[1];

    if (b.base != UNDEFINED_VALUE) {
        if (b.base == a.base) {
            // Same base => alternate between counting one or the other (effect on strand bias)
            if (pair_toggle) {
                score_single (a, param_length, counts);
            } else {
                score_single (b, param_length, counts);
            }
            pair_toggle = !pair_toggle;

        } else {
            // Different bases => count both
            score_single (a, param_length, counts);
            score_single (b, param_length, counts);
        }

    } else {
        // Only one base
        score_single (a, param_length, counts);
    }
}


int bam2R_pileup_function (const bam_pileup1_t *pileups_ptr,
                           int pos,
                           int n_pileups,
                           NTTable &nttable) {
    if (!(pos >= nttable.params.beg && pos < nttable.params.end)) {
        return 2; // out of range
    }

    // Collate alleles by read pair
    std::unordered_map<std::string, BaseInfoPair> qname_map;
    for (int pileup_i = 0; pileup_i < n_pileups; pileup_i++) {
        const bam_pileup1_t htspile = *(pileups_ptr + pileup_i);
        auto pinfo = PileupReadInfo::from_pileup (htspile);
        try {
            collate_alleles (nttable.params, pinfo, qname_map);
        } catch (std::exception &e) {
            throw;
            // return 1; // fail
        }
    }

    // Count
    int *counts = nttable_get_counts (&nttable, pos);
    BalancedPairCounter scr;
    for (auto &[qname, bpair] : qname_map) {
        scr.score_pair (bpair, (uint64_t)nttable.params.len(), counts);
    }

    return 0;
}
