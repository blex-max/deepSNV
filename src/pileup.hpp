#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "bounds.hpp"
#include "htslib/sam.h"

#include "const.hpp"
#include "structs.hpp"

// TODO: use the 4 additional bits in the base value for the SET flag
inline constexpr uint8_t FLAG_UNSET = 0;
inline constexpr uint8_t FLAG_POS_FAIL = (1 << 0); // Position fail
inline constexpr uint8_t FLAG_QUAL_FAIL = (1 << 1); // Quality fail
inline constexpr uint8_t FLAG_REV = (1 << 2); // Reverse orientation
inline constexpr uint8_t FLAG_FDEL =
    (1 << 3); // Followed by a deletion
inline constexpr uint8_t FLAG_FINS =
    (1 << 4); // Followed by an insertion
inline constexpr uint8_t FLAG_HEAD = (1 << 5); // Head
inline constexpr uint8_t FLAG_TAIL = (1 << 6); // Tail
inline constexpr uint8_t FLAG_IS_DEL = (1 << 7); // Is a deleted base

// htslib 4-bit-encoding values
inline constexpr uint8_t base_to_count_field[16] = {
    COUNT_N, COUNT_A, COUNT_C, COUNT_N, COUNT_G, COUNT_N,
    COUNT_N, COUNT_N, COUNT_T, COUNT_N, COUNT_N, COUNT_N,
    COUNT_N, COUNT_N, COUNT_N, COUNT_N};


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
        uint8_t nt = static_cast<uint8_t> (
            safe_size (bam_seqi (bam_get_seq (p.b), p.qpos),
                       {0, 15,
                        "unexpected result when accessing base at "
                        "pileup postion"}));
        // clang-format off
        return PileupReadInfo{p.qpos,
                              std::string (bam_get_qname (p.b)),
                              p.b->core.l_qseq,
                              p.b->core.qual,
                              nt,
                              bam_get_qual (p.b)[p.qpos],
                              p.indel,
                              bam_is_rev (p.b),
                              p.is_del != 0,
                              p.is_head != 0,
                              p.is_tail != 0};
        // clang-format on
    }
};

inline uint8_t get_pileup_flag (const count_params &params,
                                const PileupReadInfo &p) {
    /* LOOKUP TABLES */
    // Indexed as: map[is_del][is_head][is_tail]
    constexpr uint8_t EVENT_TO_FLAG[2][2][2] = {
        // is_del = 0
        {// is_head = 0
         {FLAG_UNSET, FLAG_TAIL},
         // is_head = 1
         {FLAG_HEAD, FLAG_HEAD | FLAG_TAIL}},
        // is_del = 1
        {// is_head = 0
         {FLAG_IS_DEL, FLAG_IS_DEL | FLAG_TAIL},
         // is_head = 1
         {FLAG_IS_DEL | FLAG_HEAD,
          FLAG_IS_DEL | FLAG_HEAD | FLAG_TAIL}}};

    constexpr uint8_t INDEL_TO_FLAG[3] = {
        [0] = FLAG_FDEL, // negative
        [1] = FLAG_UNSET, // zero
        [2] = FLAG_FINS // positive
    };
    constexpr uint8_t POS_FAIL_TO_FLAG[2] = {FLAG_UNSET,
                                             FLAG_POS_FAIL};
    constexpr uint8_t QUAL_FAIL_TO_FLAG[2] = {FLAG_UNSET,
                                              FLAG_QUAL_FAIL};
    constexpr uint8_t REV_TO_FLAG[2] = {FLAG_UNSET, FLAG_REV};

    return EVENT_TO_FLAG[p.is_del][p.is_head][p.is_tail] |
        REV_TO_FLAG[p.rev] |
        INDEL_TO_FLAG[(p.indel > 0) + (p.indel >= 0)] |
        QUAL_FAIL_TO_FLAG[p.base_q <= params.min_baseq] |
        // Position fail (should the test be separate for forward and
        // reverse?)
        POS_FAIL_TO_FLAG[(
            p.qpos < params.clip_bound ||
            (p.base_q && p.qlen - p.qpos < params.clip_bound))];
}

struct BaseInfo {
    uint8_t base = UNDEFINED_VALUE;
    uint8_t base_quality;
    uint8_t flag;
    uint8_t map_quality;

    void from_pinfo (const PileupReadInfo &pr,
                     const count_params &params) {
        auto pr_flag = get_pileup_flag (params, pr);
        uint8_t input_base;
        if ((pr_flag & (FLAG_QUAL_FAIL | FLAG_POS_FAIL)) !=
            0) { // squash to ambig if fail
            input_base = UNDEFINED_VALUE;
        } else {
            input_base = pr.base_nt16i;
        }
        this->base = input_base;
        this->base_quality = pr.base_q;
        this->flag = pr_flag;
        this->map_quality = pr.map_q;
    }
};

struct BasePairInfo {
    BaseInfo baseinfo[2];
};

inline void base_set (BaseInfo &b,
                      const count_params &params,
                      const PileupReadInfo &pri) {
    b.base = pri.base_nt16i;
    b.flag = get_pileup_flag (params, pri);
    b.map_quality = pri.map_q;
    b.base_quality = pri.base_q;
}


class AlleleEventCounter {
  private:
    const count_params params;
    std::vector<int> &counts;

  public:
    AlleleEventCounter (const count_params params_,
                        std::vector<int> &counts_)
        : params (params_),
          counts (counts_) {}

    void _collate_alleles (const count_params &par,
                           const PileupReadInfo &pir,
                           std::unordered_map<std::string,
                                              BasePairInfo> &m) {
        // first seen goes into [0], second into [1]
        // n.b. BaseInfoPair ctor inits .base to UNDEFINED_VALUE
        auto emp = m.emplace (
            pir.qname, BasePairInfo{}); // could be more efficient
        auto kv = emp.first;
        // if there was already a key, emplace fails and nothing
        // inserted.
        bool qname_new_to_map = emp.second;
        BasePairInfo &bp = kv->second;

        auto b0 = bp.baseinfo[0].base;
        auto b1 = bp.baseinfo[1].base;

        int to_set;
        if (!qname_new_to_map) { // qname seen before
            if (b0 == UNDEFINED_VALUE || b1 != UNDEFINED_VALUE) {
                throw std::runtime_error ("pair map malformed! " +
                                          pir.qname);
            }
            to_set = 1;
        } else {
            to_set = 0;
        }
        // fill new
        base_set (bp.baseinfo[to_set], par, pir);
    }

    void _score_single (const BaseInfo b,
                        const size_t pos_offset) {
        // field accessor that compiler should inline
        constexpr auto make_idx = [] (const size_t block_offset) {
            return [block_offset] (const size_t field) -> size_t {
                return block_offset + field;
            };
        };
        auto field =
            make_idx ((pos_offset * N_FIELDS_PER_OBS) +
                      ((b.flag & FLAG_REV) ? RSTRAND_OFFSET : 0));


        counts[field (COUNT_HEAD)] +=
            (b.flag & FLAG_HEAD) != FLAG_UNSET;
        counts[field (COUNT_TAIL)] +=
            (b.flag & FLAG_TAIL) != FLAG_UNSET;

        if (b.flag & FLAG_POS_FAIL) {
            counts[field (COUNT_N)]++;
        } else {
            if (b.flag & FLAG_IS_DEL) {
                counts[field (COUNT_IS_DEL)]++;
            } else {
                if (b.flag & FLAG_QUAL_FAIL) {
                    counts[field (COUNT_N)]++;
                } else {
                    // ASSUMPTION: base is 4 bit (in [0, 15])
                    counts[field (base_to_count_field[b.base])]++;
                }

                // NOTE: what about multi-base deletions (is_del
                // follwed by negative indel?)?
                counts[field (COUNT_DEL)] +=
                    (b.flag & FLAG_FDEL) != 0;
                counts[field (COUNT_INS)] +=
                    (b.flag & FLAG_FINS) != 0;
            }
            counts[field (COUNT_MAPQ)] +=
                b.map_quality; // not assessed to be positive, but not
                               // really important for our needs right
                               // now
        }
    }

    void _score_pair (const BasePairInfo &bp,
                      size_t pos_offset,
                      int &pair_toggle) {
        // NOTE: the first item is ALWAYS set, because they are set in
        // order of appearence
        const BaseInfo &a = bp.baseinfo[0];
        const BaseInfo &b = bp.baseinfo[1];

        // if bases different and b is defined
        // count both bases toggling between
        // first and second bases
        bool bundef = (b.base == UNDEFINED_VALUE);
        bool bases_differ = (b.base != a.base);

        if (bundef) {
            _score_single (a, pos_offset);
            return;
        }

        if (bases_differ) {
            _score_single (a, pos_offset);
            _score_single (b, pos_offset);
            return;
        }

        _score_single (pair_toggle ? a : b, pos_offset);
        pair_toggle = !pair_toggle;
    }

    void count_pileup (const bam_pileup1_t *pileups_ptr,
                       const size_t pos_block_offset,
                       const size_t n_reads) {
        // Collate alleles by read pair
        std::unordered_map<std::string, BasePairInfo> qname_map;
        for (size_t i = 0; i < n_reads; ++i) {
            const bam_pileup1_t htspile = *(pileups_ptr + i);
            auto pinfo = PileupReadInfo::from_pileup (htspile);
            _collate_alleles (params, pinfo, qname_map);
        }

        // Count
        int toggle = 0;
        for (auto &[qname, bpair] : qname_map) {
            _score_pair (bpair, pos_block_offset,
                         toggle); // rlen is the wrong thing to pass
        }
    }
};
