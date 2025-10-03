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
#include "htslib/khash.h"
#include "htslib/sam.h"
#include <cstdint>
#include <stdint.h>

// htslib 4-bit-encoding values

static const uint8_t base_to_count_field[16] = {
	[0] = COUNT_N,
	[NT_A] = 0,
	[NT_T] = 1,
	[3] = COUNT_N,
	[NT_C] = 2,
	[5] = COUNT_N,
	[6] = COUNT_N,
	[7] = COUNT_N,
	[NT_G] = 3,
	[9]  = COUNT_N,
	[10] = COUNT_N,
	[11] = COUNT_N,
	[12] = COUNT_N,
	[13] = COUNT_N,
	[14] = COUNT_N,
	[15] = COUNT_N
};



#define PILEUP_PAIR_STRIDE 2
#define PILEUP_PAIR_FLAG 0
#define PILEUP_PAIR_BASE 1

/* LOOKUP TABLES */

// Indexed as: map[is_del][is_head][is_tail]
static const uint8_t bam_pileup_to_flag[2][2][2] = {
	// is_del = 0
	{
		// is_head = 0
		{
			FLAG_UNSET,
			FLAG_TAIL
		},
		// is_head = 1
		{
			FLAG_HEAD,
			FLAG_HEAD | FLAG_TAIL
		}
	},
	// is_del = 1
	{
		// is_head = 0
		{
			FLAG_IS_DEL,
			FLAG_IS_DEL | FLAG_TAIL
		},
		// is_head = 1
		{
			FLAG_IS_DEL | FLAG_HEAD,
			FLAG_IS_DEL | FLAG_HEAD | FLAG_TAIL
		}
	}
};

static const uint8_t indel_to_flag[3] = {
	[0] = FLAG_FDEL,  // negative
	[1] = FLAG_UNSET,         // zero
	[2] = FLAG_FINS   // positive
};

static const uint8_t position_fail_to_flag[2] = {FLAG_UNSET, FLAG_POS_FAIL};
static const uint8_t quality_fail_to_flag[2]  = {FLAG_UNSET, FLAG_QUAL_FAIL};
static const uint8_t reverse_to_flag[2] = {FLAG_UNSET, FLAG_REV};

/* HELPER FUNCTIONS */


static int nttable_get_offset(const NTTable *nttable, const int pos) {
	// TODO: boundary check
	return pos - nttable->params.beg;
}

static int* nttable_get_counts(const NTTable *nttable, const int pos) {
	return nttable->counts + nttable_get_offset(nttable, pos);
}

uint8_t get_pileup_flag(const NTParams &params, const PileupReadInfo &p) {
	return
		bam_pileup_to_flag[p.is_del][p.is_head][p.is_tail] |
		reverse_to_flag[p.rev] |
		indel_to_flag[(p.indel > 0) + (p.indel >= 0)] |
		quality_fail_to_flag[p.base_q <= params.bq_bound] |
		// Position fail (should the test be separate for forward and reverse?)
		position_fail_to_flag[(
			p.qpos < params.head_clip_bound ||
			(p.base_q && p.qlen - p.qpos < params.head_clip_bound)
		)];
}

uint8_t get_pileup_base(const bam_pileup1_t &p) {
	return bam_seqi(bam_get_seq(p.b), p.qpos);
}

uint8_t get_pileup_base_quality(const bam_pileup1_t &p) {
	return bam_get_qual(p.b)[p.qpos];
}

void base_set(
	BaseInfo &b,
	const NTParams &params,
	const PileupReadInfo &p
) {
	b.base = p.base_nt16i;
	b.flag = get_pileup_flag(params, p);
	b.map_quality = p.map_q;
	b.base_quality = p.base_q;
}

int collate_alleles(const NTParams &params, const PileupReadInfo &p, khash_t(strh) *t) {
	// Update read pair summary hash map
	int put_rc;  // return code from put
	const khiter_t i = kh_put(strh, t, p.qname.c_str(), &put_rc);  // n.b. khash does not copy the string, so p must not die
	switch (put_rc) {
		case 0:
			// qname seen => set second read
			base_set(kh_val(t, i).bases[1], params, p);
			break;
		case 1:
		case 2:
			// new qname => set first read
			base_set(kh_val(t, i).bases[0], params, p);
			kh_val(t, i).bases[1].base = UNDEFINED_VALUE;
			break;
		case -1:
			fprintf(stderr, "Failed to put key into khash!\n");
			return 1;
		default:
			fprintf(stderr, "Unknown khash return code: %d!\n", put_rc);
			return 1;
	}
	return 0;
}

void score_single(
	const BaseInfo b,
	const uint64_t c_offset,
	int *counts
) {
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

int score_pair_biased(const BaseInfoPair info, const uint64_t param_length, int *counts) {
	// NOTE: the first item is ALWAYS set, because they are set in order of appearence
	const BaseInfo a = info.bases[0];
	const BaseInfo b = info.bases[1];

	score_single(a, param_length, counts);

	if (b.base != UNDEFINED_VALUE) {

		if (b.base == a.base) {
			return 0;
		}

		score_single(b, param_length, counts);
	}

	return 0;
}




// static int wrap_scorer(const BaseInfoPair info, const uint64_t param_length, int *counts) {
// 	return score_pair_biased(info, param_length, counts);
// 	// return count_balanced(info, param_length, counts);
// }

/* CALLBACK */

int bam2R_pileup_function(const bam_pileup1_t *pileups_ptr, int pos, int n_pileups, NTTable &nttable) {
	if (!(pos >= nttable.params.beg && pos < nttable.params.end)) {
		return 2;  // out of range
	}

	// Collate alleles by read pair
	// TODO: consider persisting the hash map across positions
	//  That would save memory allocations but require clean-up.
	khash_t(strh) *collated_pileup = kh_init(strh);
	for (int pileup_i = 0; pileup_i < n_pileups; pileup_i++) {
		const bam_pileup1_t htspile = *(pileups_ptr + pileup_i);
		auto pinfo = PileupReadInfo::from_pileup(htspile);
		if (collate_alleles(nttable.params, pinfo, collated_pileup)) {
			kh_destroy(strh, collated_pileup);
			return 1;  // fail
		}
	}

	// Count
	int *counts = nttable_get_counts(&nttable, pos);
	BalancedPairCounter scr;
	for (khint_t i = kh_begin(collated_pileup); i != kh_end(collated_pileup); ++i) {
		if (kh_exist(collated_pileup, i)) {
			scr.score_pair(kh_val(collated_pileup, i), (uint64_t)nttable.params.len(), counts);
		}
	}

	kh_destroy(strh, collated_pileup);
	return 0;
}
