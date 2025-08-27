// TODO/NOTE: the throwaway bugs are NOT fixed in this version, I'm trying to keep the original logic but separate it into testable components
// NOTE: the only thing that the R-side interface relies on is the counts array - that's the only thing it is necessary to preserve the structure of
// TODO: NOW TEST

#include "bam2r-pileup.hpp"

// char NUCLEOTIDES[] = {'A','T','C','G','*','N','+','-','^','$','Q'};
static constexpr int COUNT_FIELD(char c) {
    switch (c) {
        case 'A': return 0;
        case 'T': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case '*': return 4;
        case 'N': return 5;
        case '+': return 6;
        case '-': return 7;
        case '^': return 8;
        case '$': return 9;
        case 'Q': return 10;
        default:  return -1;
    }
}

// make the data struct in the hot loop without doing costly allocations etc.
// const the output of this in usage to freeze
static PileupRead unsafe_hot_make(const bam_pileup1_t& htspile) {
    const bam1_t* b = htspile.b;  // local alias
    PileupRead out;

    out.qpos       = htspile.qpos;
    out.qname      = bam_get_qname(b);  // borrowed pointer, so this will break if b disappears
    out.qlen       = b->core.l_qseq;
    out.base_nt16i = bam_seqi(bam_get_seq(b), htspile.qpos);
    out.indel      = htspile.indel;
    out.base_q     = bam_get_qual(b)[htspile.qpos];
    out.map_q      = b->core.qual;
    out.rev        = bam_is_rev(b);
    out.is_del     = htspile.is_del;
    out.is_head    = htspile.is_head;
    out.is_tail    = htspile.is_tail;

    return out;
}


// no static, exposed for testing
// TODO: test this
int score_pile(
  const PileupRead& pile,
  int* counts, // ptr to position in counts array where result data should be recorded ( nttable.counts + (int)pos - nttable.beg)
  const NTParams& params,
  khash_t(strh)* overlap_table
) {
  int strand_offset = pile.rev ? params.len() * N_COUNTS_FIELD: 0;
	int put_ret;
  khiter_t kht_i = kh_put(strh, overlap_table, pile.qname, &put_ret);
	uint8_t prev_base;

	if (put_ret == 0) { //Read already processed to get base processed (we only increment if base is different between overlapping read pairs)
		kht_i = kh_get(strh, overlap_table, pile.qname);
		prev_base = kh_val(overlap_table, kht_i);
	} else {
		//Add the value to the hash
		kh_value(overlap_table, kht_i) = pile.base_nt16i;
	}

  {
		if(put_ret == 0 && prev_base == pile.base_nt16i) return -1;
    if (pile.is_tail) counts[strand_offset + params.len() * COUNT_FIELD('$')]++;
    else if (pile.is_head) counts[strand_offset + params.len() * COUNT_FIELD('^')]++;

    if (pile.qpos < params.head_clip_bound || (pile.rev && pile.qlen - pile.qpos < params.head_clip_bound)) {
      counts[strand_offset + params.len() * COUNT_FIELD('N')]++;
    } else {
      if (!pile.is_del) {
        char base_ch = seq_nt16_str[pile.base_nt16i];
        if (pile.base_q > params.bq_bound) {
          counts[strand_offset + params.len() * COUNT_FIELD(base_ch)]++;
        } else {
          counts[strand_offset + params.len() * COUNT_FIELD('N')]++;
        }

        if (pile.indel > 0)
          counts[strand_offset + params.len() * COUNT_FIELD('+')]++;
        else if (pile.indel < 0)
          counts[strand_offset + params.len() * COUNT_FIELD('-')]++;

      } else {
        counts[strand_offset + params.len() * COUNT_FIELD('*')]++;
      }
      counts[strand_offset + params.len() * COUNT_FIELD('Q')] += pile.map_q;
    }
  }

  return 0;
}

void bam2R_pileup_function(const bam_pileup1_t* pileups_ptr, int pos, int n_pileups, NTTable& nttable)
{
  int pileup_i;
	khash_t(strh) *kh_ptr = kh_init(strh);

  if (pos >= nttable.params.beg && pos < nttable.params.end)
  {
    int* counts = nttable.counts + pos - nttable.params.beg;
    for (pileup_i=0; pileup_i<n_pileups; pileup_i++)
    {
      const bam_pileup1_t* htspile = pileups_ptr + pileup_i;
      const PileupRead pile = unsafe_hot_make(*htspile);
      score_pile(pile, counts, nttable.params, kh_ptr);
    }
  }
	kh_destroy(strh, kh_ptr);
}

