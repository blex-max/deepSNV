/**********************************************************************
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015-2018 drjsanger@github
 ***********************************************************************/

#include "bam2r_pileup.hpp"

#ifdef STANDALONE_TEST  // placeholder for IDE, not actually relevant to this file
#include <cstdio>
#include <cstdlib>

#define Rprintf std::printf
#define Rf_error(...)                   \
  do {                                  \
    std::fprintf(stderr, __VA_ARGS__);  \
    std::fputc('\n', stderr);           \
    std::exit(1);                       \
  } while (0)
#define Rf_warning(...)                 \
  do {                                  \
    std::fprintf(stderr, __VA_ARGS__);  \
    std::fputc('\n', stderr);           \
  } while (0)

using DL_FUNC = void (*)();
struct DllInfo {
  int _unused;
};
struct R_CMethodDef {
  const char *name;
  DL_FUNC fun;
  int numArgs;
};
inline void R_registerRoutines(DllInfo *, R_CMethodDef *, void *, void *,
                               void *) {}
inline void R_useDynamicSymbols(DllInfo *, int) {}

#else
#define R_NO_REMAP
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#endif


static inline int64_t getNM(const bam1_t *b, unsigned long long& count)
{
	const uint8_t *nm = bam_aux_get(b, "NM");
	if (nm)
		return bam_aux2i(nm);
	else {
		count++;
		return 0;  // Dummy NM value that always passes the filter
	}
}

extern "C" {

int bam2R(char** bamfile, char** ref, int* beg, int* end, int* counts, int* q, int* mq, int* s, 
          int* head_clip, int* maxdepth, int* verbose, int* mask, int *keepflag, int *maxmismatches )
{

	bam_plp_t buf = NULL;
	bam1_t *b = NULL;
	bam_hdr_t *head = NULL;

  const NTParams params{*beg-1, *end, *q, *head_clip};
	NTTable nttable{params, counts, hts_open(*bamfile, "r")};
	// nttable.s = *s; //Strand (2=both) - does nothing
	// nttable.i = 0;  // does nothing

	int64_t maxNM = (*maxmismatches != -1)? *maxmismatches : INT64_MAX;
	unsigned long long no_NM_count = 0;

	if (nttable.in == 0) {
		Rf_error("Fail to open input BAM/CRAM file %s\n", *bamfile);
		return 1;
	}

	buf = bam_plp_init(0,(void *)&nttable); // initialize pileup
	bam_plp_set_maxcnt(buf,*maxdepth);
	b = bam_init1();
	//get header
	head = sam_hdr_read(nttable.in);
	//int mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY;
  int tid, pos, n_plp = -1;
	const bam_pileup1_t *pl;

	if (strcmp(*ref, "") == 0) { // if a region is not specified
		//Replicate sampileup functionality (uses above mask without supplementary)
		int ret;
		while((ret = sam_read1(nttable.in, head, b)) >= 0){
			if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && getNM(b, no_NM_count) <= maxNM) {
					bam_plp_push(buf, b);
            };
			while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
				bam2R_pileup_function(pl,pos,n_plp,nttable);
			}
		}
	}
	else {
		int tid;
		hts_idx_t *idx;
		idx = sam_index_load(nttable.in,*bamfile); // load BAM index
		if (idx == 0) {
			Rf_error("BAM/CRAM index file is not available.\n");
			return 1;
		}
		tid = bam_name2id(head, *ref);
		if (tid < 0) {
			Rf_error("Invalid sequence %s\n", *ref);
			return 1;
		}

		if(*verbose)
			Rprintf("Reading %s, %s:%d-%d\n", *bamfile, *ref, nttable.params.beg+1, nttable.params.end);

		//Implement a fetch style iterator
		hts_itr_t *iter = sam_itr_queryi(idx, tid, nttable.params.beg, nttable.params.end);
		int result;
		while ((result = sam_itr_next(nttable.in, iter, b)) >= 0) {
			if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && getNM(b, no_NM_count) <= maxNM) {
				bam_plp_push(buf, b);
			};
			while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
				bam2R_pileup_function(pl, pos, n_plp, nttable);
			}
		}
    if(result < -1){
      Rf_error("Error code (%d) encountered reading sam iterator.\n", result);
			return 1;
    }
		sam_itr_destroy(iter);
		hts_idx_destroy(idx);
	}

	bam_plp_push(buf,0); // finalize pileup

  while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
    bam2R_pileup_function(pl, pos, n_plp, nttable);
  }

	if (*maxmismatches != -1 && no_NM_count > 0) {
		Rf_warning("%llu reads did not have NM tags; max.mismatches filter was not applied to them.\n", no_NM_count);
	}

	bam_destroy1(b);
	bam_hdr_destroy(head);
	bam_plp_destroy(buf);
	hts_close(nttable.in);
	return 0;
}

R_CMethodDef cMethods[] = {
		{"bam2R", (DL_FUNC) &bam2R, 12}
};

void R_init_bam2R(DllInfo *info) {
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

} // extern "C"
