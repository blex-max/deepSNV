/**********************************************************************
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015-2018 drjsanger@github
 ***********************************************************************/

#include "R_ext/Print.h"
#include "bam2r-pileup.hpp"
#include "deepsnv-prototypes.h"

#define R_NO_REMAP
#include <Rinternals.h>

static inline int64_t getNM (const bam1_t *b,
                             unsigned long long &count) {
    const uint8_t *nm = bam_aux_get (b, "NM");
    if (nm)
        return bam_aux2i (nm);
    else {
        count++;
        return 0; // Dummy NM value that always passes the filter
    }
}

extern "C" {

void bam2R (char **bamfile,
            char **ref,
            int *beg,
            int *end,
            int *counts,
            int *q,
            int *mq,
            int *s,
            int *head_clip,
            int *maxdepth,
            int *verbose,
            int *mask,
            int *keepflag,
            int *maxmismatches) {
    bam_plp_t buf = NULL;
    bam1_t *b = NULL;
    bam_hdr_t *head = NULL;

    const NTParams params{*beg - 1, *end, *q, *head_clip};
    NTTable nttable{params, counts, hts_open (*bamfile, "r")};
    // nttable.s = *s; //Strand (2=both) - does nothing
    // nttable.i = 0;  // does nothing

    int64_t maxNM = (*maxmismatches != -1) ? *maxmismatches : INT64_MAX;
    unsigned long long no_NM_count = 0;

    if (nttable.in == 0) {
        Rf_error ("Fail to open input BAM/CRAM file %s\n", *bamfile);
    }

    buf = bam_plp_init (0, (void *)&nttable); // initialize pileup ALEX: why is nttable passed here
    bam_plp_set_maxcnt (buf, *maxdepth);
    b = bam_init1();
    head = sam_hdr_read (nttable.in);
    if (head == NULL) {
        Rf_error ("failed to get header from alignment file");
    }
    // int mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP |
    // BAM_FSUPPLEMENTARY;
    int tid, pos, n_plp = -1;
    const bam_pileup1_t *pl;
    if (strcmp (*ref, "") == 0) { // if a region is not specified
        // Replicate sampileup functionality (uses above mask without supplementary)
        int ret;
        while ((ret = sam_read1 (nttable.in, head, b)) >= 0) {
            if ((b->core.flag & *mask) == 0 &&
                 b->core.qual >= *mq) {  // as 1.27.1 if these conds only
                // (b->core.flag & *keepflag) == *keepflag &&
                // getNM (b, no_NM_count) <= maxNM) {
                bam_plp_push (buf, b);
            };
            while ((pl = bam_plp_next (buf, &tid, &pos, &n_plp)) != 0) {
                int rc = bam2R_pileup_function (pl, pos, n_plp, nttable);
                if (rc == 1) {
                    Rf_error ("pileup callback failed!");
                }
            }
        }
    } else {
        int tid;
        hts_idx_t *idx;
        idx = sam_index_load (nttable.in, *bamfile); // load BAM index
        if (idx == 0) {
            Rf_error ("BAM/CRAM index file is not available.\n");
        }
        tid = bam_name2id (head, *ref);
        if (tid < 0) {
            Rf_error ("Invalid sequence %s\n", *ref);
        }

        if (*verbose)
            Rprintf ("Reading %s, %s:%d-%d\n", *bamfile, *ref, nttable.params.beg + 1,
                     nttable.params.end);

        // Implement a fetch style iterator
        hts_itr_t *iter = sam_itr_queryi (idx, tid, nttable.params.beg, nttable.params.end);
        int result;
        while ((result = sam_itr_next (nttable.in, iter, b)) >= 0) {
            if ((b->core.flag & *mask) == 0 &&
                 b->core.qual >= *mq) {  // as 1.27.1 if these conds only
                // (b->core.flag & *keepflag) == *keepflag &&
                // getNM (b, no_NM_count) <= maxNM) {
                bam_plp_push (buf, b);
            };
            while ((pl = bam_plp_next (buf, &tid, &pos, &n_plp)) != 0) {
                int rc = bam2R_pileup_function (pl, pos, n_plp, nttable);
                if (rc == 1) {
                    Rf_error ("pileup callback failed!");
                }
            }
        }
        if (result < -1) {
            Rf_error ("Error code (%d) encountered reading sam iterator.\n", result);
        }
        sam_itr_destroy (iter);
        hts_idx_destroy (idx);
    }

    bam_plp_push (buf, 0); // finalize pileup

    while ((pl = bam_plp_next (buf, &tid, &pos, &n_plp)) != 0) {
        bam2R_pileup_function (pl, pos, n_plp, nttable);
    }

    if (*maxmismatches != -1 && no_NM_count > 0) {
        Rf_warning ("%llu reads did not have NM tags; max.mismatches filter was not "
                    "applied to them.\n",
                    no_NM_count);
    }

    bam_destroy1 (b);
    bam_hdr_destroy (head);
    bam_plp_destroy (buf);
    hts_close (nttable.in);
}

} // extern "C"
