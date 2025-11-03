/**********************************************************************
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015-2018 drjsanger@github
 ***********************************************************************/

#include "bam2r-pileup.hpp"
#include <cstring>
#include <stdexcept>

#include "bam2R.hpp"
#include "bounds.hpp"


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


int *bam2R (htsFile *aln_read,
            std::string aln_fp,
            int tid,
            int64_t beg,
            int64_t end,
            int q,
            int mq,
            int head_clip,
            int maxdepth,
            int exclude_flag,
            int keep_flag,
            int maxmismatches) {
    bam_plp_t buf = NULL;
    bam1_t *b = NULL;
    bam_hdr_t *head = NULL;

    const NTParams params{beg - 1, end, q, head_clip};

    safe_size_opts sso; // I wish we were on >C++17
    sso.msg = "error calculating size of counts array";
    size_t rsize = safe_size ((end - beg + 1) * 11 * 2, sso);
    // must use new if to return this array
    int *counts = new int[rsize]; // I'd rather use a vector, but for
                                  // the sake of less change

    NTTable nttable{params, counts, aln_read};

    int64_t maxNM = (maxmismatches != -1) ? maxmismatches : INT64_MAX;
    unsigned long long no_NM_count = 0;

    buf = bam_plp_init (
        0,
        static_cast<void *> (&nttable)); // initialize pileup ALEX:
                                         // why is nttable passed here
    bam_plp_set_maxcnt (buf, maxdepth);
    b = bam_init1();
    // int mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
    // | BAM_FSUPPLEMENTARY;
    int pos, n_plp = -1;
    const bam_pileup1_t *pl;
    hts_idx_t *idx;
    idx = sam_index_load (nttable.in,
                          aln_fp.c_str()); // load BAM index
    if (idx == 0) {
        throw std::runtime_error (
            "BAM/CRAM index file is not available.\n");
    }

    // Implement a fetch style iterator
    hts_itr_t *iter = sam_itr_queryi (idx, tid, nttable.params.beg,
                                      nttable.params.end);
    int result;
    while ((result = sam_itr_next (nttable.in, iter, b)) >= 0) {
        if ((b->core.flag & exclude_flag) == 0 &&
            b->core.qual >= mq) { // as 1.27.1 if these conds only
            // (b->core.flag & *keepflag) == *keepflag &&
            // getNM (b, no_NM_count) <= maxNM) {
            bam_plp_push (buf, b);
        };
        while ((pl = bam_plp_next (buf, &tid, &pos, &n_plp)) != 0) {
            int rc = bam2R_pileup_function (pl, pos, n_plp, nttable);
            if (rc == 1) {
                throw std::runtime_error ("pileup callback failed!");
            }
        }
    }
    if (result < -1) {
        throw std::runtime_error ("Error reading sam iterator.\n");
    }
    sam_itr_destroy (iter);
    hts_idx_destroy (idx);

    bam_plp_push (buf, 0); // finalize pileup

    while ((pl = bam_plp_next (buf, &tid, &pos, &n_plp)) != 0) {
        bam2R_pileup_function (pl, pos, n_plp, nttable);
    }

    if (maxmismatches != -1 && no_NM_count > 0) {
        printf ("%llu reads did not have NM tags; max.mismatches "
                "filter was not "
                "applied to them.\n",
                no_NM_count);
    }

    bam_destroy1 (b);
    bam_hdr_destroy (head);
    bam_plp_destroy (buf);
    hts_close (nttable.in);

    return counts;
}
