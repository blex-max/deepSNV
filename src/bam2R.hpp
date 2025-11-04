#include <cstdint>
#include <cstring>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdexcept>
#include <vector>

#include "bounds.hpp"
#include "pileup.hpp"
#include "structs.hpp"


// static inline int64_t getNM (const bam1_t *b,
//                              unsigned long long &count) {
//     const uint8_t *nm = bam_aux_get (b, "NM");
//     if (nm)
//         return bam_aux2i (nm);
//     else {
//         count++;
//         return 0; // Dummy NM value that always passes the filter
//     }
// }

// bam2R
inline void read_and_count (htsFile *aln_fr,
                     hts_idx_t *aln_idx,
                     const hts_region &reg,
                     const count_params &params,
                     std::vector<int> &counts
                     // int keep_flag,
                     // int maxmismatches
) {
    bam_plp_t buf = NULL;
    bam1_t *b = NULL;
    bam_hdr_t *head = NULL;
    AlleleEventCounter aev (reg, params, counts);

    // int64_t maxNM = (maxmismatches != -1) ? maxmismatches :
    // INT64_MAX; unsigned long long no_NM_count = 0;

    buf = bam_plp_init (0,
                        NULL); // initialize pileup
    bam_plp_set_maxcnt (buf, params.max_depth);
    b = bam_init1();
    // int mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
    // | BAM_FSUPPLEMENTARY;
    int plp_tid = -1;
    int64_t plp_pos = -1;
    int n_plp = -1;
    const bam_pileup1_t *pl;

    // Implement a fetch style iterator
    hts_itr_t *iter =
        sam_itr_queryi (aln_idx, reg.rid, reg.start, reg.end);
    int result;
    while ((result = sam_itr_next (aln_fr, iter, b)) >= 0) {
        if ((b->core.flag & params.exclude_flag) == 0 &&
            b->core.qual >=
                params.min_mapq) { // as 1.27.1 if these conds only
            // (b->core.flag & *keepflag) == *keepflag &&
            // getNM (b, no_NM_count) <= maxNM) {
            bam_plp_push (buf, b);
        };
        while ((pl = bam_plp64_next (buf, &plp_tid, &plp_pos,
                                     &n_plp)) != NULL) {
            if (n_plp < 0 || plp_tid < 0 || plp_pos < 0) {
                throw std::runtime_error ("pileup failed");
            }
            aev.count_position (pl, plp_pos, safe_size (n_plp));
        }
    }
    if (result < -1) {
        throw std::runtime_error ("Error reading sam iterator.\n");
    }
    sam_itr_destroy (iter);

    bam_plp_push (buf, 0); // finalize pileup
    while ((pl = bam_plp64_next (buf, &plp_tid, &plp_pos, &n_plp)) !=
           NULL) {
        if (n_plp < 0) {
            throw std::runtime_error ("pileup flush failed");
        }
        aev.count_position (pl, plp_pos, safe_size (n_plp));
    }

    // if (maxmismatches != -1 && no_NM_count > 0) {
    //     printf ("%llu reads did not have NM tags; max.mismatches "
    //             "filter was not "
    //             "applied to them.\n",
    //             no_NM_count);
    // }

    bam_destroy1 (b);
    bam_hdr_destroy (head);
    bam_plp_destroy (buf);
}
