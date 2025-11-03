#include <htslib/hts.h>
#include <string>


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
            int maxmismatches);
