#include <htslib/hts.h>
#include <string>


std::pair<size_t,
          int *>
    bam2R (htsFile *aln_read,
           std::string aln_fp,
           int tid,
           int64_t beg,
           int64_t end,
           int q,
           int mq,
           int head_clip,
           int maxdepth,
           int exclude_flag);
