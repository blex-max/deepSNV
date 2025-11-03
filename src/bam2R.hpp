#include <string>


void bam2R (std::string bamfile,
            std::string contig,
            int beg,
            int end,
            int *counts,
            int q,
            int mq,
            int head_clip,
            int maxdepth,
            int verbose,
            int mask,
            int keepflag,
            int maxmismatches);
