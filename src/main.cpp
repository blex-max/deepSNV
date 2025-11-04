#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "bounds.hpp"
#include "const.hpp"
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
inline void count (htsFile *aln_fh,
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

    // fetch all reads overlapping the query region;
    // then do a pileup per base for the total region
    // covered by those retrieved reads;
    // then count events on those pileups which overlap
    // the original query region.
    hts_itr_t *iter =
        sam_itr_queryi (aln_idx, reg.rid, reg.start, reg.end);
    int result;
    while ((result = sam_itr_next (aln_fh, iter, b)) >= 0) {
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
            if (!(plp_pos >= reg.start && plp_pos < reg.end)) {
                continue;
            }
            aev.count_pileup (pl, safe_size (n_plp));
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
        if (!(plp_pos >= reg.start && plp_pos < reg.end)) {
            continue;
        }
        aev.count_pileup (pl, safe_size (n_plp));
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


int main (int argc,
          char *argv[]) {
    namespace fs = std::filesystem;

    fs::path aln_path;
    std::string region_str;
    htsFile *aln_in;
    bam_hdr_t *head;
    hts_region reg;
    count_params cp;

    // defaults
    cp.min_mapq = 25;
    cp.min_baseq = 30;
    cp.exclude_flag = 3844;
    cp.max_depth = 1000000;
    cp.clip_bound = 0;
    // int keep_flag = 0;
    // int max_mismatch = 0; // ???

    try {
        cxxopts::Options options (
            "count-alleles",
            "c++ implementation of bam2R\n\n"
            "Where reference names contain colons, surround in curly "
            "braces like {HLA-DRB1*12:17}:<start>-<end>\n\n"

            "chr1:100 is treated as the single base pair region "
            "chr1:100-100.\n"
            "chr1:-100 is shorthand for chr1:1-100 and chr1:100- is "
            "ch1:100-<end>\n.");

        // clang-format off
        options.add_options()
            ("aln", "", cxxopts::value<fs::path>())  // positional
            ("region", "", cxxopts::value<std::string>())

            // parameters
            ("b,baseq",
             "Minimum base quality to treat base as unambiguous. (default 30)",
             cxxopts::value<int>())
            ("m,mapq",
             "Minimum mapping quality to include read (default 25)",
             cxxopts::value<int>())
            ("c,clip",
             "Treat bases within <clip> bases of read edges as ambiguous. (default 0)",
             cxxopts::value<int>())
            ("e,exclude",
             "Exclude reads with any bits set in sam flag. Provide flag as integer. (default 3844)",
             cxxopts::value<int>())
            ("d,depth",
             "Maximum read depth (default 1000000)",
             cxxopts::value<int>())

            ("h,help", "Print usage");
        // clang-format on

        options.parse_positional ({"aln", "region"});
        options.positional_help ("<.BAM/.CRAM> chr:start-end");
        auto parsed_args = options.parse (argc, argv);

        if ((!parsed_args.count ("aln")) ||
            (!parsed_args.count ("region"))) {
            std::cout << "incorrect usage: all postional arguments "
                         "required. Try --help"
                      << std::endl;
            return 1;
        }

        if (parsed_args.count ("help")) {
            std::cout << options.help() << std::endl;
            return 0; // nothing given nothing done
        }

        aln_path = parsed_args["aln"].as<fs::path>();
        region_str = parsed_args["region"].as<std::string>();

        if (region_str.empty())
            throw std::runtime_error (
                "region string appears to be empty");

        if (parsed_args.count ("baseq")) {
            cp.min_baseq = parsed_args["baseq"].as<int>();
        }
        if (parsed_args.count ("mapq")) {
            cp.min_mapq = parsed_args["mapq"].as<int>();
        }
        if (parsed_args.count ("clip")) {
            cp.clip_bound = parsed_args["clip"].as<int>();
        }
        if (parsed_args.count ("exclude")) {
            cp.exclude_flag = parsed_args["exclude"].as<int>();
        }
        if (parsed_args.count ("depth")) {
            cp.max_depth = parsed_args["depth"].as<int>();
        }

    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what()
                  << std::endl;
        return 1;
    }

    // NOTE/BUG: there's a very good chance I introduced an off by
    // one, check carefully
    hts_idx_t *idx;
    int tid = -3;
    int64_t start, end;
    std::vector<int> result;
    try {
        aln_in = hts_open (aln_path.c_str(), "r");
        head = sam_hdr_read (aln_in);
        if (head == NULL) {
            throw std::runtime_error (
                "failed to get header from alignment file");
        }

        printf ("%s\n", region_str.c_str());
        auto rp =
            sam_parse_region (head, region_str.c_str(), &tid, &start,
                              &end, HTS_PARSE_ONE_COORD);
        if (rp == NULL) {
            std::string msg;
            switch (tid) {
                case -2:
                    msg = "memory error";
                    break;
                case -1:
                    msg = "could not parse contig";
                    break;
                default:
                    msg = "specified range could not be parsed";
            }
            throw std::runtime_error (
                "parse failed for input region " + region_str +
                " - " + msg);
        }
        // start - 1 cargo culted from bam2R...
        reg = hts_region::by_end (tid, start - 1, end);

        idx = sam_index_load (aln_in, aln_path.c_str());
        if (!idx) {
            throw std::runtime_error ("failed to load index file");
        }

        safe_size_opts sso;
        sso.msg =
            "error in calculating cells needed for storing result";
        result.resize (safe_size (static_cast<int64_t> (
                                      reg.rlen * N_FIELDS_PER_OBS),
                                  sso),
                       0);

    } catch (std::exception &e) {
        std::cerr << "Error during setup: " << e.what() << std::endl;
        return 1;
    }

    try {
        count (aln_in, idx, reg, cp, result);
    } catch (std::exception &e) {
        std::cerr << "Error during calculation: " << e.what()
                  << std::endl;
        return 1;
    }

    try {
        for (size_t i = 0; i < result.size(); i += N_FIELDS_PER_OBS) {
            size_t j = 0;
            while (j < (N_FIELDS_PER_OBS - 1)) {
                std::cout << result[i + j] << ",";
                ++j;
            }
            std::cout << result[i + j + 1] << "\n";
        }
    } catch (std::exception &e) {
        std::cerr << "Error during write: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
