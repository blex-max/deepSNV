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

constexpr std::string_view HEADER =
    "A,T,C,G,-,N,INS,DEL,HEAD,TAIL,QUAL,a,t,c,g,_,n,ins,del,head,"
    "tail,qual";

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


// nothing but C please
struct pf_capture {
    htsFile *fh;
    hts_itr_t *it;
    const count_params *p;
};
int pileup_func (void *data,
                 bam1_t *b) {
    pf_capture *d = static_cast<pf_capture *> (data);
    int ret;
    // find the next good read
    while (1) {
        ret = sam_itr_next (d->fh, d->it, b);
        if (ret < 0) {
            break; // EOF/err
        }
        if (!(b->core.flag & d->p->exclude_flag) &&
            b->core.qual >= d->p->min_mapq) {
            break; // found good read
        };
    }
    return ret;
};
// end nothing but C

// bam2R
inline void count (htsFile *aln_fh,
                   hts_idx_t *aln_idx,
                   const hts_region reg,
                   const count_params params,
                   std::vector<int> &counts
                   // int keep_flag,
                   // int maxmismatches
) {
    bam_plp_t buf = NULL;
    bam1_t *b = NULL;
    bam_hdr_t *head = NULL;
    AlleleEventCounter aev (params, counts);

    std::cerr << "reg start: " << std::to_string (reg.start)
              << std::endl;
    std::cerr << "reg end: " << std::to_string (reg.end) << std::endl;
    safe_size_opts sso_plp_pos;
    sso_plp_pos.msg = "error translating htslib pileup position into "
                      "appropriate index for results array";

    // int64_t maxNM = (maxmismatches != -1) ? maxmismatches :
    // INT64_MAX; unsigned long long no_NM_count = 0;

    // fetch a read overlapping the query region;
    // then do a pileup per base for the total region
    // covered by the retrieved read;
    // then count events on those pileups which overlap
    // the original query region.
    hts_itr_t *iter =
        sam_itr_queryi (aln_idx, reg.rid, reg.start, reg.end);

    pf_capture pfc{aln_fh, iter, &params};
    buf = bam_plp_init (pileup_func,
                        &pfc); // initialize pileup
    bam_plp_set_maxcnt (buf, params.max_depth);

    int64_t plp_pos = -1;
    int plp_tid = -1, n_plp = -1;
    const bam_pileup1_t *pl;
    size_t pos_offset;
    while ((pl = bam_plp64_auto (buf, &plp_tid, &plp_pos, &n_plp)) !=
           0) {
        std::cerr << "plp_pos: " << std::to_string (plp_pos)
                  << std::endl;
        if (n_plp < 0 || plp_tid < 0 || plp_pos < 0) {
            throw std::runtime_error ("pileup failed");
        }
        if (!(plp_pos >= reg.start && plp_pos < reg.end)) {
            std::cerr << "skipping" << std::endl;
            continue;
        }
        pos_offset = safe_size (plp_pos - reg.start, sso_plp_pos);
        std::cerr << "counting plp_pos: " << std::to_string (plp_pos)
                  << " to offset start: "
                  << std::to_string (pos_offset) << std::endl;
        std::cerr << "n_plp: " << std::to_string (n_plp) << std::endl;
        aev.count_pileup (pl, pos_offset, safe_size (n_plp));
    }

    // if (maxmismatches != -1 && no_NM_count > 0) {
    //     printf ("%llu reads did not have NM tags; max.mismatches "
    //             "filter was not "
    //             "applied to them.\n",
    //             no_NM_count);
    // }

    sam_itr_destroy (iter);
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
    bool print_head = false;

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
            "standalone implementation of bam2R"
            "\n\n"
            "Where reference names contain colons, surround in"
            "\n"
            "curly braces like {HLA-DRB1*12:17}:<start>-<end>"
            "\n\n"

            "chr1:100 is treated as the single base pair region"
            "\n"
            "chr1:100-100. chr1:-100 is shorthand for chr1:1-100"
            "\n"
            "and chr1:100- is ch1:100-<end>."
            "\n\n"

            "A result matrix with end-start rows and 22 columns"
            "\n"
            "of event counters (see --head) is printed to stdout"
            "\n"
            "as a csv. Logging is printed to stderr."
            "\n");

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

            ("head", "Print header")
            ("h,help", "Print usage");
        // clang-format on

        options.parse_positional ({"aln", "region"});
        options.positional_help ("<.BAM/.CRAM> chr:start-end");
        auto parsed_args = options.parse (argc, argv);

        if (parsed_args.count ("help")) {
            std::cout << options.help() << std::endl;
            return 0; // nothing given nothing done
        }

        if ((!parsed_args.count ("aln")) ||
            (!parsed_args.count ("region"))) {
            std::cout << "incorrect usage: all postional arguments "
                         "required. Try --help"
                      << std::endl;
            return 1;
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
        if (parsed_args.count ("head")) {
            print_head = true;
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
        // bam2R did start-1 and I don't know why
        reg = hts_region::by_end (tid, start, end);

        idx = sam_index_load (aln_in, aln_path.c_str());
        if (!idx) {
            throw std::runtime_error ("failed to load index file");
        }

        safe_size_opts sso;
        sso.msg =
            "error in calculating cells needed for storing result";
        size_t n_cells = safe_size (
            static_cast<int64_t> (reg.rlen * N_FIELDS_PER_OBS));
        std::cerr << "region length: " << std::to_string (reg.rlen)
                  << std::endl;
        std::cerr << "n_cells: " << std::to_string (n_cells)
                  << std::endl;
        result.resize (n_cells, 0);

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
        if (print_head)
            std::cout << HEADER << "\n";
        size_t i = 0;
        while (i < result.size()) {
            size_t j = 0;
            while (j < (N_FIELDS_PER_OBS - 1)) {
                std::cout << result[i + j] << ",";
                ++j;
            }
            std::cout << result[i + j] << "\n";
            i += N_FIELDS_PER_OBS;
        }
    } catch (std::exception &e) {
        std::cerr << "Error during write: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
