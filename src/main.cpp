#include "bounds.hpp"
#include <cstddef>
#include <cstdint>
#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "bam2R.hpp"
#include "structs.hpp"

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
            (!parsed_args.count ("region")) ||
            parsed_args.count ("help")) {
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
    int tid = -3;
    int64_t start, end;
    try {
        aln_in = hts_open (aln_path.c_str(), "r");
        head = sam_hdr_read (aln_in);
        if (head == NULL) {
            throw std::runtime_error (
                "failed to get header from alignment file");
        }

        printf ("%s\n", region_str.c_str());
        auto rp = sam_parse_region (head, region_str.c_str(), &tid,
                                    &start, &end, HTS_PARSE_ONE_COORD);
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
        reg = hts_region::by_end(tid, start - 1, end);

    } catch (std::exception &e) {
        std::cerr << "Error during setup: " << e.what() << std::endl;
        return 1;
    }

    std::pair<size_t, int *> result;
    try {
    } catch (std::exception &e) {
        std::cerr << "Error during calculation: " << e.what()
                  << std::endl;
        return 1;
    }

    // NOTE/BUG:
    // Given a 1D vector,
    // R translates data to a matrix in column major style
    // i.e. it writes top to bottom in column 0,
    // then fills column 1]
    // Hence this data is column major.
    // At present, the output is just an unformatted
    // stream of comma separated values
    // translate into matrix per R
    // then write out that matrix line by line
    try {
        for (size_t i = 0; i < (result.first - 1); ++i)
            std::cout << result.second[i] << ",";
        // flush last result without the comma
        std::cout << result.second[result.first - 1];
    } catch (std::exception &e) {
        std::cerr << "Error during write: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
