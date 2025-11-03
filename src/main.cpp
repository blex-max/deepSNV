#include <cstddef>
#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <memory>
#include <string>

#include "bam2R.hpp"

struct hts_region {
    std::string contig;
    int64_t start;
    int64_t stop;
};

int main (int argc,
          char *argv[]) {
    namespace fs = std::filesystem;

    fs::path aln_path;
    std::string region_str;

    // defaults
    int mq = 25;
    int bq = 30;
    int exclude_flag = 3844;
    int max_depth = 100000;
    int head_clip = 0;
    int keepflag = 0;
    int max_mismatch = 0; // ???

    try {
        // clang-format off
        cxxopts::Options options ("count-alleles" "cxx only implementation of bam2R");

        options.add_options()
            ("aln", "", cxxopts::value<fs::path>())  // positional

            ("r,region", "region specification as string chr:start-stop", cxxopts::value<std::string>())

            // parameters
            ("m,mapq",
             "Minimum mapping quality",
             cxxopts::value<double>())

            ("h,help", "Print usage");
        // clang-format on

        options.parse_positional ({"aln", "region"});
        options.positional_help ("<ALN> region-str");
        auto result = options.parse (argc, argv);

        if ((!result.count ("aln")) || (!result.count ("region")) ||
            result.count ("help")) {
            std::cout << options.help() << std::endl;
            return 0; // nothing given nothing done
        }

        aln_path = result["aln"].as<fs::path>();
        region_str = result["region"].as<std::string>();

        if (region_str.empty())
            throw std::runtime_error (
                "region string appears to be empty");

    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what()
                  << std::endl;
        return 1;
    }

    htsFile *aln_in;
    bam_hdr_t *head;
    int tid;
    int64_t beg;
    int64_t end;
    try {
        aln_in = hts_open (aln_path.c_str(), "r");
        head = sam_hdr_read (aln_in);
        if (head == NULL) {
            throw std::runtime_error (
                "failed to get header from alignment file");
        }

        printf ("%s\n", region_str.c_str());
        auto rp = sam_parse_region (head, region_str.c_str(), &tid,
                                    &beg, &end, 0);
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
        // std::cout << std::to_string(tid) + " " +
        // std::to_string(beg) + " " + std::to_string(end) <<
        // std::endl;


    } catch (std::exception &e) {
        std::cerr << "Error during setup: " << e.what() << std::endl;
        return 1;
    }

    try {
        auto result =
            bam2R (aln_in, aln_path, tid, beg, end, bq, mq, head_clip,
                   max_depth, exclude_flag, keepflag, max_mismatch);
    } catch (std::exception &e) {
        std::cerr << "Error during calculation: " << e.what()
                  << std::endl;
        return 1;
    }

    // TODO: write array

    return 0;
}
