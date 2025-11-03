#pragma once

#include "bounds.hpp"
#include <cstdint>
#include <format>
#include <htslib/hts.h>
#include <limits>
#include <stdexcept>


struct hts_region {
    int32_t rid;
    int64_t gstart;
    int64_t gend;
    size_t rlen;

    // default invalid constructor
    hts_region () noexcept
        : rid (-1),
          gstart (-1),
          gend (-1),
          rlen (0) {}

    static hts_region by_end (int32_t rid_,
                              int64_t gstart_,
                              int64_t gend_) {
        hts_region r;
        r.rid = rid_;
        r.gstart = gstart_;
        r.gend = gend_;
        if (!r.valid_rid() || !r.valid_span())
            throw std::invalid_argument (std::format (
                "hts_region::by_end invalid parameters\nrid "
                "{}\ngstart {}\ngend {}",
                r.rid, r.gstart, r.gend));
        r.rlen = safe_size (
            r.gend - r.gstart,
            {.msg = std::format (
                 "hts_region::by_end - span too large {} {}",
                 r.gstart, r.gend)});

        return r;
    }

    static hts_region by_len (int32_t rid_,
                              int64_t gstart_,
                              size_t rlen_) {
        hts_region r;
        r.rid = rid_;
        r.gstart = gstart_;
        r.rlen = rlen_;
        if (!r.valid_rid() || !r.valid_rlen())
            throw std::invalid_argument (
                std::format ("hts_region: invalid parameters\nrid "
                             "{}\n gstart {}\nrlen {}",
                             r.rid, r.gstart, r.rlen));
        r.gend = r.gstart + static_cast<int64_t> (r.rlen);
        return r;
    }

    bool valid_rid () const noexcept { return rid >= 0; }
    bool valid_span () const noexcept { return gend > gstart; }
    bool valid_rlen () const noexcept {
        return (std::numeric_limits<int64_t>::max() -
                    static_cast<int64_t> (rlen) >
                gstart);
    }
};
