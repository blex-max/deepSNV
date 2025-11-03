#include <cstdint>
#include <format>
#include <limits>
#include <stdexcept>

#include "bounds.hpp"

bool in_bounds (int64_t i,
                int64_t lower,
                int64_t upper) {
    return lower <= i && i < upper;
}

int int_diff (int64_t i,
              int64_t sub) {
    int64_t diff = i - sub;
    if (diff < 0 || diff > std::numeric_limits<int>::max())
        throw std::out_of_range ("length out of int range");
    return static_cast<int> (diff);
}


size_t safe_size (int64_t i,
                  safe_size_opts opts) {
    try {
        if (i < 0)
            throw std::out_of_range (
                std::format ("size {} would be negative", i));

        // Cast fine since we know non negative now
        if (static_cast<uint64_t> (i) < opts.lower)
            throw std::out_of_range (std::format (
                "size {} would be below lower bound {}", i, opts.lower));

        if (static_cast<uint64_t> (i) > opts.upper)
            throw std::out_of_range (std::format (
                "size {} would exceed upper bound {}", i, opts.upper));

        // convert in safety
        return static_cast<size_t> (i);
    } catch (std::exception &e) {
        throw std::runtime_error (
            std::format ("{}\n{}", opts.msg, e.what()));
    }
}
