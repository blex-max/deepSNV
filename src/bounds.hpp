#pragma once

#include <cstdint>
#include <limits>
#include <string>

bool in_bounds (int64_t i,
                int64_t lower = 0,
                int64_t upper = std::numeric_limits<int64_t>::max());

int int_diff (int64_t start,
              int64_t end);


struct safe_size_opts {
    size_t lower = 0;
    size_t upper = std::numeric_limits<size_t>::max();
    std::string msg = "";
};

size_t safe_size (int64_t i,
                  safe_size_opts opts=safe_size_opts{});
