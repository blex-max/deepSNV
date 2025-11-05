#pragma once

#include <cstddef>
#include <cstdint>

inline constexpr uint8_t HTS_NT_A = 1;
inline constexpr uint8_t HTS_NT_C = 2;
inline constexpr uint8_t HTS_NT_G = 4;
inline constexpr uint8_t HTS_NT_T = 8;

inline constexpr uint8_t COUNT_A = 0;
inline constexpr uint8_t COUNT_T = 1;
inline constexpr uint8_t COUNT_C = 2;
inline constexpr uint8_t COUNT_G = 3;
inline constexpr uint8_t COUNT_IS_DEL = 4; // *
inline constexpr uint8_t COUNT_N = 5; // N
inline constexpr uint8_t COUNT_INS = 6; // +
inline constexpr uint8_t COUNT_DEL = 7; // -
inline constexpr uint8_t COUNT_HEAD = 8; // ^
inline constexpr uint8_t COUNT_TAIL = 9; // $
inline constexpr uint8_t COUNT_MAPQ = 10; // Q

inline constexpr size_t N_FIELDS_PER_STRAND = 11;
inline constexpr size_t N_STRAND = 2;
inline constexpr size_t N_FIELDS_PER_OBS = N_STRAND * N_FIELDS_PER_STRAND;
inline constexpr size_t RSTRAND_OFFSET = N_FIELDS_PER_STRAND;

inline constexpr uint8_t UNDEFINED_VALUE = UINT8_MAX;
