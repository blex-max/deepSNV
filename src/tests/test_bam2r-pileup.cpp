#include <catch2/catch_test_macros.hpp>
#include "../bam2r-pileup.hpp"
#include "catch2/catch_message.hpp"
#include "htslib/khash.h"


const NTParams p1 {
  0,
  1,  // len 1 counts array, simple indexing
  10,
  5
};
// allocate counts via vector, but use pointer for the function via .data()
// 1 pos, 2 strands of data, 11 fields per strand - 1 row of a 3D arr, so 2D arr
std::vector<int> counts(N_COUNTS_FIELD*2, 0);


// NOTE: also test head_clip, and different base call
TEST_CASE("score_pile() first member good", "[Overlap Behaviour Test]") {
  khash_t(strh)* kh = kh_init(strh);
  auto test_counts = counts;

  // get qname into overlap table, data into counts table
  static constexpr PileupRead r1{
    6,  // above threshold
    "r1",
    100,
    10,
    1,  // "A"
    11,  // above threshold
    0,
    false, false, false, false
  };
  int ret = score_pile(r1, test_counts.data(), p1, kh);
  REQUIRE( ret == 0 );  // check success
  CAPTURE( test_counts );
  REQUIRE( test_counts[0] == 1 );  // base recorded as A in counts table
  REQUIRE( test_counts[10] == 10 );  // mapq correct

  auto test_snapshot = test_counts;
  assert(test_counts == test_snapshot);  // sanity

  // since qname is already in overlap table with the same base
  // counts shouldn't change as it should be rejected for having the same qname
  auto r2 = r1;
  r2.base_q = 9;  // below threshold
  score_pile(r2, counts.data(), p1, kh);
  REQUIRE( test_counts == test_snapshot );

  kh_destroy(strh, kh);
}


TEST_CASE("score_pile() first member bad qual", "[Overlap Behaviour Test]") {
  khash_t(strh)* kh = kh_init(strh);
  auto test_counts = counts;

  // get qname into overlap table, data into counts table
  static constexpr PileupRead r1{
    6,  // above threshold
    "r1",
    100,
    10,
    1,  // "A"
    9,  // below threshold
    0,
    false, false, false, false
  };
  int ret = score_pile(r1, test_counts.data(), p1, kh);
  REQUIRE( ret == 0 );  // check success
  CAPTURE( test_counts );
  REQUIRE( test_counts[5] == 1 );  // base recorded as N in counts table
  REQUIRE( test_counts[10] == 10 );  // mapq correct

  auto test_snapshot = test_counts;
  assert(test_counts == test_snapshot);  // sanity

  // counts should update despite qname presence
  // for better 2nd member of read pair
  auto r2 = r1;
  r2.base_q = 11;  // below threshold
  score_pile(r2, counts.data(), p1, kh);
  REQUIRE( test_counts != test_snapshot );  // counts should not match snapshot (fails, indicating BUG)

  kh_destroy(strh, kh);
}
