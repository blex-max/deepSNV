#include <catch2/catch_test_macros.hpp>
#include "../bam2r_pileup.hpp"
#include "htslib/khash.h"


static constexpr NTParams p1 {
  0,
  100,
  10,
  0
};

static constexpr PileupRead r1{
  0,
  "r1",
  100,
  10,
  1,  // "A"
  32,
  0,
  false, false, false, false
};


static khash_t(strh)* prefill_khash() {
  int key_in_hash;
  khash_t(strh)* kh_overlap_ptr = kh_init(strh);
  khiter_t khi = kh_put(strh, kh_overlap_ptr, r1.qname, &key_in_hash);

  kh_value(kh_overlap_ptr, khi) = r1.base_nt16i;
  return kh_overlap_ptr;  // caller must destroy
}


// allocate counts via vector, but use pointer for the function via .data()
// 1 pos, 2 strands of data, 11 fields per strand - 1 row of a 3D arr, so 2D arr
std::vector<int> counts(N_COUNTS_FIELD*2, 0);
// snapshot before test
std::vector<int> counts_snapshot = counts;

TEST_CASE("score pile overlap bug confirm" "[idk what these strings do]") {
  auto kh = prefill_khash();

  score_pile(r1, counts.data(), p1, kh);

  // since we're processing the same struct twice, counts shouldn't change
  // as it should be rejected for having the same qname
  REQUIRE( counts == counts_snapshot );

  kh_destroy(strh, kh);
}

