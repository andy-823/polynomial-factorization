#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <factorization/concepts.hpp>
#include <factorization/galois_field/log_based_field.hpp>

using namespace factorization::internal;  // NOLINT

enum QueryType { kAdd, kNegative, kMultiply, kInverse, kPow };

template <std::integral Int = uint32_t>
struct Test {
  QueryType type;
  Int first;
  Int second;
  Int expected;
};

template <typename NewInt, typename OldInt>
std::vector<Test<NewInt>> TransformTests(const std::vector<OldInt>& tests) {
  std::vector<Test<NewInt>> result(tests.size());
  for (size_t i = 0; i < tests.size(); ++i) {
    result[i].type = tests[i].type;
    result[i].first = tests[i].first;
    result[i].second = tests[i].second;
    result[i].expected = tests[i].expected;
  }
  return result;
}

template <GaloisField GaloisField, typename Int>
void RunTests(const std::vector<Test<Int>>& tests) {
  GaloisField field;

  for (const Test<Int>& test : tests) {
    switch (test.type) {
      case QueryType::kAdd:
        REQUIRE(field.Add(test.first, test.second) == test.expected);
        break;
      case QueryType::kNegative:
        REQUIRE(field.Negative(test.first) == test.expected);
        break;
      case QueryType::kMultiply:
        REQUIRE(field.Multiply(test.first, test.second) == test.expected);
        break;
      case QueryType::kInverse:
        REQUIRE(field.Inverse(test.first) == test.expected);
        break;
      case QueryType::kPow:
        REQUIRE(field.Pow(test.first, test.second) == test.expected);
        break;
      default:
        throw std::runtime_error("Unexpected query type");
    }
  }
}

TEST_CASE("LogBaseGaloisField") {
  SECTION("GF8") {
    std::vector<Test<int64_t>> tests = {
      {QueryType::kMultiply, 0, 0, 0},
      {QueryType::kMultiply, 0, 1, 0},
      {QueryType::kMultiply, 0, 2, 0},
      {QueryType::kMultiply, 0, 3, 0},
      {QueryType::kMultiply, 0, 4, 0},
      {QueryType::kMultiply, 0, 5, 0},
      {QueryType::kMultiply, 0, 6, 0},
      {QueryType::kMultiply, 0, 7, 0},

      {QueryType::kMultiply, 1, 0, 0},
      {QueryType::kMultiply, 1, 1, 1},
      {QueryType::kMultiply, 1, 2, 2},
      {QueryType::kMultiply, 1, 3, 3},
      {QueryType::kMultiply, 1, 4, 4},
      {QueryType::kMultiply, 1, 5, 5},
      {QueryType::kMultiply, 1, 6, 6},
      {QueryType::kMultiply, 1, 7, 7},
  
      {QueryType::kMultiply, 2, 0, 0},
      {QueryType::kMultiply, 2, 1, 2},
      {QueryType::kMultiply, 2, 2, 4},
      {QueryType::kMultiply, 2, 3, 6},
      {QueryType::kMultiply, 2, 4, 3},
      {QueryType::kMultiply, 2, 5, 1},
      {QueryType::kMultiply, 2, 6, 7},
      {QueryType::kMultiply, 2, 7, 5},

      {QueryType::kMultiply, 3, 0, 0},
      {QueryType::kMultiply, 3, 1, 3},
      {QueryType::kMultiply, 3, 2, 6},
      {QueryType::kMultiply, 3, 3, 5},
      {QueryType::kMultiply, 3, 4, 7},
      {QueryType::kMultiply, 3, 5, 4},
      {QueryType::kMultiply, 3, 6, 1},
      {QueryType::kMultiply, 3, 7, 2},

      {QueryType::kMultiply, 4, 0, 0},
      {QueryType::kMultiply, 4, 1, 4},
      {QueryType::kMultiply, 4, 2, 3},
      {QueryType::kMultiply, 4, 3, 7},
      {QueryType::kMultiply, 4, 4, 6},
      {QueryType::kMultiply, 4, 5, 2},
      {QueryType::kMultiply, 4, 6, 5},
      {QueryType::kMultiply, 4, 7, 1},

      {QueryType::kMultiply, 5, 0, 0},
      {QueryType::kMultiply, 5, 1, 5},
      {QueryType::kMultiply, 5, 2, 1},
      {QueryType::kMultiply, 5, 3, 4},
      {QueryType::kMultiply, 5, 4, 2},
      {QueryType::kMultiply, 5, 5, 7},
      {QueryType::kMultiply, 5, 6, 3},
      {QueryType::kMultiply, 5, 7, 6},

      {QueryType::kMultiply, 6, 0, 0},
      {QueryType::kMultiply, 6, 1, 6},
      {QueryType::kMultiply, 6, 2, 7},
      {QueryType::kMultiply, 6, 3, 1},
      {QueryType::kMultiply, 6, 4, 5},
      {QueryType::kMultiply, 6, 5, 3},
      {QueryType::kMultiply, 6, 6, 2},
      {QueryType::kMultiply, 6, 7, 4},

      {QueryType::kMultiply, 7, 0, 0},
      {QueryType::kMultiply, 7, 1, 7},
      {QueryType::kMultiply, 7, 2, 5},
      {QueryType::kMultiply, 7, 3, 2},
      {QueryType::kMultiply, 7, 4, 1},
      {QueryType::kMultiply, 7, 5, 6},
      {QueryType::kMultiply, 7, 6, 4},
      {QueryType::kMultiply, 7, 7, 3},

      {QueryType::kNegative, 1, -1, 1},
      {QueryType::kNegative, 2, -1, 2},
      {QueryType::kNegative, 4, -1, 4},
      {QueryType::kNegative, 7, -1, 7},
      {QueryType::kAdd, 0, 2, 2},
      {QueryType::kAdd, 3, 5, 6},
      {QueryType::kAdd, 4, 4, 0},
      {QueryType::kAdd, 1, 6, 7},

      {QueryType::kInverse, 1, -1, 1},
      {QueryType::kInverse, 2, -1, 5},
      {QueryType::kInverse, 3, -1, 6},
      {QueryType::kInverse, 4, -1, 7},
      {QueryType::kInverse, 5, -1, 2},
      {QueryType::kInverse, 6, -1, 3},
      {QueryType::kInverse, 7, -1, 4},

      {QueryType::kPow, 2, 0, 1},
      {QueryType::kPow, 2, 1, 2},
      {QueryType::kPow, 2, 2, 4},
      {QueryType::kPow, 2, 3, 3},
      {QueryType::kPow, 2, 4, 6},
      {QueryType::kPow, 2, 5, 7},
      {QueryType::kPow, 2, 6, 5},
      {QueryType::kPow, 2, 7, 1},
    };

    using GaloisField = LogBasedField<2, 3, {1, 1, 0, 1}>;
    RunTests<GaloisField>(tests);
  }

  SECTION("GF9") {


  }
}
