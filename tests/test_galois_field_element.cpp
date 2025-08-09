#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <factorization/concepts.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>

using namespace factorization;       // NOLINT
using Catch::Matchers::RangeEquals;  // NOLINT

enum QueryType { kAdd, kNegative, kMultiply, kInverse, kPow };

template <std::integral Int = uint32_t>
struct Test {
  QueryType type;
  Int first;
  Int second;
  Int expected;
};

template <concepts::GaloisFieldElement Element, typename Value>
void RunTests(const std::vector<Test<Value>>& tests) {
  for (const auto& test : tests) {
    Element first(test.first);
    Element second;
    Element expected(test.expected);

    switch (test.type) {
      case QueryType::kAdd:
        second = Element(test.second);

        REQUIRE((Element(first) += second) == expected);
        REQUIRE(first + second == expected);
        REQUIRE(first + Element::Zero() == first);

        REQUIRE((Element(expected) -= first) == second);
        REQUIRE(expected - first == second);
        REQUIRE(first - Element::Zero() == first);
        break;

      case QueryType::kNegative:
        REQUIRE(-first == expected);
        break;

      case QueryType::kMultiply:
        second = Element(test.second);

        REQUIRE((Element(first) *= second) == expected);
        REQUIRE(first * second == expected);
        REQUIRE(first * Element::One() == first);

        if (expected != Element::Zero()) {
          REQUIRE((Element(expected) /= first) == second);
          REQUIRE(expected / first == second);
        }
        break;

      case QueryType::kInverse:
        REQUIRE(first.Inverse() == expected);
        break;

      case QueryType::kPow:
        REQUIRE(first.Pow(test.second) == expected);
        break;

      default:
        throw std::runtime_error("Unexpected query type");
    }
  }
}

TEST_CASE("FieldElementWrapper") {
  std::vector<Test<int64_t>> tests = {
      {QueryType::kMultiply, 0, 0, 0},
      {QueryType::kMultiply, 0, 1, 0},
      {QueryType::kMultiply, 0, 2, 0},
      {QueryType::kMultiply, 0, 7, 0},
      {QueryType::kMultiply, 0, 7, 0},
      {QueryType::kMultiply, 0, 8, 0},
      {QueryType::kMultiply, 0, 12, 0},
      {QueryType::kMultiply, 0, 13, 0},
      {QueryType::kMultiply, 0, 14, 0},

      {QueryType::kMultiply, 1, 0, 0},
      {QueryType::kMultiply, 1, 1, 1},
      {QueryType::kMultiply, 1, 2, 2},
      {QueryType::kMultiply, 1, 6, 6},
      {QueryType::kMultiply, 1, 7, 7},
      {QueryType::kMultiply, 1, 8, 8},
      {QueryType::kMultiply, 1, 12, 12},
      {QueryType::kMultiply, 1, 13, 13},
      {QueryType::kMultiply, 1, 14, 14},

      {QueryType::kMultiply, 2, 0, 0},
      {QueryType::kMultiply, 2, 1, 2},
      {QueryType::kMultiply, 2, 2, 1},
      {QueryType::kMultiply, 2, 6, 12},
      {QueryType::kMultiply, 2, 7, 14},
      {QueryType::kMultiply, 2, 8, 13},
      {QueryType::kMultiply, 2, 12, 6},
      {QueryType::kMultiply, 2, 13, 8},
      {QueryType::kMultiply, 2, 14, 7},

      {QueryType::kMultiply, 6, 0, 0},
      {QueryType::kMultiply, 6, 1, 6},
      {QueryType::kMultiply, 6, 2, 12},
      // x * x = x + 1
      {QueryType::kMultiply, 6, 6, 7},
      // x * (x + 1) = x^2 + x = 2x + 1
      {QueryType::kMultiply, 6, 7, 13},
      // x * (x + 2) = x^2 + 2x = 1
      {QueryType::kMultiply, 6, 8, 1},
      // x * (2x) = 2x^2 = 2x + 2
      {QueryType::kMultiply, 6, 12, 14},
      // x * (2x + 1) = 2x + 2 + x = 2
      {QueryType::kMultiply, 6, 13, 2},
      // x * (2x + 2) = 2x + 2 + 2x = x + 2
      {QueryType::kMultiply, 6, 14, 8},

      {QueryType::kAdd, 2, 0, 2},
      {QueryType::kAdd, 6, 1, 7},
      {QueryType::kAdd, 8, 2, 7},
      {QueryType::kAdd, 13, 6, 1},
      {QueryType::kAdd, 14, 7, 0},
      {QueryType::kAdd, 6, 8, 14},
      {QueryType::kAdd, 2, 12, 14},
      {QueryType::kAdd, 6, 13, 1},
      {QueryType::kAdd, 1, 14, 12},
  };

  // x^2 = x + 1
  using GaloisField = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
  using Element = galois_field::FieldElementWrapper<GaloisField>;
  using Int = Element::Value;

  STATIC_REQUIRE(Element::FieldBase() == 3);
  STATIC_REQUIRE(Element::FieldPower() == 2);

  std::vector<Int> elements({
      Int{0}, Int{1}, Int{2}, Int{6}, Int{7},
      Int{8}, Int{12}, Int{13}, Int{14}
    });
  CHECK_THAT(Element::AllFieldElements(), RangeEquals(elements));

  RunTests<Element>(tests);
}
