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

// template <concepts::GaloisFieldElement Element, std::ranges::range Range>
// std::vector<Element> MakeElements(Range&& range) {
//   std::vector<Element> result(range.begin(), range.end());
//   return result;
// }

TEST_CASE("FieldElementWrapper") {
  std::vector<Test<int64_t>> tests = {
    {QueryType::kMultiply, 0, 0, 0},
    {QueryType::kMultiply, 0, 1, 0},
    {QueryType::kMultiply, 0, 2, 0},
    {QueryType::kMultiply, 0, 3, 0},
    {QueryType::kMultiply, 0, 4, 0},
    {QueryType::kMultiply, 0, 5, 0},
    {QueryType::kMultiply, 0, 6, 0},
    {QueryType::kMultiply, 0, 7, 0},
    {QueryType::kMultiply, 0, 8, 0},

    {QueryType::kMultiply, 1, 0, 0},
    {QueryType::kMultiply, 1, 1, 1},
    {QueryType::kMultiply, 1, 2, 2},
    {QueryType::kMultiply, 1, 3, 3},
    {QueryType::kMultiply, 1, 4, 4},
    {QueryType::kMultiply, 1, 5, 5},
    {QueryType::kMultiply, 1, 6, 6},
    {QueryType::kMultiply, 1, 7, 7},
    {QueryType::kMultiply, 1, 8, 8},

    {QueryType::kMultiply, 2, 0, 0},
    {QueryType::kMultiply, 2, 1, 2},
    {QueryType::kMultiply, 2, 2, 1},
    {QueryType::kMultiply, 2, 3, 6},
    {QueryType::kMultiply, 2, 4, 8},
    {QueryType::kMultiply, 2, 5, 7},
    {QueryType::kMultiply, 2, 6, 3},
    {QueryType::kMultiply, 2, 7, 5},
    {QueryType::kMultiply, 2, 8, 4},

    {QueryType::kMultiply, 3, 0, 0},
    {QueryType::kMultiply, 3, 1, 3},
    {QueryType::kMultiply, 3, 2, 6},
    // x * x = x + 1
    {QueryType::kMultiply, 3, 3, 4},
    // x * (x + 1) = x^2 + x = 2x + 1
    {QueryType::kMultiply, 3, 4, 7},
    // x * (x + 2) = x^2 + 2x = 1
    {QueryType::kMultiply, 3, 5, 1},
    // x * (2x) = 2x^2 = 2x + 2
    {QueryType::kMultiply, 3, 6, 8},
    // x * (2x + 1) = 2x + 2 + x = 2
    {QueryType::kMultiply, 3, 7, 2},
    // x * (2x + 2) = 2x + 2 + 2x = x + 2
    {QueryType::kMultiply, 3, 8, 5},

    {QueryType::kAdd, 2, 0, 2},
    {QueryType::kAdd, 3, 1, 4},
    {QueryType::kAdd, 5, 2, 4},
    {QueryType::kAdd, 7, 3, 1},
    {QueryType::kAdd, 8, 4, 0},
    {QueryType::kAdd, 4, 5, 6},
    {QueryType::kAdd, 2, 6, 8},
    {QueryType::kAdd, 3, 7, 1},
    {QueryType::kAdd, 1, 8, 6},
  };

  // x^2 = x + 1
  using GaloisField = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
  using Element = galois_field::FieldElementWrapper<GaloisField>;
  using Int = Element::Value;

  STATIC_REQUIRE(Element::FieldBase() == 3);
  STATIC_REQUIRE(Element::FieldPower() == 2);

  std::vector<Int> elements({
      Int{0}, Int{1}, Int{2}, Int{3}, Int{4},
      Int{5}, Int{6}, Int{7}, Int{8}
    });
  CHECK_THAT(Element::AllFieldElements(), RangeEquals(elements));

  RunTests<Element>(tests);
}
