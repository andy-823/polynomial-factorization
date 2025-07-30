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
    {QueryType::kMultiply, 0, 4, 0},
    {QueryType::kMultiply, 0, 5, 0},
    {QueryType::kMultiply, 0, 6, 0},
    {QueryType::kMultiply, 0, 8, 0},
    {QueryType::kMultiply, 0, 9, 0},
    {QueryType::kMultiply, 0, 10, 0},

    {QueryType::kMultiply, 1, 0, 0},
    {QueryType::kMultiply, 1, 1, 1},
    {QueryType::kMultiply, 1, 2, 2},
    {QueryType::kMultiply, 1, 8, 8},
    {QueryType::kMultiply, 1, 9, 9},
    {QueryType::kMultiply, 1, 10, 10},
    {QueryType::kMultiply, 1, 16, 16},
    {QueryType::kMultiply, 1, 17, 17},
    {QueryType::kMultiply, 1, 18, 18},

    {QueryType::kMultiply, 2, 0, 0},
    {QueryType::kMultiply, 2, 1, 2},
    {QueryType::kMultiply, 2, 2, 1},
    {QueryType::kMultiply, 2, 8, 16},
    {QueryType::kMultiply, 2, 9, 18},
    {QueryType::kMultiply, 2, 10, 17},
    {QueryType::kMultiply, 2, 16, 8},
    {QueryType::kMultiply, 2, 17, 10},
    {QueryType::kMultiply, 2, 18, 9},

    {QueryType::kMultiply, 8, 0, 0},
    {QueryType::kMultiply, 8, 1, 8},
    {QueryType::kMultiply, 8, 2, 16},
    // x * x = x + 1
    {QueryType::kMultiply, 8, 8, 9},
    // x * (x + 1) = x^2 + x = 2x + 1
    {QueryType::kMultiply, 8, 9, 17},
    // x * (x + 2) = x^2 + 2x = 1
    {QueryType::kMultiply, 8, 10, 1},
    // x * (2x) = 2x^2 = 2x + 2
    {QueryType::kMultiply, 8, 16, 18},
    // x * (2x + 1) = 2x + 2 + x = 2
    {QueryType::kMultiply, 8, 17, 2},
    // x * (2x + 2) = 2x + 2 + 2x = x + 2
    {QueryType::kMultiply, 8, 18, 10},

    {QueryType::kAdd, 2, 0, 2},
    {QueryType::kAdd, 8, 1, 9},
    {QueryType::kAdd, 10, 2, 9},
    {QueryType::kAdd, 17, 8, 1},
    {QueryType::kAdd, 18, 9, 0},
    {QueryType::kAdd, 9, 10, 16},
    {QueryType::kAdd, 2, 16, 18},
    {QueryType::kAdd, 8, 17, 1},
    {QueryType::kAdd, 1, 18, 16},
  };

  // x^2 = x + 1
  using GaloisField = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
  using Element = galois_field::FieldElementWrapper<GaloisField>;
  using Int = Element::Value;

  STATIC_REQUIRE(Element::FieldBase() == 3);
  STATIC_REQUIRE(Element::FieldPower() == 2);

  std::vector<Int> elements({
      Int{0}, Int{1}, Int{2}, Int{8}, Int{9},
      Int{10}, Int{16}, Int{17}, Int{18}
    });
  CHECK_THAT(Element::AllFieldElements(), RangeEquals(elements));

  RunTests<Element>(tests);
}
