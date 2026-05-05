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

template <std::integral Int, uint32_t kPower>
struct Test {
  QueryType type;
  std::array<Int, kPower> first;
  std::array<Int, kPower> second;
  std::array<Int, kPower> expected;
};

template <concepts::GaloisFieldElement Element, typename Value, uint32_t kPower>
void RunTests(const std::vector<Test<Value, kPower>>& tests) {
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

        // case QueryType::kPow:
        //   REQUIRE(first.Pow(test.second) == expected);
        //   break;

      default:
        throw std::runtime_error("Unexpected query type");
    }
  }
}

TEST_CASE("FieldElementWrapper") {
  std::vector<Test<int64_t, 2>> tests = {
      {QueryType::kMultiply, {0, 0}, {0, 0}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {1, 0}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {2, 0}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {0, 1}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {1, 1}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {2, 1}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {0, 2}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {1, 2}, {0, 0}},
      {QueryType::kMultiply, {0, 0}, {2, 2}, {0, 0}},

      {QueryType::kMultiply, {1, 0}, {0, 0}, {0, 0}},
      {QueryType::kMultiply, {1, 0}, {1, 0}, {1, 0}},
      {QueryType::kMultiply, {1, 0}, {2, 0}, {2, 0}},
      {QueryType::kMultiply, {1, 0}, {0, 1}, {0, 1}},
      {QueryType::kMultiply, {1, 0}, {1, 1}, {1, 1}},
      {QueryType::kMultiply, {1, 0}, {2, 1}, {2, 1}},
      {QueryType::kMultiply, {1, 0}, {0, 2}, {0, 2}},
      {QueryType::kMultiply, {1, 0}, {1, 2}, {1, 2}},
      {QueryType::kMultiply, {1, 0}, {2, 2}, {2, 2}},

      {QueryType::kMultiply, {2, 0}, {0, 0}, {0, 0}},
      {QueryType::kMultiply, {2, 0}, {1, 0}, {2, 0}},
      {QueryType::kMultiply, {2, 0}, {2, 0}, {1, 0}},
      {QueryType::kMultiply, {2, 0}, {0, 1}, {0, 2}},
      {QueryType::kMultiply, {2, 0}, {1, 1}, {2, 2}},
      {QueryType::kMultiply, {2, 0}, {2, 1}, {1, 2}},
      {QueryType::kMultiply, {2, 0}, {0, 2}, {0, 1}},
      {QueryType::kMultiply, {2, 0}, {1, 2}, {2, 1}},
      {QueryType::kMultiply, {2, 0}, {2, 2}, {1, 1}},

      {QueryType::kMultiply, {0, 1}, {0, 0}, {0, 0}},
      {QueryType::kMultiply, {0, 1}, {1, 0}, {0, 1}},
      {QueryType::kMultiply, {0, 1}, {2, 0}, {0, 2}},
      // x * x = x + 1
      {QueryType::kMultiply, {0, 1}, {0, 1}, {1, 1}},
      // x * (x + 1) = x^2 + x = 2x + 1
      {QueryType::kMultiply, {0, 1}, {1, 1}, {1, 2}},
      // x * (x + 2) = x^2 + 2x = 1
      {QueryType::kMultiply, {0, 1}, {2, 1}, {1, 0}},
      // x * (2x) = 2x^2 = 2x + 2
      {QueryType::kMultiply, {0, 1}, {0, 2}, {2, 2}},
      // x * (2x + 1) = 2x + 2 + x = 2
      {QueryType::kMultiply, {0, 1}, {1, 2}, {2, 0}},
      // x * (2x + 2) = 2x + 2 + 2x = x + 2
      {QueryType::kMultiply, {0, 1}, {2, 2}, {2, 1}},

      {QueryType::kAdd, {2, 0}, {0, 0}, {2, 0}},
      {QueryType::kAdd, {0, 1}, {1, 0}, {1, 1}},
      {QueryType::kAdd, {2, 1}, {2, 0}, {1, 1}},
      {QueryType::kAdd, {1, 2}, {0, 1}, {1, 0}},
      {QueryType::kAdd, {2, 2}, {1, 1}, {0, 0}},
      {QueryType::kAdd, {0, 1}, {2, 1}, {2, 2}},
      {QueryType::kAdd, {2, 0}, {0, 2}, {2, 2}},
      {QueryType::kAdd, {0, 1}, {1, 2}, {1, 0}},
      {QueryType::kAdd, {1, 0}, {2, 2}, {0, 2}},
  };

  // x^2 = x + 1
  using GaloisField = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
  using Element = galois_field::FieldElementWrapper<GaloisField>;
  using Int = Element::Coefficient;

  STATIC_REQUIRE(Element::FieldBase() == 3);
  STATIC_REQUIRE(Element::FieldPower() == 2);

  std::vector<Element> elements({
      Element({0, 0}),
      Element({1, 0}),
      Element({2, 0}),
      Element({0, 1}),
      Element({1, 1}),
      Element({2, 1}),
      Element({0, 2}),
      Element({1, 2}),
      Element({2, 2}),
  });
  CHECK_THAT(Element::AllFieldElements(), RangeEquals(elements));

  RunTests<Element>(tests);
}
