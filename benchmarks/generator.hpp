#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

namespace factorization {

template <concepts::GaloisFieldElement Element>
constexpr auto AllElements() {
  constexpr int kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());
  auto elements_range = Element::AllFieldElements();
  std::array<Element, kFieldSize> result{};
  std::copy(elements_range.begin(), elements_range.end(), result.begin());
  return result;
}

template <concepts::GaloisFieldElement Element, typename RandomGen>
Element GenElement(RandomGen& gen) {
  using T = typename Element::Coefficient;
  constexpr auto kFieldBase = Element::FieldBase();
  constexpr auto kFieldPower = Element::FieldPower();

  std::array<T, kFieldPower> result;
  for (auto& c : result) {
    c = gen() % kFieldBase;
  }
  return Element(result);
}

template <concepts::Polynom Poly, typename RandomGen>
Poly GenPoly(RandomGen& gen, size_t size) {
  using Element = typename Poly::Element;

  std::vector<Element> elements(size + 1);
  for (auto& element : elements) {
    element = GenElement<Element>(gen);
  }
  elements.back() = Element::One();
  return Poly(std::move(elements));
}

}  // namespace factorization
