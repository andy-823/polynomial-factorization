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
  return Element(result);
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

enum GenMode { kFixed, kRandom };
// template because different polynomials may want different sizes
template <typename Poly, size_t kMaxSize, GenMode kMode, typename RandomGen>
size_t GenSize(RandomGen& gen) {
  if constexpr (kMode == kFixed) {
    return kMaxSize;
  } else {
    return gen() % kMaxSize;
  }
}

template <concepts::Polynom Poly, size_t kMaxSize = 128,
          GenMode kMode = kRandom, typename RandomGen>
Poly GenPoly(RandomGen& gen) {
  using Element = typename Poly::Element;

  do {
    std::vector<Element> elements(GenSize<Poly, kMaxSize, kMode>(gen));
    for (auto& element : elements) {
      element = GenElement<Element>(gen);
    }
    Poly result(std::move(elements));
    if (!result.IsZero()) {
      return result;
    }
  } while (true);
}

}  // namespace factorization
