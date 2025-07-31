#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

namespace factorization {

template <concepts::GaloisFieldElement Element>
constexpr auto AllElements() {
  constexpr int kFieldSize = utils::BinPow(Element::FieldBase(),
                                           Element::FieldPower());
  auto elements_range = Element::AllFieldElements();
  std::array<Element, kFieldSize> result{};
  std::copy(elements_range.begin(), elements_range.end(), result.begin());
  return result;
}

template <concepts::GaloisFieldElement Element, typename RandomGen>
Element GenElement(RandomGen& gen) {
  // TODO: use std discrete distribution
  // assume field to be relatively small
  constexpr static auto elements = AllElements<Element>();
  size_t index = gen() %  elements.size();
  return elements[index];
}

template <concepts::Polynom Poly, typename RandomGen>
Poly GenPoly(RandomGen& gen, size_t size) {
  using Element = typename Poly::Element;

  do {
    std::vector<Element> elements(size);
    for (auto& element : elements) {
      element = GenElement<Element>(gen);
    }
    Poly result(std::move(elements));
    if (!result.IsZero()) {
      return result;
    }
  } while (true);
}

}  // namespace factoriation