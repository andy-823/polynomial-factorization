// MIT License
//
// Copyright (c) 2025 Andrei Ishutin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <concepts>
#include <ranges>
#include <vector>

namespace factorization::concepts {

template <typename Field>
concept GaloisField = requires (const Field& field,
                                typename Field::Value value,
                                typename Field::Coefficient coef) {
  typename Field::Coefficient;
  requires std::integral<typename Field::Coefficient>;
  // return addition and multiplication neutral elements respectively
  { field.Zero() } -> std::same_as<typename Field::Value>;
  { field.One() } -> std::same_as<typename Field::Value>;
  // perform arithmetic operations over galois field
  { field.Add(value, value) } -> std::same_as<typename Field::Value>;
  { field.Sub(value, value) } -> std::same_as<typename Field::Value>;
  { field.Negative(value) } -> std::same_as<typename Field::Value>;
  { field.Multiply(value, value) } -> std::same_as<typename Field::Value>;
  { field.Divide(value, value) } -> std::same_as<typename Field::Value>;
  { field.Inverse(value) } -> std::same_as<typename Field::Value>;
  { field.Pow(value, 0) } -> std::same_as<typename Field::Value>;
  // getters for field base and dimension
  { Field::FieldBase() } -> std::integral;
  { Field::FieldPower() } -> std::integral;

  { field.Encode(std::array<typename Field::Coefficient, Field::FieldPower()>())}
    -> std::same_as<typename Field::Value>;
  { field.Encode(coef)} -> std::same_as<typename Field::Value>;
  { field.Decode(value) }
    -> std::same_as<std::array<typename Field::Coefficient, Field::FieldPower()>>;
};

template <typename Element>
concept GaloisFieldElement = requires (Element element,
                                       typename Element::Coefficient coeff) {
  requires std::semiregular<Element>;
  requires std::totally_ordered<Element>;

  requires std::constructible_from<std::array<typename Element::Coefficient, Element::FieldPower()>>;
  requires std::constructible_from<typename Element::Coefficient>;

  { Element::One() } -> std::same_as<Element>;
  { Element::Zero() } -> std::same_as<Element>;
  // get raw value of galois field element
  { element.Get() }
    -> std::same_as<std::array<typename Element::Coefficient, Element::FieldPower()>>;

  // arithmetic operations that can be performed
  { element += element } -> std::same_as<Element&>;
  { element -= element } -> std::same_as<Element&>;
  { element *= element } -> std::same_as<Element&>;
  { element /= element } -> std::same_as<Element&>;

  { element + element } -> std::same_as<Element>;
  { element - element } -> std::same_as<Element>;
  { element * element } -> std::same_as<Element>;
  { element / element } -> std::same_as<Element>;

  { -element } -> std::same_as<Element>;
  // C++ doesn't have operator for getting inverse value or power
  { element.Inverse() } -> std::same_as<Element>;
  { element.Pow(0) } -> std::same_as<Element>;

  { Element::FieldBase() } -> std::integral;
  { Element::FieldPower() } -> std::integral;

  { Element::AllFieldElements() } -> std::ranges::range;
};

template <typename Poly>
concept Polynom = requires(const Poly& poly, typename Poly::Element value) {
  // This requirement is needed since sometimes we need
  // to perform arithmetic outside of polynom class
  requires GaloisFieldElement<typename Poly::Element>;

  requires std::constructible_from<Poly, std::vector<typename Poly::Element>>;
  requires std::semiregular<Poly>;
  requires std::totally_ordered<Poly>;

  // Here methods are used instead of operators
  // It was made for better flexibility
  { poly.Add(poly) } -> std::same_as<Poly>;
  { poly.Sub(poly) } -> std::same_as<Poly>;
  { poly.Mul(poly) } -> std::same_as<Poly>;
  { poly.Div(poly) } -> std::same_as<Poly>;
  { poly.Rem(poly) } -> std::same_as<Poly>;
  // <quotient, remainder>
  { poly.DivRem(poly) } -> std::same_as<std::pair<Poly, Poly>>;
  { poly.Gcd(poly) } -> std::same_as<Poly>;
  { poly.MakeMonic() } -> std::same_as<Poly>;
  { poly.Derivative() } -> std::same_as<Poly>;

  { poly.Add(value) } -> std::same_as<Poly>;
  { poly.Sub(value) } -> std::same_as<Poly>;
  { poly.Mul(value) } -> std::same_as<Poly>;
  { poly.Div(value) } -> std::same_as<Poly>;

  // This method has to follow this invariant
  //   a[0] + a[1] x + a[2] x^2 + ... + a[n] x^n
  // From lower power to higher
  { poly.Get() } -> std::same_as<std::vector<typename Poly::Element>>;
  // return polynom power + 1 if polynom is nonzero
  // otherwise return zero
  { poly.Size() } -> std::integral;
  { poly.IsOne() } -> std::same_as<bool>;
  { poly.IsZero() } -> std::same_as<bool>;
};

}  // namespace factorization::concepts