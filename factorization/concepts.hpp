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
concept GaloisField = requires (Field field, Field::Value value) {
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

  { field.FieldValueFromConstant(value) } -> std::same_as<typename Field::Value>;

  // need to iterate over all field elements
  { field.FirstFieldValue() } -> std::same_as<typename Field::Value>;
  { field.NextFieldValue(value) } -> std::same_as<typename Field::Value>;
  { field.LastFieldValue() } -> std::same_as<typename Field::Value>;
};

template <typename Element>
concept GaloisFieldElement = requires (Element element, Element::Value value) {
  std::copy_constructible<Element>;
  { element = element } -> std::same_as<Element&>;

  { Element::One() } -> std::same_as<Element>;
  { Element::Zero() } -> std::same_as<Element>;
  { Element::AsPolyConstant(value) } -> std::same_as<Element>;
  // get raw value of galois field element
  { element.Get() } -> std::same_as<typename Element::Value>;

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
concept Polynom = requires(Poly poly, Poly::Element value) {
  // This requirement is needed since sometimes we need
  // to perform arithmetic outside of polynom class
  GaloisFieldElement<typename Poly::Element>;

  // TODO: make construction from range
  std::constructible_from<Poly, std::vector<typename Poly::Element>>;
  std::copy_constructible<Poly>;
  { poly = poly } -> std::same_as<Poly&>;

  // need to put it into std::set or simply compare
  { poly == poly } -> std::same_as<bool>;
  { poly < poly } -> std::same_as<bool>;

  { poly += poly } -> std::same_as<Poly&>;
  { poly -= poly } -> std::same_as<Poly&>;
  { poly *= poly } -> std::same_as<Poly&>;
  { poly /= poly } -> std::same_as<Poly&>;
  { poly %= poly } -> std::same_as<Poly&>;

  { -poly } -> std::same_as<Poly>;

  { poly + poly } -> std::same_as<Poly>;
  { poly - poly } -> std::same_as<Poly>;
  { poly * poly } -> std::same_as<Poly>;
  { poly / poly } -> std::same_as<Poly>;
  { poly % poly } -> std::same_as<Poly>;

  // this specialization is needed since we probably don't want
  // to construct Polynom from single element just for this
  { poly += value } -> std::same_as<Poly&>;
  { poly -= value } -> std::same_as<Poly&>;
  { poly /= value } -> std::same_as<Poly&>;
  { poly *= value } -> std::same_as<Poly&>;

  { poly + value } -> std::same_as<Poly>;
  { value + poly } -> std::same_as<Poly>;
  { poly - value } -> std::same_as<Poly>;
  { value - poly } -> std::same_as<Poly>;
  { poly * value } -> std::same_as<Poly>;
  { value * poly } -> std::same_as<Poly>;
  { poly / value } -> std::same_as<Poly>;

  // This method has to follow this invariant
  //   a[0] + a[1] x + a[2] x^2 + ... + a[n] x^n
  // From lower power to higher
  { poly.GetElements() } -> std::same_as<std::vector<typename Poly::Element>>;
  // return polynom power + 1 if polynom is nonzero
  // otherwise return zero
  { poly.Size() } -> std::integral;
  { poly.MakeMonic() } -> std::same_as<void>;
  { poly.Derivative() } -> std::same_as<Poly>;

  // Need to check if polynomial is equal to one or zero without any construction
  { poly.IsOne() } -> std::same_as<bool>;
  { poly.IsZero() } -> std::same_as<bool>;
};

}  // namespace factorization::concepts