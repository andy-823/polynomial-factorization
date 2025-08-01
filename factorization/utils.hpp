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

#pragma  once

namespace factorization::utils {

template <class Iterator>
class IteratorRange {
 public:
 inline IteratorRange(Iterator begin, Iterator end) : begin_(begin), end_(end) {
  }

  inline Iterator begin() const {  // NOLINT
    return begin_;
  }
  inline Iterator end() const {  // NOLINT
    return end_;
  }

 private:
  Iterator begin_;
  Iterator end_;
};

template <typename Value, typename Power = Value>
constexpr inline Value BinPow(Value base, Power power) {
  Value result{1};
  while (power > 0) {
    if (power % 2 != 0) {
      result = result * base;
    }
    base = base * base;
    power /= 2;
  }
  return result;
}

} // namespace factorization::utils
