// MIT License
//
// Copyright (c) 2026 Andrei Ishutin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stack>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>
#include <factorization/polynomial/common.hpp>
#include <factorization/utils.hpp>

namespace factorization::ddf {

template <concepts::Polynom Polynom>
struct DistinctDegreeFactor {
  Polynom factor;
  int degree;
};

namespace naive {

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  using Element = typename Poly::Element;
  constexpr static auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

  std::vector<DistinctDegreeFactor<Poly>> result;

  Poly h = x;
  size_t degree = 1;
  auto mod = poly.BuildModulus(2 * poly.Size());

  while (2 * degree <= poly.Size() - 1) {
    h = polynomial::BinPowMod(std::move(h), kFieldSize, mod);
    Poly factor = poly.Gcd(h.Sub(x));
    if (!factor.IsOne()) {
      result.emplace_back(factor, static_cast<int>(degree));
      poly = std::move(poly).Div(factor).MakeMonic();
      // factorization is complete
      if (poly.IsOne()) {
        break;
      }
      mod = poly.BuildModulus(2 * poly.Size());
      h = std::move(h).Rem(mod);
    }
    ++degree;
  }
  if (!poly.IsOne()) {
    const int factor_degree = static_cast<int>(poly.Size() - 1);
    result.emplace_back(std::move(poly), factor_degree);
  }
  return result;
}

}  // namespace naive

enum StepsMode { kExactNtl, kSmallField };

template <concepts::Polynom Poly>
int Degree(const Poly& value) {
  return value.IsZero() ? 0 : static_cast<int>(value.Size()) - 1;
}

namespace ntl_like {

template <concepts::Polynom Poly>
int Degree(const Poly& value) {
  return value.IsZero() ? 0 : static_cast<int>(value.Size()) - 1;
}

template <concepts::Polynom Poly, StepsMode kMode = kSmallField>
class DistinctDegreeFactorizer {
  using Element = typename Poly::Element;
  using Modulus = typename Poly::Modulus;
  static constexpr auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

  // Each factor represents interval ((j - 1)l, jl].
  class GcdBuffer {
    constexpr static int kBufferSize = 4;

   public:
    explicit GcdBuffer(int interval_size)
        : interval_size_(interval_size) {
    }

    bool Empty() const {
      return size_ == 0;
    }

    bool Full() const {
      return size_ == kBufferSize;
    }

    void Add(int j, Poly factor) {
      assert(!Full());
      buf_[size_] = {j, std::move(factor)};
      ++size_;
    }

    bool Proceed(Poly& poly, std::vector<Poly>& F, const Modulus& mod) {
      if (Empty()) {
        return false;
      }
      Poly product = buf_[0].second;
      for (int i = 1; i < size_; ++i) {
        product = std::move(product).Mul(buf_[i].second).Rem(mod);
      }
      product = std::move(product).Gcd(poly);
      if (product.IsOne()) {
        size_ = 0;
        return false;
      }

      poly = std::move(poly).Div(product).MakeMonic();
      bool changed = false;
      int i = 0;
      int j = buf_[0].first;
      int min_degree = (j - 1) * interval_size_ + 1;
      while (i < size_ - 1 && 2 * min_degree <= Degree(product)) {
        Poly factor = buf_[i].second.Gcd(product);
        if (!factor.IsOne()) {
          F[j] = std::move(factor); 
          product = std::move(product).Div(F[j]).MakeMonic();
          changed = true;
        }
        ++i;
        ++j;
        min_degree += interval_size_;
      }

      if (!product.IsOne()) {
        if (i != size_ - 1) {
          j = (Degree(product) + interval_size_ - 1) / interval_size_;
        }
        F[j] = std::move(product);
        changed = true;
      }

      size_ = 0;
      return changed;
    }

   private:
    int size_ = 0;
    int interval_size_;
    std::array<std::pair<int, Poly>, kBufferSize> buf_;
  };

 public:
  explicit DistinctDegreeFactorizer(Poly poly)
      : poly_(std::move(poly)) {
    n = Degree(this->poly_);
    if constexpr (kMode == kExactNtl) {
      l = std::floor(std::sqrt(n / 2.0));
    } else {
      auto q = kFieldSize;
      l = std::floor(std::pow(n, 0.75L) / std::sqrt(3.0L * std::log2(q)));
    }
    if (l == 0) {
      l = 1;
    }
    m = (n + 2 * l - 1) / (2 * l);
  }

  // One-shot
  std::vector<DistinctDegreeFactor<Poly>> Run() {
    return RunWithObserver([](const char*, bool) {});
  }

  template <typename Observer>
  std::vector<DistinctDegreeFactor<Poly>> RunWithObserver(Observer&& observer) {
    observer("BuildModulus", true);
    auto mod = poly_.BuildModulus(2 * poly_.Size());
    observer("BuildModulus", false);
    observer("GenerateBabySteps", true);
    GenerateBabySteps(mod);
    observer("GenerateBabySteps", false);
    observer("GenerateGiantSteps", true);
    GenerateGiantSteps(mod);
    observer("GenerateGiantSteps", false);
    observer("GiantRefine", true);
    GiantRefine();
    observer("GiantRefine", false);
    observer("BabyRefine", true);
    BabyRefine();
    observer("BabyRefine", false);
    return std::move(result_);
  }

 private:
  void GenerateBabySteps(const Modulus& mod) {
    const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

    h.resize(l + 1);
    h[0] = x;
    h[1] = polynomial::BinPowMod(x, kFieldSize, mod);

    if constexpr (kMode == kExactNtl) {
      int t = std::floor(std::sqrt(n));
      if (t == 0) {
        t = 1;
      }
      auto matrix = polynomial::BuildCompModMatrix(h[1], t, mod);
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::CompMod(h[i - 1], matrix, mod);
      }
    } else {
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::BinPowMod(h[i - 1], kFieldSize, mod);
      }
    }
  }

  void GenerateGiantSteps(const Modulus& mod) {
    // H[0] is never used
    H.resize(m + 1);
    H[1] = h[l];

    int t = std::floor(std::sqrt(n));
    if (t == 0) {
      t = 1;
    }
    auto matrix = polynomial::BuildCompModMatrix(H[1], t, mod);
    for (int j = 2; j <= m; ++j) {
      // H_j = H_{j-1} (H_1) (mod f)
      H[j] = polynomial::CompMod(H[j - 1], matrix, mod);
    }
  }

  void GiantRefine() {
    // F[0] is never used btw
    F.assign(m + 1, Poly(Element::One()));
    std::vector<Poly> hh = h;

    Modulus mod = poly_.BuildModulus(2 * poly_.Size());
    GcdBuffer buffer(l);
    for (int j = 1; j <= m; ++j) {
      Poly HH = H[j].Rem(mod);  // NOLINT(readability-identifier-naming)
      // I = (H_j - h_0) * ... * (H_j - h_{l-1})
      Poly I = HH.Sub(hh[0]);  // NOLINT(readability-identifier-naming)
      for (int i = 1; i < l; ++i) {
        I = std::move(I).Mul(HH.Sub(hh[i])).Rem(mod);
      }
      buffer.Add(j, std::move(I));
      if (buffer.Full() && buffer.Proceed(poly_, F, mod)) {
        if (!poly_.IsOne() && Degree(poly_) >= 2 * j * l) {
          mod = poly_.BuildModulus(2 * poly_.Size());
          for (int i = 1; i <= l; ++i) {
            hh[i] = std::move(hh[i]).Rem(mod);
          }
        }
      }
      if (Degree(poly_) < 2 * j * l) {
        break;
      }
    }
    buffer.Proceed(poly_, F, mod);
    if (!poly_.IsOne()) {
      const int d = Degree(poly_);
      result_.emplace_back(std::move(poly_), d);
    }
  }

  void BabyRefine() {
    for (int j = 1; j < F.size(); ++j) {
      IntervalRefine(j, std::move(F[j]));
    }
  }

  void IntervalRefine(int j,
                      Poly F_j) {  // NOLINT(readability-identifier-naming)
    Modulus mod = F_j.BuildModulus(n);
    for (int i = l; i-- > 0;) {
      const int factor_degree = j * l - i;
      if (Degree(F_j) < 2 * factor_degree) {
        if (!F_j.IsOne()) {
          const int d = Degree(F_j);
          result_.emplace_back(std::move(F_j), d);
        }
        break;
      }
      // factor = Gcd(H_j - h_i, F_j)
      // bu since deg F_j is low we use rem at first
      Poly factor = H[j].Sub(h[i]).Rem(mod).Gcd(F_j);
      if (!factor.IsOne()) {
        F_j = std::move(F_j).Div(factor).MakeMonic();
        mod = F_j.BuildModulus(n);
        result_.emplace_back(std::move(factor), factor_degree);
      }
    }
  }

  Poly poly_;
  int n = 0;            // NOLINT(readability-identifier-naming)
  int l = 1;            // NOLINT(readability-identifier-naming)
  int m = 1;            // NOLINT(readability-identifier-naming)
  std::vector<Poly> h;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> H;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> F;  // NOLINT(readability-identifier-naming)
  std::vector<DistinctDegreeFactor<Poly>> result_;
};

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const size_t n = Degree(poly);
  if (n == 1) {
    return {{std::move(poly), 1}};
  }
  return DistinctDegreeFactorizer<Poly>(std::move(poly)).Run();
}

}  // namespace ntl_like

namespace own_lazy {

template <concepts::Polynom Poly, StepsMode kMode = kSmallField>
class DistinctDegreeFactorizer {
  using Element = typename Poly::Element;
  using Modulus = typename Poly::Modulus;
  static constexpr auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

 public:
  explicit DistinctDegreeFactorizer(Poly poly)
      : poly_(std::move(poly)) {
    n = Degree(this->poly_);
    if constexpr (kMode == kExactNtl) {
      l = std::floor(std::sqrt(n / 2.0));
    } else {
      auto q = kFieldSize;
      l = std::floor(std::pow(n, 0.75L) / std::sqrt(3.0L * std::log2(q)));
    }
    if (l == 0) {
      l = 1;
    }
    m = (n + 2 * l - 1) / (2 * l);
  }

  // One-shot
  std::vector<DistinctDegreeFactor<Poly>> Run() {
    Modulus mod = poly_.BuildModulus(2 * poly_.Size());

    GenerateBabySteps(mod);
    Poly H = h[l];  // NOLINT(readability-identifier-naming)

    int t = std::floor(std::sqrt(n));
    t = std::max(t, 1);
    auto matrix = polynomial::BuildCompModMatrix(H, t, mod);

    for (int j = 1; j <= m; ++j) {
      Poly I = H.Sub(h[0]);  // NOLINT(readability-identifier-naming)
      for (int i = 1; i < l; ++i) {
        I = std::move(I).Mul(H.Sub(h[i])).Rem(mod);
      }
      // NOLINTNEXTLINE(readability-identifier-naming)
      Poly F = std::move(I).Gcd(poly_);
      bool is_one = F.IsOne();
      if (!is_one) {
        poly_ = std::move(poly_).Div(F).MakeMonic();
        IntervalRefine(j, std::move(F), H);
      }
      if (Degree(poly_) < 2 * j * l) {
        break;
      }
      if (!is_one) {
        mod = poly_.BuildModulus(2 * poly_.Size());
        t = std::floor(std::sqrt(Degree(poly_)));
        t = std::max(t, 1);
        polynomial::UpdateCompModMatrix<Poly>(matrix, t, mod);
        H = std::move(H).Rem(mod);
        for (int i = 1; i <= l; ++i) {
          h[i] = std::move(h[i]).Rem(mod);
        }
      }
      H = polynomial::CompMod(H, matrix, mod);
    }
    if (!poly_.IsOne()) {
      const int d = Degree(poly_);
      result_.emplace_back(std::move(poly_), d);
    }
    return std::move(result_);
  }

 private:
  void GenerateBabySteps(const Modulus& mod) {
    const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

    h.resize(l + 1);
    h[0] = x;
    h[1] = polynomial::BinPowMod(x, kFieldSize, mod);

    if constexpr (kMode == kExactNtl) {
      int t = std::floor(std::sqrt(n));
      if (t == 0) {
        t = 1;
      }
      auto matrix = polynomial::BuildCompModMatrix(h[1], t, mod);
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::CompMod(h[i - 1], matrix, mod);
      }
    } else {
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::BinPowMod(h[i - 1], kFieldSize, mod);
      }
    }
  }

  void IntervalRefine(int j, Poly F,    // NOLINT(readability-identifier-naming)
                      const Poly& H) {  // NOLINT(readability-identifier-naming)
    Modulus mod = F.BuildModulus(n);
    for (int i = l; i-- > 0;) {
      const int factor_degree = j * l - i;
      if (Degree(F) < 2 * factor_degree) {
        if (!F.IsOne()) {
          const int d = Degree(F);
          result_.emplace_back(std::move(F), d);
        }
        break;
      }
      // factor = Gcd(H_j - h_i, F_j)
      // bu since deg F_j is low we use rem at first
      Poly factor = H.Sub(h[i]).Rem(mod).Gcd(F);
      if (!factor.IsOne()) {
        F = std::move(F).Div(factor).MakeMonic();
        mod = F.BuildModulus(n);
        result_.emplace_back(std::move(factor), factor_degree);
      }
    }
  }

  Poly poly_;
  int n = 0;            // NOLINT(readability-identifier-naming)
  int l = 1;            // NOLINT(readability-identifier-naming)
  int m = 1;            // NOLINT(readability-identifier-naming)
  std::vector<Poly> h;  // NOLINT(readability-identifier-naming)
  std::vector<DistinctDegreeFactor<Poly>> result_;
};

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const size_t n = Degree(poly);
  if (n == 1) {
    return {{std::move(poly), 1}};
  }
  return DistinctDegreeFactorizer<Poly>(std::move(poly)).Run();
}

}  // namespace own_lazy

namespace own_tree {

template <concepts::Polynom Poly, StepsMode kMode = kSmallField>
class DistinctDegreeFactorizer {
  using Element = typename Poly::Element;
  using Modulus = typename Poly::Modulus;
  static constexpr auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

  // one-shot tree
  class ComputationTree {
    struct Node {
      Poly forward;
      Poly backward;
      // (from, to]
      int from;
      int to;

      Node* left = nullptr;
      Node* right = nullptr;
    };

   public:
    ComputationTree(int interval_size, int interval_count, Modulus mod)
        : interval_size_(interval_size),
          nodes_(2 * interval_count),
          mod_(std::move(mod)) {
      for (auto& node : nodes_) {
        pool_.push(&node);
      }
    }

    // add leave
    void Add(int interval, Poly poly) {
      Node* node = GetNode();
      leaves_.push_back(node);

      node->forward = std::move(poly);
      node->from = interval - 1;
      node->to = interval;
      node->left = nullptr;
      node->right = nullptr;
    }

    Poly Build() {
      // nodes don't intersect
      auto compare_nodes = [](Node* lhs, Node* rhs) {
        return lhs->from < rhs->from;
      };
      std::sort(leaves_.begin(), leaves_.end(), compare_nodes);
      std::vector<Node*> current_level = leaves_;
      while (current_level.size() > 1) {
        std::vector<Node*> next_level;
        int k = current_level.size() / 2;
        for (int i = 0; i < k; ++i) {
          Node* next = GetNode();
          Merge(current_level[2 * i], current_level[2 * i + 1], next);
          next_level.push_back(next);
        }
        if (current_level.size() % 2 != 0) {
          next_level.push_back(current_level.back());
        }
        current_level = std::move(next_level);
      }
      root_ = current_level[0];
      return std::move(root_->forward);
    }

    void Run(Poly root_val) {
      root_->backward = std::move(root_val);
      BackwardVisitNode(root_);
    }

    std::vector<std::pair<int, Poly>> Extract() {
      std::vector<std::pair<int, Poly>> result;
      result.reserve(leaves_.size());
      for (Node* leave : leaves_) {
        result.emplace_back(leave->to, std::move(leave->backward));
      }
      return result;
    }

   private:
    void BackwardVisitNode(Node* node) {
      if (node->left == nullptr || node->right == nullptr) {
        return;
      }
      Split(node->left, node->right, node);
      BackwardVisitNode(node->left);
      BackwardVisitNode(node->right);
    }

    Node* GetNode() {
      Node* result = pool_.top();
      pool_.pop();
      return result;
    }

    void RetireNode(Node* node) {
      pool_.push(node);
    }

    void Merge(Node* left, Node* right, Node* root) {
      root->left = left;
      root->right = right;

      root->from = left->from;
      root->to = right->to;
      root->forward = std::move(right->forward).Mul(left->forward).Rem(mod_);
    }

    void Split(Node* left, Node* right, Node* root) {
      if (root->backward.IsOne()) {
        left->backward = root->backward;
        right->backward = root->backward;
        return;
      }
      int d = interval_size_ * root->from + 1;
      // that means that backward is single divisor
      if (Degree(root->backward) < 2 * d) {
        int border = interval_size_ * right->from + 1;
        if (Degree(root->backward) < border) {
          left->backward = std::move(root->backward);
          right->backward = Poly(Element::One());
        } else {
          right->backward = std::move(root->backward);
          left->backward = Poly(Element::One());
        }
        return;
      }
      left->backward = root->backward.Gcd(std::move(left->forward));
      right->backward = root->backward.Div(left->backward);
    }

    const int interval_size_;
    std::vector<Node> nodes_;
    std::stack<Node*> pool_;
    std::vector<Node*> leaves_;

    const Modulus mod_;
    Node* root_;
  };

 public:
  explicit DistinctDegreeFactorizer(Poly poly)
      : poly_(std::move(poly)) {
    n = Degree(this->poly_);
    if constexpr (kMode == kExactNtl) {
      l = std::floor(std::sqrt(n / 2.0));
    } else {
      auto q = kFieldSize;
      l = std::floor(std::pow(n, 0.75L) / std::sqrt(3.0L * std::log2(q)));
    }
    if (l == 0) {
      l = 1;
    }
    m = (n + 2 * l - 1) / (2 * l);
  }

  // One-shot
  std::vector<DistinctDegreeFactor<Poly>> Run() {
    return RunWithObserver([](const char*, bool) {});
  }

  template <typename Observer>
  std::vector<DistinctDegreeFactor<Poly>> RunWithObserver(Observer&& observer) {
    observer("BuildModulus", true);
    Modulus mod = poly_.BuildModulus(2 * poly_.Size());
    observer("BuildModulus", false);
    observer("GenerateBabySteps", true);
    GenerateBabySteps(mod);
    observer("GenerateBabySteps", false);
    observer("GenerateGiantSteps", true);
    GenerateGiantSteps(mod);
    observer("GenerateGiantSteps", false);
    observer("GiantRefine", true);
    GiantRefine();
    observer("GiantRefine", false);
    observer("BabyRefine", true);
    BabyRefine();
    observer("BabyRefine", false);
    return std::move(result_);
  }

 private:
  void GenerateBabySteps(const Modulus& mod) {
    const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

    h.resize(l + 1);
    h[0] = x;
    h[1] = polynomial::BinPowMod(x, kFieldSize, mod);

    if constexpr (kMode == kExactNtl) {
      int t = std::floor(std::sqrt(n));
      if (t == 0) {
        t = 1;
      }
      auto matrix = polynomial::BuildCompModMatrix(h[1], t, mod);
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::CompMod(h[i - 1], matrix, mod);
      }
    } else {
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::BinPowMod(h[i - 1], kFieldSize, mod);
      }
    }
  }

  void GenerateGiantSteps(const Modulus& mod) {
    // H[0] is never used
    H.resize(m + 1);
    H[1] = h[l];

    int t = std::floor(std::sqrt(n));
    if (t == 0) {
      t = 1;
    }
    auto matrix = polynomial::BuildCompModMatrix(H[1], t, mod);
    for (int j = 2; j <= m; ++j) {
      // H_j = H_{j-1} (H_1) (mod f)
      H[j] = polynomial::CompMod(H[j - 1], matrix, mod);
    }
  }

  void GiantRefine() {
    // F[0] is never used btw
    F.assign(m + 1, Poly(Element::One()));
    Modulus mod = poly_.BuildModulus(2 * poly_.Size());
    ComputationTree tree(l, m, mod);
    for (int j = 1; j <= m; ++j) {
      // I = (H_j - h_0) * ... * (H_j - h_{l-1})
      Poly I = H[j].Sub(h[0]);  // NOLINT(readability-identifier-naming)
      for (int i = 1; i < l; ++i) {
        I = std::move(I).Mul(H[j].Sub(h[i])).Rem(mod);
      }
      tree.Add(j, std::move(I));
    }
    Poly product = tree.Build();
    product = std::move(product).Gcd(poly_);
    poly_ = std::move(poly_).Div(product);
    tree.Run(std::move(product));

    for (auto [j, factor] : tree.Extract()) {
      F[j] = std::move(factor);
    }
    if (!poly_.IsOne()) {
      const int d = Degree(poly_);
      result_.emplace_back(std::move(poly_), d);
    }
  }

  void BabyRefine() {
    for (int j = 1; j < F.size(); ++j) {
      IntervalRefine(j, std::move(F[j]));
    }
  }

  void IntervalRefine(int j,
                      Poly F_j) {  // NOLINT(readability-identifier-naming)
    Modulus mod = F_j.BuildModulus(n);
    for (int i = l; i-- > 0;) {
      const int factor_degree = j * l - i;
      if (Degree(F_j) < 2 * factor_degree) {
        if (!F_j.IsOne()) {
          const int d = Degree(F_j);
          result_.emplace_back(std::move(F_j), d);
        }
        break;
      }
      // factor = Gcd(H_j - h_i, F_j)
      // bu since deg F_j is low we use rem at first
      Poly factor = H[j].Sub(h[i]).Rem(mod).Gcd(F_j);
      if (!factor.IsOne()) {
        F_j = std::move(F_j).Div(factor).MakeMonic();
        mod = F_j.BuildModulus(n);
        result_.emplace_back(std::move(factor), factor_degree);
      }
    }
  }

  Poly poly_;
  int n = 0;            // NOLINT(readability-identifier-naming)
  int l = 1;            // NOLINT(readability-identifier-naming)
  int m = 1;            // NOLINT(readability-identifier-naming)
  std::vector<Poly> h;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> H;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> F;  // NOLINT(readability-identifier-naming)
  std::vector<DistinctDegreeFactor<Poly>> result_;
};

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const size_t n = Degree(poly);
  if (n == 1) {
    return {{std::move(poly), 1}};
  }
  return DistinctDegreeFactorizer<Poly>(std::move(poly)).Run();
}

}  // namespace own_tree

}  // namespace factorization::ddf
