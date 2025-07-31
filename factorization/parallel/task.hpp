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

#include <memory>
#include <utility>

namespace factorization::parallel {

namespace internal {

class UniqueFunction {
 private:
  // For type erasure
  struct IRunnable {
    virtual ~IRunnable() = default;

    virtual void Run() = 0;
  };

  template <typename F>
  struct Runnable : public IRunnable {
    explicit Runnable(F body)
        : body_(std::move(body)) {
    }
    virtual void Run() {
      body_();
    }

    F body_;
  };

 public:
  template <typename F>
  explicit UniqueFunction(F&& body) {
    body_ = std::unique_ptr<IRunnable>(
        new Runnable<std::remove_cvref_t<F>>(std::forward<F>(body)));
  }

  // Movable
  UniqueFunction(UniqueFunction&&) = default;
  UniqueFunction& operator=(UniqueFunction&&) = default;
  // Non-copyable
  UniqueFunction(const UniqueFunction&) = delete;
  UniqueFunction& operator=(const UniqueFunction&) = delete;

  void operator()() {
    body_->Run();
  }

 private:
  std::unique_ptr<IRunnable> body_;
};

}  // namepsace internal

struct ITask {
  virtual void Run() noexcept = 0;

 protected:
  virtual ~ITask() = default;
};

class Task : public ITask {
  using Body = internal::UniqueFunction;

 public:
  template <typename F>
  explicit Task(F body)
      : body_(std::move(body)) {
  }

  virtual void Run() noexcept {
    try {
      body_();
    } catch (...) {
    }
    delete this;
  }

 private:
  Body body_;
};


}  // factorization::parallel