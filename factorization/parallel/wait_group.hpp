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

#include <mutex>
#include <condition_variable>

namespace factorization::parallel {

class WaitGroup {
 public:
  void Add(size_t count) {
    std::unique_lock lock(jobs_done_mutex_);
    n_jobs_ += count;
  }

  void Done() {
    std::unique_lock lock(jobs_done_mutex_);
    if (--n_jobs_ == 0 && n_wait_ != 0) {
      jobs_done_.notify_all();
    }
  }

  void Wait() {
    std::unique_lock lock(jobs_done_mutex_);
    ++n_wait_;
    while (n_jobs_ != 0) {
      jobs_done_.wait(lock);
    }
    --n_wait_;
  }

 private:
  size_t n_wait_{0};
  size_t n_jobs_{0};

  std::mutex jobs_done_mutex_;
  std::condition_variable jobs_done_;
};
 
}  // factorization::parallel