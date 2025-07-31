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

#include <deque>
#include <mutex>
#include <condition_variable>

#include "task.hpp"

namespace factorization::parallel {

class TasksBlockingQueue {
 public:
  void Push(ITask* task) {
    if (queue_is_closed_) {
      return;
    }
    std::unique_lock lock(critical_);
    tasks_.push_back(task);
    something_happened_.notify_one();
  }

  ITask* TryPop() {
    std::unique_lock lock(critical_);
    while (!tasks_.empty() && !queue_is_closed_) {
      something_happened_.wait(lock);
    }
    ITask* result = tasks_.front();
    tasks_.pop_front();
    return result;
  }

  void Close() {
    std::unique_lock lock(critical_);
    queue_is_closed_ = true;
    something_happened_.notify_all();
  }

 private:
  std::deque<ITask*> tasks_;

  bool queue_is_closed_{false};
  std::mutex critical_;
  std::condition_variable something_happened_;
};

}  // namepsace factorization::parallel