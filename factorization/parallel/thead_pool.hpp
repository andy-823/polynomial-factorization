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

#include <cassert>
#include <vector>
#include <thread>

#include "task.hpp"
#include "queue.hpp"

namespace factorization::parallel {

struct IScheduler {
  virtual void Submit(ITask*) = 0;

 protected:
  ~IScheduler() = default;
};

template <typename F>
void SubmitTask(IScheduler* scheduler, F task) {
  scheduler->Submit(new Task(std::move(task)));
}

class ThreadPool final : public IScheduler {
 public:
  explicit ThreadPool(size_t threads) : workers_(threads) {
  }

  ~ThreadPool() {
    assert(workers_.empty());
  }

  // Non-copyable
  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  // Non-movable
  ThreadPool(ThreadPool&&) = delete;
  ThreadPool& operator=(ThreadPool&&) = delete;

  void Start() {
    for (auto& worker : workers_) {
      worker = std::thread([&] {
        while (auto task = tasks_.TryPop()) {
          task->Run();
        }
      });
    }
  }

  // task::IScheduler
  void Submit(ITask* task) override {
    tasks_.Push(task);
  }

  void Stop() {
    tasks_.Close();
    for (auto& worker : workers_) {
      worker.join();
    }
    workers_.resize(0);
  }

 private:
  std::vector<std::thread> workers_;
  TasksBlockingQueue tasks_;
};

}  // namespace factorization::parallel