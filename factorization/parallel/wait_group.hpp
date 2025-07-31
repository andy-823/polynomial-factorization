
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