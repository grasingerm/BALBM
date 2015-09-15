#ifndef MEM_HELPERS_HH
#define MEM_HELPERS_HH

// Complex flow simulator using lattice Boltzmann method
// Copyright (C) 2015 Matthew Grasinger
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// A copy of the GNU General Public License is at the root directory of
// this program.  If not, see <http://www.gnu.org/licenses/>

#include <vector>

class SimpleMemPool {
public:
  SimpleMemPool(const std::size_t bytes = 1) : allocated_(0) {
    mem_.reserve(bytes);
    internal_ptr_ = raw_ptr();
  }
  inline void *raw_ptr() { return static_cast<void *>(mem_.data()); }
  inline std::size_t capacity() const { return mem_.capacity(); }
  template <typename T, typename... Args> T *allocate(Args...);

private:
  std::size_t allocated_;
  std::vector<char> mem_;
  void *internal_ptr_;
};

//! Allocate object in the memory pool
//!
//! \param ... args Constructor arguments
//! \return Pointer to object constructed into memory pool
template <typename T, typename... Args>
T *SimpleMemPool::allocate(Args... args) {
  std::size_t to_alloc = sizeof(T);
  if (to_alloc + allocated_ > capacity())
    return nullptr;

  T *temp_ptr = static_cast<T *>(internal_ptr_);
  new (temp_ptr) T(args...);
  T *result = temp_ptr;
  internal_ptr_ = static_cast<void *>(++temp_ptr);
  allocated_ += to_alloc;

  return result;
}

#endif // MEM_HELPERS_HH
