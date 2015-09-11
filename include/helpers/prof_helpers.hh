#ifndef PROF_HELPERS_HH
#define PROF_HELPERS_HH

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

#include <iostream>
#include <chrono>
#include <ctime>
#include <tuple>
#include <functional>

namespace baprof {

static std::chrono::time_point<std::chrono::high_resolution_clock> tic_start__ =
    std::chrono::time_point<std::chrono::high_resolution_clock>::min();

//! start stop watch
inline void tic() { tic_start__ = std::chrono::high_resolution_clock::now(); }

//! stop stop watch. tell time
void toc() {
  if (tic_start__ ==
      std::chrono::time_point<std::chrono::high_resolution_clock>::min())
    throw std::logic_error("toc() called before tic()");
  const std::chrono::duration<double> elapsed =
      std::chrono::high_resolution_clock::now() - tic_start__;
  std::cout << "Time elapsed: " << elapsed.count() << " seconds.\n";
  tic_start__ = std::chrono::high_resolution_clock::now();
}

//! Profile a simple function given a function pointer
//!
//! \param f Function pointer
//! \param args Function arguments
//! \return (function result, function duration)
template <typename Result, typename... Sig>
std::tuple<Result, std::chrono::duration<double>> profile(Result (*f)(Sig...),
                                                          Sig... args) {
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::high_resolution_clock::now();
  Result result = f(args...);
  end = std::chrono::high_resolution_clock::now();

  return std::tuple<Result, std::chrono::duration<double>>(result, end - start);
}

//! Profile a simple function given a functor
//!
//! \param f Functor object
//! \param args Function arguments
//! \return (functor result, function duration)
template <typename Result, typename... Sig>
std::tuple<Result, std::chrono::duration<double>>
profile(const std::function<Result(Sig...)> &f, Sig... args) {
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::high_resolution_clock::now();
  Result result = f(args...);
  end = std::chrono::high_resolution_clock::now();

  return std::tuple<Result, std::chrono::duration<double>>(result, end - start);
}

//! Profile a simple function given a function pointer
//!
//! \param f Function pointer
//! \param args Function arguments
//! \return function duration
template <typename... Sig>
std::chrono::duration<double> profile_void(void (*f)(Sig...), Sig... args) {
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::high_resolution_clock::now();
  f(args...);
  end = std::chrono::high_resolution_clock::now();

  return end - start;
}

//! Profile a simple function given a functor
//!
//! \param f Functor object
//! \param args Function arguments
//! \return function duration
template <typename... Sig>
std::chrono::duration<double> profile_void(const std::function<void(Sig...)> &f,
                                           Sig... args) {
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::high_resolution_clock::now();
  f(args...);
  end = std::chrono::high_resolution_clock::now();

  return end - start;
}

} // namespace baprof

#endif // PROF_HELPERS_HH
