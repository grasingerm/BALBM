#ifndef LATTICE_HH
#define LATTICE_HH

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

// TODO: consider using a different data structure for 2D and 3D arrays
// TODO: consider overloading [] for indexing???

#include "balbm_config.hh"
#include "helpers/mem_helpers.hh"
#include "node_desc.hh"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exception>
//#include <iosfwd>
#include <memory>
#include <sstream>
#include <vector>

namespace balbm {

namespace d2q9 {

// TODO: "More code (ALWAYS) runs slower" -- John Lakos --

//! \class Lattice
//!
//! \brief  lattice for the lattice Boltzmann method
class Lattice {
public:
  // constructors and assignment
  // TODO: make more constructors, initializers, and factories
  Lattice() : ni_(0), nj_(0), spf_(nullptr), spftemp_(nullptr) {}
  Lattice(const unsigned ni, const unsigned nj, const double rho = 1.0)
      : ni_(ni), nj_(nj), spf_(new double[ni * nj * num_k()]),
        spftemp_(new double[ni * nj * num_k()]), node_descs_(ni * nj),
        mem_pool_(max_node_desc_size() * ni * nj) {
    init_f_(rho);
  }
  Lattice(const Lattice &);
  Lattice &operator=(const Lattice &);
  Lattice(Lattice &&);
  Lattice &operator=(Lattice &&);
  ~Lattice() {
    try {
      for (auto &pnode_desc : node_descs_)
        pnode_desc->~AbstractNodeDesc();
    } catch (...) {
    }
  }

  // accessors
  // static and constant expression values
  static constexpr double dx() { return 1.0; }
  static constexpr double dt() { return 1.0; }
  static constexpr double c() { return dx() / dt(); }
  static constexpr double cs() { return c() / sqrt(3.0); }
  static constexpr double cssq() { return cs() * cs(); }
  inline unsigned num_i() const { return ni_; }
  inline unsigned num_j() const { return nj_; }
  static constexpr unsigned num_k() { return nk_; }
  inline const double *pf() const noexcept { return spf_.get(); }
  inline double f(unsigned i, unsigned j, unsigned k) const noexcept {
    assert(k < 9 && static_cast<int>(k) >= 0 &&
           "index `k` out of bounds in Lattice::pc");
    assert(in_bounds(i, j) && "out of bounds in Lattice::f");
    return spf_[ni_ * (i * nj_ + j) + k];
  }
  inline const double *pftemp() const noexcept { return spftemp_.get(); }
  inline double ftemp(unsigned i, unsigned j, unsigned k) const noexcept {
    assert(k < 9 && static_cast<int>(k) >= 0 &&
           "index `k` out of bounds in Lattice::ftemp");
    assert(in_bounds(i, j) && "out of bounds in Lattice::ftemp");
    return spftemp_[ni_ * (i * nj_ + j) + k];
  }
  inline double *pf(const unsigned i, const unsigned j) {
    assert(in_bounds(i, j) && "out of bounds in Lattice::pf");
    return &(spf_[(i * nj_ + j) * num_k()]);
  }
  inline double &f(const unsigned i, const unsigned j, const unsigned k) {
    assert(k < 9 && static_cast<int>(k) >= 0 &&
           "index `k` out of bounds in Lattice::f");
    assert(in_bounds(i, j) && "out of bounds in Lattice::f");
    return *(pf(i, j) + k);
  }
  inline double *pft(const unsigned i, const unsigned j) {
    assert(in_bounds(i, j) && "out of bounds in Lattice::pft");
    return &(spftemp_[(i * nj_ + j) * num_k()]);
  }
  inline double &ft(const unsigned i, const unsigned j, const unsigned k) {
    assert(k < 9 && static_cast<int>(k) >= 0 &&
           "index `k` out of bounds in Lattice::ft");
    assert(in_bounds(i, j) && "out of bounds in Lattice::ft");
    return *(pft(i, j) + k);
  }
  inline const std::vector<AbstractNodeDesc *> &node_descs() const noexcept {
    return node_descs_;
  }
  inline const AbstractNodeDesc &node_desc(const unsigned i,
                                           const unsigned j) const {
    assert(in_bounds(i, j) && "out of bounds in Lattice::node_desc");
    return *(node_descs_[nj_ * i + j]);
  }
  template <typename Node, typename... Args>
  inline void set_node_desc(const unsigned i, const unsigned j, Args... args) {
#ifndef NDEBUG
    assert(in_bounds(i, j) && "out of bounds in Lattice::set_node_desc");
    AbstractNodeDesc *pnd = mem_pool_.allocate<Node>(args...);
    assert(pnd != nullptr);
    node_descs_[nj_ * i + j] = pnd;
#else
    node_descs_[nj_ * i + j] = mem_pool_.allocate<Node>(args...);
#endif
  }
  inline const double *pc(const unsigned k) const noexcept {
    assert(k <= 9 && static_cast<int>(k) >= 0 &&
           "index `k` out of bounds in Lattice::pc");
    return lat_vecs_[k];
  }
  inline double c(const unsigned k, const unsigned c) const noexcept {
    assert(k <= 9 && static_cast<int>(k) >= 0 &&
           "index `k` out of bounds in Lattice::c");
    assert(c <= 2 && static_cast<int>(c) >= 0 &&
           "index `c` out of bounds in Lattice::c");
    return *(pc(k) + c);
  }
  inline double w(const unsigned k) const noexcept {
    assert(k <= 9 && static_cast<int>(k) >= 0 &&
           "index `k` out of bounds in Lattice::w");
    return w_[k];
  }

  // mutators
  // stream
  inline void stream(const unsigned i, const unsigned j) {
    assert(in_bounds(i, j) && "out of bounds in Lattice::stream");
    node_desc(i, j).stream(*this, i, j);
  }
  void stream(const unsigned bi, const unsigned ei, const unsigned bj,
              const unsigned ej) {
    for (unsigned i = bi; i <= ei; ++i)
      for (unsigned j = bj; j <= ej; ++j)
        stream(i, j);
  }
  inline void stream() { stream(0, ni_ - 1, 0, nj_ - 1); }
  void stream(const std::vector<std::array<unsigned, 4>> &);

  // collide
  inline void collide_and_bound(IncompFlowMultiscaleMap &mmap,
                                const IncompFlowCollisionManager &cman,
                                const unsigned i, const unsigned j) {
    assert(in_bounds(i, j) && "out of bounds in Lattice::collide_and_bound");
    node_desc(i, j).collide_and_bound(*this, mmap, cman, i, j);
  }
  void collide_and_bound(IncompFlowMultiscaleMap &mmap,
                         const IncompFlowCollisionManager &cman,
                         const unsigned bi, const unsigned ei,
                         const unsigned bj, const unsigned ej) {
    for (unsigned i = bi; i <= ei; ++i)
      for (unsigned j = bj; j <= ej; ++j)
        collide_and_bound(mmap, cman, i, j);
  }
  inline void collide_and_bound(IncompFlowMultiscaleMap &mmap,
                                const IncompFlowCollisionManager &cman) {
    collide_and_bound(mmap, cman, 0, ni_ - 1, 0, nj_ - 1);
  }
  void collide_and_bound(IncompFlowMultiscaleMap &,
                         const IncompFlowCollisionManager &,
                         const std::vector<std::array<unsigned, 4>> &);

  inline void swap_f_ptrs() { spf_.swap(spftemp_); }

  // bounds checking
  inline bool in_bounds(const int i, const int j) const noexcept {
    return (i < static_cast<int>(ni_) && i >= 0 && j < static_cast<int>(nj_) &&
            j >= 0);
  }
  bool check_bounds(const unsigned, const unsigned) const
      throw(std::out_of_range);

private:
  static constexpr unsigned nk_ = 9;
  static const double lat_vecs_[nk_][2];
  static const double w_[nk_];
  unsigned ni_;
  unsigned nj_;
  std::unique_ptr<double[]> spf_;
  std::unique_ptr<double[]> spftemp_;
  std::vector<AbstractNodeDesc *> node_descs_;
  SimpleMemPool mem_pool_;

  void init_f_(const double);
};

} // namespace d2q9

} // namespace balbm

#endif
