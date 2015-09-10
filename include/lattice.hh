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

#include "balbm_config.hh"
#include "helpers/mem_helpers.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace balbm
{

namespace d2q9
{

//TODO: "More code (ALWAYS) runs slower" -- John Lakos --

//! \class Lattice
//!
//! \brief  lattice for the lattice Boltzmann method
//!
//! 6     2     5
//!   \   |   /
//!    \  |  /
//!     \ | /
//! 3 --- 0 --- 1
//!     / | \
//!    /  |  \
//!   /   |   \
//! 7     4     8
//!
class Lattice
{
  friend class NodeDesc;
public:
  // constructors and assignment
  // TODO: make more constructors, initializers, and factories
  Lattice() : ni_(0), nj_(0), f_(nullptr), ftemp_(nullptr) {}
  Lattice(const unsigned ni, const unsigned nj, const double rho = 1.0) 
    : nx_(ni), ny_(nj),
      f_(std::make_unique(new double[ni * nj * num_k()])), 
      ftemp_(std::make_unique(new double[ni * nj * num_k()])),
      node_descs_(std::vector<AbstractNodeDesc*>(ni * nj)),
      mem_pool_(SimpleMemPool(max_node_desc_size() * ni * nj))
    { init_f_(rho); }
  Lattice(const Lattice&);
  Lattice& operator=(const Lattice&);
  Lattice(Lattice&&);
  Lattice& operator=(Lattice&&);
  ~Lattice() 
  { 
    try 
    { 
      for (auto& node_desc : node_descs_) node_desc->::~AbstractNodeDesc(); 
    }
    catch(...) {} 
  }

  // accessors
    // static and constant expression values
  static constexpr double dx() const            { return dx_; }
  static constexpr double dt() const            { return dt_; }
  static constexpr double c() const             { return dx_ / dt_; }
  static constexpr double cs() const            { return c() / sqrt(3.0); }
  static constexpr double cssq() const          { return cs() * cs(); }
  inline unsigned num_i() const                 { return ni_; }
  inline unsigned num_j() const                 { return nj_; }
  static constexpr unsigned num_k() const       { return 9; }
  inline const double* pf() const noexcept      { return spf_.get(); }
  inline double f(unsigned i, unsigned j, unsigned k) const noexcept     
                                                { return f_(i, j, k); }
  inline const double* pftemp() const noexcept  { return spftemp_.get(); }
  inline double ftemp(unsigned i, unsigned j, unsigned k) const noexcept 
                                                { return ftemp_(i, j, k); }
  inline const std::vector<std::unique_ptr<NodeDesc>>& node_descs ()
    const noexcept                              { return node_descs_; }
  inline const AbstractNodeDesc& node_desc(const unsigned i, const unsigned j) 
    const { return *(node_descs_[nj_ * i + j]); }
  inline AbstractNodeDesc& node_desc(const unsigned i, const unsigned j)
    { return *(node_descs_[nj_ * i + j]); }
  inline const double* pc(const unsigned k) const noexcept 
                                                { return (&lat_vecs_[2 * k]); }
  inline double c(const unsigned k, const unsigned c) const noexcept 
                                                { return *(pc(k) + c); }
  inline double w(const unsigned k) const noexcept
                                                { return w_[k]; }
    
  // mutators
    // stream
  inline void stream(const unsigned i, const unsigned j) 
    { node_desc(i, j).stream(*this, i, j); }
  inline void stream(const unsigned bi, const unsigned ei, const unsigned bj, 
                     const unsigned ej)
    { 
      for (unsigned i = bi; i <= ei; ++i)
        for (unsigned j = bj; j <= ej; ++j)
          stream(i, j);
    }
  inline void stream() { stream(0, ni_ - 1, 0, nj_ - 1); }
  void stream(const std::vector<std::array<unsigned, 4>>&);

    // collide
  inline void collide_and_bound
    (IncompFlowMultiscaleMap& mmap, const IncompFlowCollisionManager& cman,
     const unsigned i, const unsigned j) 
    { node_desc(i, j).collide_and_bound(*this, mmap, cman, i, j); }
  void collide_and_bound 
    (IncompFlowMultiscaleMap& mmap, const IncompFlowCollisionManager& cman,
     const unsigned bi, const unsigned ei, const unsigned bj, const unsigned ej)
    { 
      for (unsigned i = bi; i <= ei; ++i)
        for (unsigned j = bj; j <= ej; ++j)
          collide_and_bound(i, j);
    }
  inline void collide_and_bound
    (IncompFlowMultiscaleMap& mmap, const IncompFlowCollisionManager& cman)
    { collide_and_bound(mmap, cman, 0, ni_ - 1, 0, nj_ - 1); }
  void collide_and_bound
    (IncompFlowMultiscaleMap&, const IncompFlowCollisionManager&
     const std::vector<std::array<unsigned, 4>>&);

  inline void swap_f_ptrs() { spf_.swap(spftemp_); }

  // bounds checking
  inline bool in_bounds(const unsigned i, const unsigned j) const noexcept 
    { return (i < ni && i >= 0 && j < nj && j >= 0); }
  void check_bounds(const unsigned i, const unsigned j) 
    const throw(std::out_of_range);

private:
  static const double[num_k()][2] lat_vecs_;
  static const double[num_k()] w_;
  static const double dx_;
  static const double dt_;
  const unsigned ni_;
  const unsigned nj_;
  std::unique_ptr<double[]> spf_;
  std::unique_ptr<double[]> spftemp_;
  std::vector<AbstractNodeDesc*> node_descs_;
  SimpleMemPool mem_pool_;

  // lattice vectors and particle distribution functions
  inline double* pf_(const unsigned i, const unsigned j) 
    { return &(spf_[(i*nj_ + j)*num_k()]); }
  inline double& f_(const unsigned i, const unsigned j, const unsigned k)
    { return *(pf_(i, j) + k); }
  inline double* pft_(const unsigned i, const unsigned j) 
    { return &(spftemp_[(i*nj_ + j)*num_k()]); }
  inline double& ft_(const unsigned i, const unsigned j, const unsigned k)
    { return *(pft_(i, j) + k); }

  void init_f_(const double rho);
};

} // namespace d2q9

} // namespace balbm

#endif
