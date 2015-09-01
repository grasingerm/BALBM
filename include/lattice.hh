#ifndef __LATTICE_HH__
#define __LATTICE_HH__

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

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include "balbm_config.hh"

namespace balbm
{

namespace d2q9
{

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
  Lattice() : nx_(0), ny_(0), f_(nullptr), ftemp_(nullptr) {}
  Lattice(unsigned nx, unsigned ny) : nx_(nx), ny_(ny),
    f_(std::make_unique(new double[nx * ny * num_k()])), 
    ftemp_(std::make_unique(new double[nx * ny * num_k()])),
    node_descs_(std::vector<std::unique_ptr<NodeDesc>>(nx * ny * num_k())) {}
  Lattice(const Lattice&);
  Lattice& operator=(const Lattice&);
  Lattice(Lattice&&);
  Lattice& operator=(Lattice&&);
  ~Lattice() {}

  // accessors
    // static and constant expression values
  static constexpr double dx() const            { return dx_; }
  static constexpr double dt() const            { return dt_; }
  static constexpr double c() const             { return dx_ / dt_; }
  static constexpr double cs() const            { return c() / sqrt(3.0); }
  static constexpr double cssq() const          { return cs() * cs(); }
  inline unsigned num_x() const                 { return nx_; }
  inline unsigned num_y() const                 { return ny_; }
  static constexpr unsigned num_k() const       { return 9; }
  inline const double* pf() const noexcept      { return spf_.get(); }
  inline double f(unsigned i, unsigned j, unsigned k) const noexcept     
                                                { return f_(i, j, k); }
  inline const double* pftemp() const noexcept  { return spftemp_.get(); }
  inline double ftemp(unsigned i, unsigned j, unsigned k) const noexcept 
                                                { return ftemp_(i, j, k); }
  inline const std::vector<std::unique_ptr<NodeDesc>>& node_descs ()
    const noexcept                            { return node_descs_; }
  inline const NodeDesc& node_desc(const unsigned i, const unsigned j) const
    { return *(node_descs_[nx_ * i + j]); }
  inline const double* pk(const unsigned k) const noexcept 
    { return (&lat_vecs_[2 * k]); }
  inline double k(const unsigned k, const unsigned c) const noexcept 
    { return *(pk(k) + c); }
    
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
  inline void stream() { stream(0, nx_ - 1, 0, nj_ - 1); }
  void stream(const std::vector<std::array<unsigned, 4>>&);

    // collide
  inline void collide_and_bound(const unsigned i, const unsigned j) 
    { node_desc(i, j).collide_and_bound(*this, i, j); }
  inline void collide_and_bound(const unsigned bi, const unsigned ei, 
                                const unsigned bj, const unsigned ej)
    { 
      for (unsigned i = bi; i <= ei; ++i)
        for (unsigned j = bj; j <= ej; ++j)
          collide_and_bound(i, j);
    }
  inline void collide_and_bound() { collide_and_bound(0, nx_ - 1, 0, nj_ - 1); }
  void collide_and_bound(const std::vector<std::array<unsigned, 4>>&);

private:
  static const double[9][2] lat_vecs_;
  static const double dx_;
  static const double dt_;
  const unsigned nx_;
  const unsigned ny_;
  std::unique_ptr<double[]> spf_;
  std::unique_ptr<double[]> spftemp_;
  std::vector<std::unique_ptr<NodeDesc>> node_descs_;

  // lattice vectors and particle distribution functions
  inline double* pf_(const unsigned i, const unsigned j) 
    { return &(spf_[i*nx_ + j]); }
  inline double& f_(const unsigned i, const unsigned j, const unsigned k)
    { return *(pf_(i, j) + k); }
  inline double* pft_(const unsigned i, const unsigned j) 
    { return &(spftemp_[i * nx_ + j]); }
  inline double& ft_(const unsigned i, const unsigned j, const unsigned k)
    { return *(pft_(i, j) + k); }
};

} // namespace d2q9

} // namespace balbm

#endif
