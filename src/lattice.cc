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

#include "lattice.hh"
#include <algorithm>
#include <array>
#include <vector>

namespace balbm {

namespace d2q9 {

//! \var static class member lat_vecs_ Vectors for each lattice direction
const double Lattice::lat_vecs_[Lattice::nk_][2] = {{0.0, 0.0},
                                                    {1.0, 0.0},
                                                    {0.0, 1.0},
                                                    {-1.0, 0.0},
                                                    {0.0, -1.0},
                                                    {1.0, 1.0},
                                                    {-1.0, 1.0},
                                                    {-1.0, -1.0},
                                                    {1.0, -1.0}};

//! \var static class member w_ Weights for each lattice direction
const double Lattice::w_[] = {4. / 9.,  1. / 9.,  1. / 9.,  1. / 9., 1. / 9.,
                              1. / 36., 1. / 36., 1. / 36., 1. / 36.};

const double Lattice::dx_ = 1.0;
const double Lattice::dt_ = 1.0;

//! Copy constructor for  lattice
//!
//! \param lat Lattice to copy
//! \return Copied lattice
Lattice::Lattice(const Lattice &lat)
    : ni_(lat.num_i()), nj_(lat.num_j()), spf_(new double[ni_ * nj_ * num_k()]),
      spftemp_(new double[ni_ * nj_ * num_k()]), node_descs_(lat.node_descs()) {
  std::copy(&lat.spf_[0], &lat.spf_[ni_ * nj_ * num_k() - 1], &spf_[0]);
  std::copy(&lat.spftemp_[0], &lat.spftemp_[ni_ * nj_ * num_k() - 1],
            &spftemp_[0]);
}

//! Assignment for  lattice
//!
//! \param lat Lattice to assign
//! \return Copied lattice
Lattice &Lattice::operator=(const Lattice &lat) {
  if (this == &lat)
    return *this;

  if (lat.num_i() * lat.num_j() > ni_ * nj_) {
    // TODO: consider writing an iterator for the lattice class
    spf_.reset(new double[lat.num_i() * lat.num_j() * num_k()]);
    spftemp_.reset(new double[lat.num_i() * lat.num_j() * num_k()]);
  }

  ni_ = lat.num_i();
  nj_ = lat.num_j();
  std::copy(&lat.spf_[0], &lat.spf_[ni_ * nj_ * num_k() - 1], &spf_[0]);
  std::copy(&lat.spftemp_[0], &lat.spftemp_[ni_ * nj_ * num_k() - 1],
            &spftemp_[0]);

  return *this;
}

//! Move constructor for  lattice
//!
//! \param lat Lattice to be moved
//! \return Moved lattice
Lattice::Lattice(Lattice &&lat)
    : ni_(lat.ni_), nj_(lat.nj_), spf_(std::move(lat.spf_)),
      spftemp_(std::move(lat.spftemp_)),
      node_descs_(std::move(lat.node_descs_)) {
  lat.spf_.reset(nullptr);
  lat.spftemp_.reset(nullptr);
}

//! Move assignment operator
//!
//! \param Lattice to be move assigned
//! \return Moved lattice
Lattice &Lattice::operator=(Lattice &&lat) {
  ni_ = lat.num_i();
  nj_ = lat.num_j();
  spf_ = std::move(lat.spf_);
  spftemp_ = std::move(lat.spftemp_);
  node_descs_ = std::move(lat.node_descs_);

  lat.f.reset(nullptr);
  lat.ftemp.reset(nullptr);

  return *this;
}

//! Stream particle distribution functions
//!
//! \param bounds Vector of arrays of {begin_i, end_i, begin_j, end_j}
void stream(const std::vector<std::array<unsigned, 4>> &bounds) {
  for (const auto &row : bounds)
    stream(row[0], row[1], row[2], row[3]);
}

//! Perform bounds checking
//!
//! \param i Index in the y-direction
//! \param j Index in the x-direction
//! \return true if in bounds, false if out of bounds
//! \throw out_of_range
bool Lattice::check_bounds(const unsigned i, const unsigned j) const
    throw(std::out_of_range) {
  bool is_in_bounds = in_bounds(i, j);
  if (is_in_bounds) {
    std::ostringstream oss;
    oss << "Ill-defined boundaries. Particles streamed out of bounds from "
        << "node (" << i << " ," << j << ") to (" << i_next << ", " << j_next
        << "). Check boundary conditions.";
    throw std::out_of_range(oss.str());
  }

  return is_in_bounds;
}

//! Initialize domain to equilibrium based on a reference density
//!
//! \param rho Reference density
void Lattice::init_f_(const double rho) {
  const unsigned ni = num_i();
  const unsigned nj = num_j();

  // should this loop be multithreaded/parallelized?
  for (unsigned i = 0; i < ni; ++i)
    for (unsigned j = 0; j < nj; ++j) {
      f_(i, j, 0) = w(0) * rho;
      f_(i, j, 1) = w(1) * rho;
      f_(i, j, 2) = w(2) * rho;
      f_(i, j, 3) = w(3) * rho;
      f_(i, j, 4) = w(4) * rho;
      f_(i, j, 5) = w(5) * rho;
      f_(i, j, 6) = w(6) * rho;
      f_(i, j, 7) = w(7) * rho;
      f_(i, j, 8) = w(8) * rho;
    }
  /*  for (unsigned k = 0; k < num_k(); ++k) // TODO: should I unroll this
     inner?
      f_(i, j, k) = w(k) * rho; */
}

} // namespace d2q9

} // namespace balbm
