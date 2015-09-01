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
#include <armadillo>

namespace balbm
{

namespace d2q9
{

//! \var static class member lat_vecs Lattice vectors for 
const arma::mat::fixed<2, 9> Lattice::lat_vecs_ ({ {   0.0,    0.0   },
                                       {   1.0,    0.0   }, {   0.0,    1.0   },
                                       {  -1.0,    0.0   }, {   0.0,   -1.0   },
                                       {   1.0,    1.0   }, {  -1.0,    1.0   },
                                       {  -1.0,   -1.0   }, {   1.0,   -1.0   }
                                                });

const double Lattice::dx_ = 1.0;
const double Lattice::dt_ = 1.0;

//! Copy constructor for  lattice
//!
//! \param lat Lattice to copy
//! \return Copied lattice
Lattice::Lattice(const Lattice& lat) 
  : nx_(lat.num_x()), ny_(lat.num_y()), f_(lat.f()),
    ftemp_(lat.ftemp()), node_descs_(lat.node_descs())
{}

//! Assignment for  lattice
//!
//! \param lat Lattice to assign
//! \return Copied lattice
Lattice& Lattice::operator=(const Lattice& lat) 
{
  if (this == &lat) return *this;

  nx_ = lat.num_x();
  ny_ = lat.num_y();

  // resize and copy node descriptors
  node_descs_.resize(_nx * _ny);
  for (unsigned j = 0; j < _ny; ++j)
    for (unsigned i = 0; i < _nx; ++i)
      node_desc(i, j) = lat.node_desc(i, j);

  f_ = lat.f();
  ftemp_ = lat.ftemp();

  return *this;
}

//! Move constructor for  lattice
//!
//! \param lat Lattice to be moved
//! \return Moved lattice
Lattice::Lattice(Lattice&& lat)
  : nx_(lat.nx_), ny_(lat.ny_), f_(std::move(lat.f_)), 
    ftemp_(std::move(lat.ftemp_)), node_descs_(std::move(lat.node_descs_))
{}

//! Move assignment operator
//!
//! \param Lattice to be move assigned
//! \return Moved lattice
Lattice& Lattice::operator=(Lattice&& lat)
{
  nx_ = lat.nx_;
  ny_ = lat.ny_;
  f_ = std::move(lat.f_);
  ftemp_ = std::move(lat.ftemp_);
  node_descs_ = std::move(lat.node_descs_);

  return *this;
}

//! Stream particle distribution functions
//!
//! \param bounds Vector of arrays of {begin_i, end_i, begin_j, end_j}
void Lattice::stream(const std::vector<std::array<unsigned, 4>>& bounds)
{
  for (const auto& row : bounds)
    stream(row[0], row[1], row[2], row[3]);
}

} // namespace d2q9

} // namespace balbm
