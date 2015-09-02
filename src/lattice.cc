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

namespace balbm
{

namespace d2q9
{

//! \var static class member lat_vecs Vectors for each lattice direction
const double Lattice::lat_vecs_[9][2] = { {   0.0,    0.0   },
                                       {   1.0,    0.0   }, {   0.0,    1.0   },
                                       {  -1.0,    0.0   }, {   0.0,   -1.0   },
                                       {   1.0,    1.0   }, {  -1.0,    1.0   },
                                       {  -1.0,   -1.0   }, {   1.0,   -1.0   }
                                     };

//! \var static class member w_ Weights for each lattice direction
const double Lattice::w_ = { 4./9.,
                             1./9.,  1./9.,  1./9.,  1./9.,
                             1./36., 1./36., 1./36., 1./36. };

const double Lattice::dx_ = 1.0;
const double Lattice::dt_ = 1.0;

//! Copy constructor for  lattice
//!
//! \param lat Lattice to copy
//! \return Copied lattice
Lattice::Lattice(const Lattice& lat) 
  : nx_(lat.num_x()), ny_(lat.num_y()), spf_(new double[nx_ * ny_ * num_k()]),
    spftemp_(new double[nx_ * ny_ * num_k()]), node_descs_(lat.node_descs())
{
  std::copy(&lat.spf_[0], &lat.spf_[nx * ny * num_k() - 1], &spf_[0]);
  std::copy(&lat.spftemp_[0], &lat.spftemp_[nx * ny * num_k() - 1], 
            &spftemp_[0]);
}

//! Assignment for  lattice
//!
//! \param lat Lattice to assign
//! \return Copied lattice
Lattice& Lattice::operator=(const Lattice& lat) 
{
  if (this == &lat) return *this;

  if (lat.num_x() * lat.num_y() > nx_ * ny_)
  {
    //TODO: consider writing an iterator for the lattice class
    spf_.reset(new double[lat.num_x() * lat.num_y() * num_k()]);
    spftemp_.reset(new double[lat.num_x() * lat.num_y() * num_k()]);
  }

  nx_ = lat.num_x();
  ny_ = lat.num_y();
  std::copy(&lat.spf_[0], &lat.spf_[nx * ny * num_k() - 1], &spf_[0]);
  std::copy(&lat.spftemp_[0], &lat.spftemp_[nx * ny * num_k() - 1], 
            &spftemp_[0]);

  return *this;
}

//! Move constructor for  lattice
//!
//! \param lat Lattice to be moved
//! \return Moved lattice
Lattice::Lattice(Lattice&& lat)
  : nx_(lat.nx_), ny_(lat.ny_), spf_(std::move(lat.spf_)), 
    spftemp_(std::move(lat.spftemp_)), 
    node_descs_(std::move(lat.node_descs_))
{
  lat.spf_.reset(nullptr);
  lat.spftemp_.reset(nullptr);
}

//! Move assignment operator
//!
//! \param Lattice to be move assigned
//! \return Moved lattice
Lattice& Lattice::operator=(Lattice&& lat)
{
  nx_ = lat.num_x();
  ny_ = lat.num_y();
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
void stream(const std::vector<std::array<unsigned, 4>>& bounds)
{
  for (const auto& row : bounds)
    stream(row[0], row[1], row[2], row[3]);
}

} // namespace d2q9

} // namespace balbm
