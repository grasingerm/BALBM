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

#include "force.hh"
#include <armadillo>

namespace balbm
{

namespace d2q9
{

//! Virtual destructor for AbstractForce base class
AbstractForce::~AbstractForce() {}

//! Transforms macroscopic velocity vector to simulate external forces
//!
//! \param lat Lattice
//! \param u Macroscopic velocity vector
//! \return Transformed velocity vector
arma::vec::fixed<2> SukopThorneForce::u_trans
  (const Lattice& lat, const arma::vec& u) const
{
  return u;
}

//! Value to add to particle distributions to simulate effect of external forces
//!
//! \param lat Lattice
//! \param omega Collision frequency
//! \param u Macroscopic velocity vector
//! \param k Index of lattice direction
double SukopThorneForce::f_col_(const Lattice& lat, const double omega, 
                                const arma::vec& u, const unsigned k)
{
  const arma::vec::fixed<2> ck(lat.cp(k), 2, false, true);
  return lat.w(k) * lat.dt() / lat.cssq() * arma::dot(F(), ck);
}

//! Transforms macroscopic velocity vector to simulate external forces
//!
//! \param lat Lattice
//! \param u Macroscopic velocity vector
//! \return Transformed velocity vector
arma::vec::fixed<2> GuoForce::u_trans
  (const Lattice& lat, const arma::vec& u) const
{
  return u + (lat.dt() / 2.0 * F());
}

//! Value to add to particle distributions to simulate effect of external forces
//!
//! \param lat Lattice
//! \param omega Collision frequency
//! \param u Macroscopic velocity vector
//! \param k Index of lattice direction
double GuoForce::f_col_(const Lattice& lat, const double omega, 
                        const arma::vec& u, const unsigned k)
{
  const arma::vec::fixed<2> ck(lat.cp(k), 2, false, true);
  return ((1 - 0.5 * omega) * lat.w(k) 
           * arma::dot((ck - u)/(lat.cssq()) 
                        - arma::dot(ck, u) / (lat.cssq()*lat.cssq()) * ck,
                       F()));
}

} // namespace d2q9

} // namespace balbm
