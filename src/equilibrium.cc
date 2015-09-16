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

#include "equilibrium.hh"
#include "lattice.hh"
#include <armadillo>

namespace balbm {

namespace d2q9 {

//! Virtual destructor definition
AbstractIncompFlowEqFunct::~AbstractIncompFlowEqFunct() {}

//! Equilibrium distribution function for incompressible flow
//!
//! \param lat Lattice
//! \param rho Density at the lattice node
//! \param u Macroscopic velocity vector
//! \param k Index of lattice direction
//! \return Equilibrium particle distribution
double IncompFlowEqFunct::f_(const Lattice &lat, const double rho,
                             const arma::vec &u, const unsigned k) const {
  const arma::vec ck(const_cast<double *>(lat.pc(k)), 2, false, true);
  const double ckdotu = arma::dot(ck, u);
  const double cssq = lat.cssq();

  return rho * lat.w(k) *
         (1.0 + ckdotu / cssq + 0.5 * (ckdotu * ckdotu) / (cssq * cssq) -
          0.5 * dot(u, u) / cssq);
}

//! Equilibrium distribution function for incompressible flow
//!
//! \param lat Lattice
//! \param rho Density at the lattice node
//! \param u Macroscopic velocity vector
//! \param k Index of lattice direction
//! \return Equilibrium particle distribution
double IncompFlowHLEqFunct::f_(const Lattice &lat, const double rho,
                               const arma::vec &u, const unsigned k) const {
  const arma::vec ck(const_cast<double *>(lat.pc(k)), 2, false, true);
  const double ckdotu = arma::dot(ck, u);
  const double cssq = lat.cssq();

  return (lat.w(k) * (rho + rho_o_ +
                      (ckdotu / cssq + 0.5 * (ckdotu * ckdotu) / (cssq * cssq) -
                       0.5 * dot(u, u) / cssq)));
}

} // namespace d2q9

} // namespace balbm
