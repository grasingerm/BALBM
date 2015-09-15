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

#include "balbm_config.hh"
#include "constitutive.hh"
#include "lattice.hh"
#include "multiscale_map.hh"
#include <armadillo>

namespace balbm {

namespace d2q9 {

//! Virtual destructor definition
AbstractConstitutiveEq::~AbstractConstitutiveEq() {}

//! Constitutive equation for a Newtonian fluid
//!
//! \param lat Lattice
//! \param mmap Multiscale map of macroscopic variables
//! \param fneq Non-equilibrium particle distribution
//! \param i Index in x-direction
//! \param j Index in y-direction
//! \return Kinematic viscosity
double NewtonianConstitutiveEq::mu_(const Lattice &lat,
                                    const AbstractMultiscaleMap &mmap,
                                    const arma::vec &fneq, const unsigned i,
                                    const unsigned j) {
  return cmu_;
}

} // d2q9

} // balbm
