#ifndef CONSTITUTIVE_HH
#define CONSTITUTIVE_HH

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
#include <armadillo>
#include <limits>

namespace balbm {

namespace d2q9 {

//! \class AbstractConstitutiveEq
//!
//! \brief Base class for constitutive equations for viscosity
//!
//! Maps particle distributions to macroscopic strain rate to viscosity
class AbstractConstitutiveEq {
public:
  virtual ~AbstractConstitutiveEq() = 0;
  inline double mu(const Lattice &lat, const AbstractMultiscaleMap &mmap,
                   const arma::vec &fneq, const unsigned i,
                   const unsigned j) const {
    return mu_(lat, mmap, fneq, i, j);
  }

private:
  virtual double mu_(const Lattice &, const AbstractMultiscaleMap &,
                     const arma::vec &, const unsigned,
                     const unsigned) const = 0;
};

//! \class NewtonianConstitutiveEq
//!
//! \brief Class for constant viscosity
//!
//! Newtonian constitutive equation
class NewtonianConstitutiveEq : public AbstractConstitutiveEq {
public:
  ~NewtonianConstitutiveEq() {}
  NewtonianConstitutiveEq(const double mu) : cmu_(mu) {}

private:
  const double cmu_;
  double mu_(const Lattice &, const AbstractMultiscaleMap &, const arma::vec &,
             const unsigned, const unsigned) const;
};

//! \class BinghamConstitutiveEq
//!
//! \brief Class for constant viscosity
//!
//! Newtonian constitutive equation
// TODO: should this be broken down further into strain rate calc and solver?
class BinghamConstitutiveEq {
public:
  ~BinghamConstitutiveEq() {}
  BinghamConstitutiveEq(
      const double mu_p, const double tau_y, const double m,
      const double gamma_min = std::numeric_limits<double>::epsilon())
      : mu_p_(mu_p), tau_y_(tau_y), m_(m), gamma_min_(gamma_min) {}

private:
  const double mu_p_;
  const double tau_y_;
  const double m_;
  const double gamma_min_;
  double mu_(const Lattice &, const AbstractMultiscaleMap &, const arma::vec &,
             const unsigned, const unsigned);
};

} // namespace d2q9

} // namespace balbm

#endif // CONSTITUTIVE_HH
