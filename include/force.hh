#ifndef FORCE_HH
#define FORCE_HH

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
#include <armadillo>
#include <array>

namespace balbm {

namespace d2q9 {

class Lattice;

//! \class AbstractForce
//!
//! \brief Abstract base class for force implementations
//!
//! Provide polymorphic behavior for forces and force implementations
class AbstractForce {
public:
  virtual ~AbstractForce() = 0;
  AbstractForce(double *F) : F_(arma::vec::fixed<2>(F)) {}
  inline arma::vec::fixed<2> u_trans(const Lattice &lat,
                                     const arma::vec::fixed<2> &u) const {
    return u_trans_(lat, u);
  }
  inline double f_col(const Lattice &lat, const double omega,
                      const arma::vec::fixed<2> &u, const unsigned k) const {
    return f_col_(lat, omega, u, k);
  }
  inline const arma::vec::fixed<2> &F() const { return F_; }

private:
  arma::vec::fixed<2> F_;
  virtual arma::vec::fixed<2> u_trans_(const Lattice &,
                                       const arma::vec::fixed<2> &) const = 0;
  virtual double f_col_(const Lattice &, const double omega,
                        const arma::vec::fixed<2> &,
                        const unsigned k) const = 0;
};

//! \class SukopThorneForce
//!
//! \brief Force implementation from Sukop and Thorne 2005
//!
//! One of the forcing implementations found in Lattice Boltzmann Method for
//! Geoscientists and Engineers by Sukop and Thorne 2005
class SukopThorneForce : public AbstractForce {
public:
  ~SukopThorneForce(){};
  SukopThorneForce(double *F) : AbstractForce(F) {}

private:
  arma::vec::fixed<2> u_trans_(const Lattice &,
                               const arma::vec::fixed<2> &) const;
  double f_col_(const Lattice &, const double, const arma::vec::fixed<2> &,
                const unsigned) const;
};

//! \class GuoForce
//!
//! \brief Force implementation from Guo 2002
//!
//! One of the forcing implementations found in Discrete lattice effects on the
//! forcing termin the lattice Boltzmann method Guo et. al. 2002
class GuoForce : public AbstractForce {
public:
  ~GuoForce(){};
  GuoForce(double *F) : AbstractForce(F) {}

private:
  arma::vec::fixed<2> u_trans_(const Lattice &,
                               const arma::vec::fixed<2> &) const;
  double f_col_(const Lattice &, const double, const arma::vec::fixed<2> &,
                const unsigned) const;
};

} // namespace d2q9

} // namespace balbm

#endif // FORCE_HH
