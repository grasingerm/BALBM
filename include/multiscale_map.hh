#ifndef MULTISCALE_MAP_HH
#define MULTISCALE_MAP_HH

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

#include <memory>
#include <armadillo>
#include "lattice.hh"

namespace balbm
{

namespace d2q9
{

//! \class AbstractMultiscaleMap
//!
//! \brief Base class for map from mesoscale to macroscale
//!
//! Maps particle distributions to macroscopic variables of interest
class AbstractMultiscaleMap
{
public:
  AbstractMultiscaleMap(const unsigned nx, const unsigned ny) 
    : nx_(nx), ny_(ny), rho_(arma::mat(nx, ny)) {}
  virtual ~AbstractMultiscaleMap()=0;
  inline double num_x() const noexcept { return nx_; }
  inline double num_y() const noexcept { return ny_; }
  inline double rho(const unsigned i, const unsigned j) const 
    { return rho_(i, j); }
  inline void map_to_macro(const Lattice& lat) { map_to_macro_(lat); }
protected:
  arma::mat rho_;
  virtual void map_to_macro_(const Lattice&);
  virtual void map_to_macro_(const Lattice&, const unsigned, const unsigned);
private:
  unsigned nx_;
  unsigned ny_;
};

//! \class DensityMultiscaleMap
//!
//! \brief Maps particle distributions to local densities
//!
//! Concrete class for density-based multiscale map
class DensityMultiscaleMap : public AbstractMultiscaleMap
{
public:
  DensityMultiscaleMap(const unsigned nx, const unsigned ny) 
    : AbstractMultiscaleMap(nx, ny) {}
  ~DensityMultiscaleMap() {}
};

class IncompFlowMultiscaleMap : public AbstractMultiscaleMap
{
public:
  IncompFlowMultiscaleMap(const unsigned nx, const unsigned ny) 
    : AbstractMultiscaleMap(nx, ny), u_(arma::cube(2, nx, ny)) {}
  ~IncompFlowMultiscaleMap() {}
  inline const arma::cube& u() const { return u_; }
  inline double u(const unsigned i, const unsigned j, const unsigned c) const 
    { return u_(c, i, j); }
private:
  void map_to_macro_(const Lattice&, const unsigned, const unsigned);
  arma::cube u_;
};

} // namespace d2q9

} // namespace balbm

#endif //MULTISCALE_MAP_HH
