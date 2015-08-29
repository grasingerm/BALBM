#ifndef __MULTISCALE_MAP_HH__
#define __MULTISCALE_MAP_HH__

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
    : _nx(nx), _ny(ny), _rho(std::unique_ptr<double[]>(new double[nx * ny])) {}
  virtual ~AbstractMultiscaleMap()=0;
  inline double num_x() const noexcept { return _nx; }
  inline double num_y() const noexcept { return _ny; }
  inline double rho(const unsigned i, const unsigned j) const 
    { return _sprho[i * _nx + j]; }
  inline void map_to_macro(const Lattice& lat) { _map_to_macro(lat); }
protected:
  inline double& _rho(const unsigned i, const unsigned j) 
    { return _sprho[i * _nx + j]; }
  virtual void _map_to_macro(const Lattice&);
  virtual void _map_to_macro(const Lattice&, const unsigned, const unsigned);
private:
  unsigned _nx;
  unsigned _ny;
  std::unique_ptr<double[]> _sprho;
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
    : AbstractMultiscaleMap(nx, ny) {}
  ~IncompFlowMultiscaleMap() {}
  inline double u(const unsigned i, const unsigned j, const unsigned c) const 
    { return _spu[2 * (i * _nx + j) + c]; }
private:
  inline double& _u(const unsigned i, const unsigned j, const unsigned c)
    { return _spu[2 * (i * _nx + j) + c]; }
  void _map_to_macro(const Lattice&, const unsigned, const unsigned);
  std::unique_ptr<double[]> _spu;
};

} // namespace d2q9

} // namespace balbm

#endif //__MULTISCALE_MAP_HH__
