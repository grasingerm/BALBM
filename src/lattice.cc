#include "lattice.hh"
#include <algorithm>
#include <array>
#include <vector>

namespace balbm
{

//! Copy constructor for D2Q9 lattice
//!
//! \param lat Lattice to copy
//! \return Copied lattice
LatticeD2Q9::LatticeD2Q9(const LatticeD2Q9& lat) 
  : nx(lat.nx), ny(lat.ny), f(new double[nx * ny * num_k()]),
    ftemp(new double[nx * ny * num_k()]), node_descs(lat.node_descs)
{
  std::copy(&lat.f[0], &lat.f[nx * ny * num_k() - 1], &f[0]);
  std::copy(&lat.ftemp[0], &lat.ftemp[nx * ny * num_k() - 1], &ftemp[0]);
}

//! Assignment for D2Q9 lattice
//!
//! \param lat Lattice to assign
//! \return Copied lattice
LatticeD2Q9& LatticeD2Q9::operator=(const LatticeD2Q9& lat) 
  : nx(lat.nx), ny(lat.ny), f(new double[nx * ny * num_k()]),
    ftemp(new double[nx * ny * num_k()]), node_descs(lat.node_descs)
{
  if (this == &lat) return *this;

  if (f.nx * f.ny > nx * ny)
  {
    delete[] f;
    delete[] ftemp;
    f = new double[nx * ny * num_k()];
    ftemp = new double[nx * ny * num_k()];
  }

  nx = lat.nx;
  ny = lat.ny;
  std::copy(&lat.f[0], &lat.f[nx * ny * num_k() - 1], &f[0]);
  std::copy(&lat.ftemp[0], &lat.ftemp[nx * ny * num_k() - 1], &ftemp[0]);

  return *this;
}

//! Move constructor for D2Q9 lattice
//!
//! \param lat Lattice to be moved
//! \return Moved lattice
LatticeD2Q9::LatticeD2Q9(LatticeD2Q9&& lat)
  : nx(lat.nx), ny(lat.ny), f(lat.f), ftemp(lat.ftemp), 
    node_descs(std::move(lat.node_descs))
{
  lat.f = nullptr;
  lat.ftemp = nullptr;
}

//! Move assignment operator
//!
//! \param Lattice to be move assigned
//! \return Moved lattice
LatticeD2Q9& LatticeD2Q9::operator=(LatticeD2Q9&& lat)
{
  nx = lat.nx;
  ny = lat.ny;
  f = lat.f;
  ftemp = lat.ftemp;
  node_descs = std::move(lat.node_descs);

  lat.f = nullptr;
  lat.ftemp = nullptr;

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
