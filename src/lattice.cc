#include "lattice.hh"
#include <algorithm>
#include <array>
#include <vector>

namespace balbm
{

namespace d2q9
{

//! \var static class member lat_vecs Lattice vectors for 
const double Lattice::lat_vecs[9][2] = { {   0.0,    0.0   },
                                       {   1.0,    0.0   }, {   0.0,    1.0   },
                                       {  -1.0,    0.0   }, {   0.0,   -1.0   },
                                       {   1.0,    1.0   }, {  -1.0,    1.0   },
                                       {  -1.0,   -1.0   }, {   1.0,   -1.0   }
                                     };

const double Lattice::_dx = 1.0;
const double Lattice::_dt = 1.0;

//! Copy constructor for  lattice
//!
//! \param lat Lattice to copy
//! \return Copied lattice
Lattice::Lattice(const Lattice& lat) 
  : nx(lat.nx), ny(lat.ny), f(new double[nx * ny * num_k()]),
    ftemp(new double[nx * ny * num_k()]), node_descs(lat.node_descs)
{
  std::copy(&lat.f[0], &lat.f[nx * ny * num_k() - 1], &f[0]);
  std::copy(&lat.ftemp[0], &lat.ftemp[nx * ny * num_k() - 1], &ftemp[0]);
}

//! Assignment for  lattice
//!
//! \param lat Lattice to assign
//! \return Copied lattice
Lattice& Lattice::operator=(const Lattice& lat) 
  : nx(lat.nx), ny(lat.ny), f(new double[nx * ny * num_k()]),
    ftemp(new double[nx * ny * num_k()]), node_descs(lat.node_descs)
{
  if (this == &lat) return *this;

  if (f.nx * f.ny > nx * ny)
  {
    f = std::make_unique(new double[nx * ny * num_k()]);
    ftemp = std::make_unique(new double[nx * ny * num_k()]);
  }

  nx = lat.nx;
  ny = lat.ny;
  std::copy(&lat.f[0], &lat.f[nx * ny * num_k() - 1], &f[0]);
  std::copy(&lat.ftemp[0], &lat.ftemp[nx * ny * num_k() - 1], &ftemp[0]);

  return *this;
}

//! Move constructor for  lattice
//!
//! \param lat Lattice to be moved
//! \return Moved lattice
Lattice::Lattice(Lattice&& lat)
  : nx(lat.nx), ny(lat.ny), f(std::move(lat.f)), ftemp(std::move(lat.ftemp)), 
    node_descs(std::move(lat.node_descs))
{
  lat.f.reset(nullptr);
  lat.ftemp.reset(nullptr);
}

//! Move assignment operator
//!
//! \param Lattice to be move assigned
//! \return Moved lattice
Lattice& Lattice::operator=(Lattice&& lat)
{
  nx = lat.nx;
  ny = lat.ny;
  f = std::move(lat.f);
  ftemp = std::move(lat.ftemp);
  node_descs = std::move(lat.node_descs);

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
