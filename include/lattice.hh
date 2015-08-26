#ifndef __LATTICE_HH__
#define __LATTICE_HH__

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace balbm
{

namespace d2q9
{

//! \class Lattice
//!
//! \brief  lattice for the lattice Boltzmann method
//!
//! 6     2     5
//!   \   |   /
//!    \  |  /
//!     \ | /
//! 3 --- 0 --- 1
//!     / | \
//!    /  |  \
//!   /   |   \
//! 7     4     8
//!
class Lattice
{
public:
  // constructors and assignment
  Lattice() : _nx(0), _ny(0), _f(nullptr), _ftemp(nullptr) {}
  Lattice(unsigned nx, unsigned ny) : _nx(nx), _ny(ny),
    _f(std::make_unique(new double[nx * ny * num_k()])), 
    _ftemp(std::make_unique(new double[nx * ny * num_k()])),
    _node_descs(std::vector<std::unique_ptr<NodeDesc>>(nx * ny * num_k())) {}
  Lattice(const Lattice&);
  Lattice& operator=(const Lattice&);
  Lattice(Lattice&&);
  Lattice& operator=(Lattice&&);
  ~Lattice() {}

  // accessors
    // static and constant expression values
  static constexpr double dx() const { return _dx; }
  static constexpr double dt() const { return _dt; }
  static constexpr double c() { return _dx / _dt; }
  static constexpr double cs() { return c() / sqrt(3.0); }
  static constexpr double cssq() { return cs() * cs(); }
  inline unsigned num_x() const { return _nx; }
  inline unsigned num_y() const { return _ny; }
  static constexpr unsigned num_k() const { return 9; }
    // lattice vectors and particle distribution functions
  inline double* pk(const unsigned k) const { return (&_lat_vecs[2 * k]); }
  inline double k(const unsigned c, const unsigned k) const 
    { return *(pk(k) + c); }
  inline double* pf(const unsigned i, const unsigned j) const 
    { return &(_f[i*_nx + j]); }
  inline double f(const unsigned i, const unsigned j, const unsigned k) const
    { return *(pf(i, j) + k); }
  inline double* pft(const unsigned i, const unsigned j) const 
    { return &(ftemp[i*_nx + j]); }
  inline double ft(const unsigned i, const unsigned j, const unsigned k) const
    { return *(pft(i, j) + k); }
  inline NodeDesc& node_desc(const unsigned i, const unsigned j) const
    { return *(_node_descs[_nx * i + j]); }

  // mutators
    // stream
  inline void stream(const unsigned i, const unsigned j) 
    { node_desc(i, j).stream(*this, i, j); }
  inline void stream(const unsigned bi, const unsigned ei, const unsigned bj, 
                     const unsigned ej)
    { 
      for (unsigned i = bi; i <= ei; ++i)
        for (unsigned j = bj; j <= ej; ++j)
          stream(i, j);
    }
  inline void stream() { stream(0, _nx - 1, 0, _nj - 1); }
  void stream(const std::vector<std::array<unsigned, 4>>&);

    // collide
  inline void collide_and_bound(const unsigned i, const unsigned j) 
    { node_desc(i, j).collide_and_bound(*this, i, j); }
  inline void collide_and_bound(const unsigned bi, const unsigned ei, 
                                const unsigned bj, const unsigned ej)
    { 
      for (unsigned i = bi; i <= ei; ++i)
        for (unsigned j = bj; j <= ej; ++j)
          collide_and_bound(i, j);
    }
  inline void collide_and_bound() { collide_and_bound(0, _nx - 1, 0, _nj - 1); }
  void collide_and_bound(const std::vector<std::array<unsigned, 4>>&);

private:
  static const double[9][2] _lat_vecs;
  static const double _dx;
  static const double _dt;
  const unsigned _nx;
  const unsigned _ny;
  std::unique_ptr<double[]> _f;
  std::unique_ptr<double[]> _ftemp;
  std::vector<std::unique_ptr<NodeDesc>> _node_descs;
  inline void finalize_step() { std::swap(_f, _ftemp); }
};

} // namespace d2q9

} // namespace balbm

#endif
