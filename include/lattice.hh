#ifndef __LATTICE_HH__
#define __LATTICE_HH__

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace balbm
{

//! \class LatticeD2Q9
//!
//! \brief D2Q9 lattice for the lattice Boltzmann method
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
class LatticeD2Q9
{
public:
  // constructors and assignment
  LatticeD2Q9() : _nx(0), _ny(0), _f(nullptr), _ftemp(nullptr) {}
  LatticeD2Q9(unsigned nx, unsigned ny) : _nx(nx), _ny(ny),
    f(new double[nx * ny * num_k()]), ftemp(new double[nx * ny * num_k()]),
    node_descs(std::vector<std::unique_ptr<NodeDescD2Q9>>(nx * ny * num_k())) {}
  LatticeD2Q9(const LatticeD2Q9&);
  LatticeD2Q9& operator=(const LatticeD2Q9&);
  LatticeD2Q9(LatticeD2Q9&&);
  LatticeD2Q9& operator=(LatticeD2Q9&&);
  ~LatticeD2Q9() { delete[] _f; delete[] _ftemp; }

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
  inline double* pk(const unsigned k) const { return (_lat_vecs + 2 * k); }
  inline double k(const unsigned c, const unsigned k) const 
    { return *(pk(k) + c); }
  inline double* pf(const unsigned i, const unsigned j) const 
    { return (_f + i*_nx + j); }
  inline double f(const unsigned i, const unsigned j, const unsigned k) const
    { return *(pf(i, j) + k); }
  inline double* pft(const unsigned i, const unsigned j) const 
    { return (ftemp + i*_ny + j); }
  inline double ft(const unsigned i, const unsigned j, const unsigned k) const
    { return *(pft(i, j) + k); }
  inline NodeDescD2Q9& node_desc(const unsigned i, const unsigned j) const
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
  inline void collide(const unsigned i, const unsigned j) 
    { node_desc(i, j).collide(*this, i, j); }
  inline void collide(const unsigned bi, const unsigned ei, const unsigned bj, 
                      const unsigned ej)
    { 
      for (unsigned i = bi; i <= ei; ++i)
        for (unsigned j = bj; j <= ej; ++j)
          collide(i, j);
    }
  inline void collide() { collide(0, _nx - 1, 0, _nj - 1); }
  void collide(const std::vector<std::array<unsigned, 4>>&);

private:
  static const double[9][2] _lat_vecs;
  static const double _dx;
  static const double _dt;
  const unsigned _nx;
  const unsigned _ny;
  double* _f;
  double* _ftemp;
  std::vector<std::unique_ptr<NodeDescD2Q9>> _node_descs;
  inline void finalize_step() { std::swap(_f, _ftemp); }
};

} // namespace balbm

#endif
