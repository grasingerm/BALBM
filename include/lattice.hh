#ifndef __LATTICE_HH__
#define __LATTICE_HH__

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

namespace balbm
{

class LatticeD2Q9
{
public:
  // constructors and assignment
  LatticeD2Q9() : nx(0), ny(0), f(nullptr), ftemp(nullptr) {}
  LatticeD2Q9(unsigned nx, unsigned ny) : nx(nx), ny(ny),
    f(new double[nx * ny * num_k()]), ftemp(new double[nx * ny * num_k()]),
    node_descs(std::vector<std::unique_ptr<NodeDesc>>(nx * ny * num_k())) {}
  LatticeD2Q9(const LatticeD2Q9&);
  LatticeD2Q9& operator=(const LatticeD2Q9&);
  LatticeD2Q9(LatticeD2Q9&&);
  LatticeD2Q9& operator=(LatticeD2Q9&&);
  ~LatticeD2Q9() { delete[] f; delete[] ftemp; }

  // accessors
  inline unsigned num_x() const { return nx; }
  inline unsigned num_y() const { return ny; }
  static constexpr unsigned num_k() const { return 9; }
  inline double* pk(const unsigned k) const { return (lat_vecs + 2 * k); }
  inline double k(const unsigned c, const unsigned k) const 
    { return *(pk(k) + c); }
  inline double* pf(const unsigned i, const unsigned j) const 
    { return (f + i*ny + j); }
  inline double f(const unsigned i, const unsigned j, const unsigned k) const
    { return *(pf(i, j) + k); }
  inline double* pft(const unsigned i, const unsigned j) const 
    { return (ftemp + i*ny + j); }
  inline double ft(const unsigned i, const unsigned j, const unsigned k) const
    { return *(pft(i, j) + k); }
  inline NodeDesc& node_desc(const unsigned i, const unsigned j) const
    { return *(node_descs[ny * i + j]); }

  // mutators
  inline void stream(const unsigned i, const unsigned j) 
    { node_desc(i,j).stream(i,j,f,ftemp); }
  inline void stream(const unsigned bi, const unsigned ei, const unsigned bj, 
                     const unsigned ej)
    { 
      for (unsigned i = bi; i <= ei; ++i)
        for (unsigned j = bj; j <= ej; ++j)
          stream(i, j);
    }
  inline void stream() { stream(0, nx - 1, 0, nj - 1); }
  void stream(const std::vector<std::array<unsigned, 4>>&);

private:
  static const double[9][2] lat_vecs;
  const unsigned nx;
  const unsigned ny;
  double* f;
  double* ftemp;
  node_descs vector<std::unique_ptr<NodeDesc>>;
  inline void finalize_step() { std::swap(f, ftemp); }
};

}

#endif
