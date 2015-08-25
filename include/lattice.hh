#ifndef __LATTICE_HH__
#define __LATTICE_HH__

#include "balbm_config.hh"
#include <vector>
#include <memory>
#include <algorithm>

namespace balbm
{

class LatticeD2Q9
{
public:
  LatticeD2Q9() : nx(0), ny(0), f(nullptr), ftemp(nullptr) {}
  LatticeD2Q9(unsigned nx, unsigned ny) : nx(nx), ny(ny),
    f(new double[nx * ny * num_k()]), ftemp(new double[nx * ny * num_k()]),
    node_descs(std::vector<std::unique_ptr<NodeDesc>>(nx * ny * num_k())) {}
  LatticeD2Q9(const LatticeD2Q9&);
  LatticeD2Q9& operator=(const LatticeD2Q9&);
  LatticeD2Q9(LatticeD2Q9&&);
  LatticeD2Q9& operator=(LatticeD2Q9&&);
  ~LatticeD2Q9() { delete[] f; delete[] ftemp; }
  inline unsigned num_x() const { return nx; }
  inline unsigned num_y() const { return ny; }
  static constexpr unsigned num_k() const { return 9; }
  inline double* pk(unsigned k) const { return (lat_vecs + 2 * k); }
  inline double k(unsigned c, unsigned k) const 
    { return *(pk(k) + c); }
  inline double* pf(unsigned i, unsigned j) const { return (f + i*ny + j); }
  inline double f(unsigned i, unsigned j, unsigned k) const
    { return *(pf(i, j) + k); }
  inline double* pft(unsigned i, unsigned j) const 
    { return (ftemp + i*ny + j); }
  inline double ft(unsigned i, unsigned j, unsigned k) const
    { return *(pft(i, j) + k); }
  inline NodeDesc& node_desc(unsigned i, unsigned j) const
    { return node_descs[ny * i + j]; }
private:
  static const double[9][2] lat_vecs;
  const unsigned nx;
  const unsigned ny;
  double* f;
  double* ftemp;
  node_descs vector<std::unique_ptr<NodeDesc>>;
  inline void finalize_step() { std::swap(f, ftemp); }
};

#endif
