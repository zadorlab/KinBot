// Compile against the SAME MESS headers and LAPACK library as the executable.
// Checks LP64/ILP64 status and pivot-array compatibility, not chemical accuracy.
#include "libmess/lapack.hh"
#include "libmess/lapack.h"
#include <cmath>
#include <cstdint>
#include <iostream>

int main() {
  using I = Lapack::int_t;
  static_assert(sizeof(I)==4 || sizeof(I)==8, "Expected LP64 or ILP64 integers");
  // Padded arguments also make the reverse mismatch safe to diagnose: a
  // 64-bit library can read/write eight bytes even with a 32-bit header.
  // Zero padding represents the same dimension for either little-endian ABI.
  const std::uint32_t endian=1;
  if (*reinterpret_cast<const unsigned char*>(&endian)!=1) {
    std::cerr << "This boundary probe requires little-endian storage.\n";
    return 2;
  }
  alignas(8) I n[8/sizeof(I)]={}, nrhs[8/sizeof(I)]={}, info[8/sizeof(I)]={};
  n[0]=2; nrhs[0]=1;
  double packed[3]={4.,1.,3.};
  // Nonzero upper bits expose a 32-bit library writing a 64-bit status.
  info[0] = sizeof(I)==8 ? I(UINT64_C(0x100000001)) : I(1);
  if(sizeof(I)==4) info[1]=0x5a5a5a5a;
  dpptrf_('U', n[0], packed, info[0]);
  bool ok = info[0]==0 && (sizeof(I)==8 || info[1]==0x5a5a5a5a) && std::abs(packed[0]-2.)<1e-14 &&
            std::abs(packed[1]-.5)<1e-14 && std::abs(packed[2]-std::sqrt(2.75))<1e-14;
  std::cout << "integer_bytes=" << sizeof(I) << " cholesky_info=" << info[0]
            << " cholesky_pass=" << ok << "\n";
  // If the status interface is wrong, do not exercise a pivoted solve.
  if(!ok) return 1;
  double a[4]={0.,2.,1.,3.}, rhs[2]={1.,5.};
  I piv[2]={-1,-1}; info[0]=1;
  dgesv_(n[0],nrhs[0],a,n[0],piv,rhs,n[0],info[0]);
  ok=info[0]==0 && piv[0]==2 && piv[1]==2 &&
      std::abs(rhs[0]-1.)<1e-14 && std::abs(rhs[1]-1.)<1e-14;
  std::cout << "solve_info=" << info[0] << " pivots=" << piv[0] << "," << piv[1]
            << " solution=" << rhs[0] << "," << rhs[1] << " solve_pass=" << ok << "\n";
  return ok?0:1;
}
