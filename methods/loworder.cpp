#include "loworder.hpp"

// Implementation of Class LowOrderScheme
LowOrderScheme::LowOrderScheme(ParFiniteElementSpace &fes_,
                               const Vector &lumpedmassmatrix_, FunctionCoefficient &inflow,
                               VectorFunctionCoefficient &velocity, ParBilinearForm &M):
   FE_Evolution(fes_, lumpedmassmatrix_, inflow, velocity, M)
{ }

void LowOrderScheme::Mult(const Vector &x, Vector &y) const
{
   ComputeLOTimeDerivatives(x, y);
}

LowOrderScheme::~LowOrderScheme()
{ }
