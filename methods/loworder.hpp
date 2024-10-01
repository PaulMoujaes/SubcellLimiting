#ifndef SUBCELL_LOWORDER
#define SUBCELL_LOWORDER

#include "feevol.hpp"


class LowOrderScheme : public FE_Evolution
{
public:
   LowOrderScheme(ParFiniteElementSpace &fes_,
                  const Vector &lumpedmassmatrix_, FunctionCoefficient &inflow,
                  VectorFunctionCoefficient &velocity, ParBilinearForm &M);

   virtual void Mult(const Vector &x, Vector &y) const override;

   virtual ~LowOrderScheme();
};

#endif