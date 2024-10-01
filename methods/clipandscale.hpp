#ifndef SUBCELL_CLIPANDSCALE
#define SUBCELL_CLIPANDSCALE

#include "feevol.hpp"

class ClipAndScale : public FE_Evolution
{
private:
   mutable Array<double> umin, umax;
   mutable Vector udot;

   virtual void ComputeBounds(const Vector &u, Array<double> &u_min,
                              Array<double> &u_max) const;

public:
   ClipAndScale(ParFiniteElementSpace &fes_,
                const Vector &lumpedmassmatrix_, FunctionCoefficient &inflow,
                VectorFunctionCoefficient &velocity, ParBilinearForm &M);

   virtual void Mult(const Vector &x, Vector &y) const override;

   virtual ~ClipAndScale();
};

#endif