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
                FunctionCoefficient &inflow,
                VectorCoefficient &velocity, ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel);

   virtual void Mult(const Vector &x, Vector &y) const override;

   virtual ~ClipAndScale();
};

#endif