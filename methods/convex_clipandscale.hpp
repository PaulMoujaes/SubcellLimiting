#ifndef SUBCELL_CONVEXCLIPANDSCALE
#define SUBCELL_CONVEXCLIPANDSCALE

#include "feevol.hpp"

class Convex_ClipAndScale : public FE_Evolution
{
private:
   mutable Array<double> umin, umax;
   mutable Vector udot;

   virtual void ComputeBounds(const Vector &u, Array<double> &u_min,
                              Array<double> &u_max) const;

public:
   Convex_ClipAndScale(ParFiniteElementSpace &fes_,
                FunctionCoefficient &inflow,
                VectorCoefficient &velocity, ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, int exec_mode_);

   virtual void Mult(const Vector &x, Vector &y) const override;

   virtual ~Convex_ClipAndScale();
};

#endif