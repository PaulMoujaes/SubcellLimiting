#ifndef SUBCELL_SUBCELLCLIPANDSCALE
#define SUBCELL_SUBCELLCLIPANDSCALE

#include "subcell_feevol.hpp"

class Subcell_ClipAndScale : public Subcell_FE_Evolution
{
private:
   mutable Array<double> umin, umax;
   mutable Vector udot;

   virtual void ComputeBounds(const Vector &u, Array<double> &u_min,
                              Array<double> &u_max) const;
public:
   Subcell_ClipAndScale(ParFiniteElementSpace &fes_, ParFiniteElementSpace *subcell_fes_,
                   FunctionCoefficient &inflow, VectorCoefficient &velocity,
                   ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction *submesh_vel, int exec_mode_);

   virtual void Mult(const Vector &x, Vector &y) const override;

   virtual ~Subcell_ClipAndScale();
};

#endif