#ifndef SUBCELL_LOWORDER
#define SUBCELL_LOWORDER

#include "feevol.hpp"


class LowOrderScheme : public FE_Evolution
{
public:
   LowOrderScheme(ParFiniteElementSpace &fes_,
                  FunctionCoefficient &inflow,
                  VectorCoefficient &velocity, ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, int exec_mode_);

   virtual void Mult(const Vector &x, Vector &y) const override;

   virtual ~LowOrderScheme();
};

#endif