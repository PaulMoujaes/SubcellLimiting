#ifndef SUBCELL_SUBCELLLOWORDER
#define SUBCELL_SUBCELLLOWORDER

#include "subcell_feevol.hpp"

class Subcell_LowOrder : public Subcell_FE_Evolution
{
private:

public:
   Subcell_LowOrder(ParFiniteElementSpace &fes_, ParFiniteElementSpace &subcell_fes_,
                   FunctionCoefficient &inflow, VectorCoefficient &velocity,
                   ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction &submesh_vel, int exec_mode_);

   virtual void Mult(const Vector &x, Vector &y) const override;

   virtual ~Subcell_LowOrder();
};

#endif