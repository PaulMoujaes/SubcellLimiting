#ifndef SUBCELL_FEEVOL
#define SUBCELL_FEEVOL

#include "mfem.hpp"
using namespace std;
using namespace mfem;

class FE_Evolution : public TimeDependentOperator 
{
protected:
   mutable Vector lumpedmassmatrix;
   Vector &v_gf;
   mutable ParBilinearForm lumpedM;
   Vector &x_now;
   const Vector &x0;
   ParFiniteElementSpace &fes;
   GroupCommunicator &gcomm;
   int *I, *J;
   VectorGridFunctionCoefficient v_mesh_coeff;
   bool remap;
   //ParLinearForm b_lumped;
   //ParGridFunction u_inflow;

   mutable DenseMatrix Ke, Me;
   mutable Vector ue, re, udote, fe, fe_star, gammae;
   mutable ConvectionIntegrator conv_int, conv_int_remap;
   mutable ConvectionIntegrator *conv;
   mutable MassIntegrator mass_int;
   mutable Vector z;

   virtual void ComputeLOTimeDerivatives(const Vector &u, Vector &udot) const;

public:
   FE_Evolution(ParFiniteElementSpace &fes_,
                   FunctionCoefficient &inflow,
                   VectorCoefficient &velocity,
                   ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, int exec_mode_);

   virtual void Mult(const Vector &x, Vector &y) const = 0;

   virtual ~FE_Evolution();
};

#endif