#ifndef SUBCELL_FEEVOL
#define SUBCELL_FEEVOL

#include "mfem.hpp"
using namespace std;
using namespace mfem;

class FE_Evolution : public TimeDependentOperator 
{
protected:
   const Vector &lumpedmassmatrix;
   //DofInfo &dofs;
   ParFiniteElementSpace &fes;
   GroupCommunicator &gcomm;
   int *I, *J;
   ParLinearForm b_lumped;
   //ParGridFunction u_inflow;

   mutable DenseMatrix Ke, Me;
   mutable Vector ue, re, udote, fe, fe_star, gammae;
   mutable ConvectionIntegrator conv_int;
   mutable MassIntegrator mass_int;
   mutable Vector z;

   virtual void ComputeLOTimeDerivatives(const Vector &u, Vector &udot) const;

public:
   FE_Evolution(ParFiniteElementSpace &fes_,
                   const Vector &lumpedmassmatrix_, FunctionCoefficient &inflow,
                   VectorFunctionCoefficient &velocity,
                   ParBilinearForm &M);

   virtual void Mult(const Vector &x, Vector &y) const = 0;

   virtual ~FE_Evolution();
};

#endif