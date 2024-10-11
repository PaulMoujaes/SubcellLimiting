#ifndef SUBCELL_SUBCELLFEEVOL
#define SUBCELL_SUBCELLFEEVOL

#include "feevol.hpp"

class Subcell_FE_Evolution : public FE_Evolution 
{
protected:
   // Vector &vsub_gf;
   Vector &xsub_now;
   ParGridFunction &vsub_gf;
   const Vector x0_sub;
   ParFiniteElementSpace &subcell_fes;
   //VectorGridFunctionCoefficient v_subcellmesh_coeff;
   Table coarse_to_fine;
   Array< Array<int>* > dofs2subcelldofs;

   mutable DenseMatrix Kse;
   mutable VectorDivergenceIntegrator *div_int;
   //mutable SparseMatrix Ke_tilde;

public:
   Subcell_FE_Evolution(ParFiniteElementSpace &fes_, ParFiniteElementSpace &subcell_fes_,
                   FunctionCoefficient &inflow, VectorCoefficient &velocity,
                   ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction &submesh_vel, int exec_mode_);

   virtual void ComputeHOTimeDerivatives(const Vector &u, Vector &udot) const; 

   virtual void Mult(const Vector &x, Vector &y) const = 0;

   virtual void BuildSubcellElementMatrix(const int e, SparseMatrix &Ke_tilde) const;

   virtual ~Subcell_FE_Evolution();
};

#endif