#include "subcell_loworder.hpp"

Subcell_LowOrder::Subcell_LowOrder(ParFiniteElementSpace &fes_, ParFiniteElementSpace *subcell_fes_,
                              FunctionCoefficient &inflow, VectorCoefficient &velocity,
                              ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction *submesh_vel, int exec_mode_) :
   Subcell_FE_Evolution(fes_, subcell_fes_, inflow, velocity, M, x0_, mesh_vel, submesh_vel, exec_mode_)
{
}

void Subcell_LowOrder::Mult(const Vector &x, Vector &y) const
{
     
    if(remap)
    { 
        const double t = GetTime();
        // since v_gf has is multiplied with -1 for the convection integrator to have the correct direction
        double mt = - t;
        add(x0, mt, v_gf, x_now);
        add(x0_sub, mt, *vsub_gf, *xsub_now);

        lumpedM.BilinearForm::operator=(0.0);
        lumpedM.Assemble();
        lumpedM.SpMat().GetDiag(lumpedmassmatrix);
        Array<double> lumpedmassmatrix_array(lumpedmassmatrix.GetData(),
                                        lumpedmassmatrix.Size());
        gcomm.Reduce<double>(lumpedmassmatrix_array, GroupCommunicator::Sum);
        gcomm.Bcast(lumpedmassmatrix_array);
    }

    y = 0.0;

    // compute low-order time derivative for high-order stabilization and local bounds
    //ComputeLOTimeDerivatives(x, udot);
    //udot = 0.0;
   //ComputeBounds(x, umin, umax);

   Array<int> dofs;
   for (int e = 0; e < fes.GetNE(); e++)
   {
      auto element = fes.GetFE(e);
      auto eltrans = fes.GetElementTransformation(e);
      fes.GetElementDofs(e, dofs);

      // assemble element mass and convection matrices
      conv->AssembleElementMatrix(*element, *eltrans, Ke);
      SparseMatrix Ke_tilde(dofs.Size());
      BuildSubcellElementMatrix(e, Ke_tilde);
      mass_int.AssembleElementMatrix(*element, *eltrans, Me);

      ue.SetSize(dofs.Size());
      re.SetSize(dofs.Size());

      x.GetSubVector(dofs, ue);
         
      re = 0.0;

      auto II = Ke_tilde.GetI();
      auto JJ = Ke_tilde.GetJ();
      auto KK = Ke_tilde.ReadData();

      for (int i = 0; i < Ke_tilde.Height(); i++)
      {  
         for(int k = II[i]; k < II[i+1]; k++)
         {
            int j = JJ[k];
            //MFEM_VERIFY(abs(KK[k] - Ke_tilde(i,j)) < 1e-15, "index wrong in sparsity pattern" )
            if( j >= i){continue;}
            double dije_tilde = max(max(KK[k], Ke_tilde(j,i)), 0.0);
            double diffusion = dije_tilde * (ue(j) - ue(i));
            re(i) += diffusion;
            re(j) -= diffusion;
         }

      }
      // add convective term
      Ke_tilde.AddMult(ue, re, -1.0);
     
      y.AddElementVector(dofs, re);
   }

   // add boundary condition (u - u_inflow) * b
   //subtract(x, u_inflow, z);
   //z *= b_lumped;
   //y += z;

   // distribute
   Array<double> y_array(y.GetData(), y.Size());
   gcomm.Reduce<double>(y_array, GroupCommunicator::Sum);
   gcomm.Bcast(y_array);

   // apply inverse lumped mass matrix
   y /= lumpedmassmatrix;

   //MFEM_ABORT("");
}


Subcell_LowOrder::~Subcell_LowOrder()
{ }
