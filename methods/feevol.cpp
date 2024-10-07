#include "feevol.hpp"

FE_Evolution::FE_Evolution(ParFiniteElementSpace &fes_,
                                 FunctionCoefficient &inflow,
                                 VectorCoefficient &velocity, ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, int exec_mode_) :
   TimeDependentOperator(M.Height()),
   lumpedM(&fes_), fes(fes_), x_now(*fes_.GetMesh()->GetNodes()), 
   gcomm(fes_.GroupComm()), I(M.SpMat().GetI()), J(M.SpMat().GetJ()),
   x0(x0_), v_gf(mesh_vel), v_mesh_coeff(&mesh_vel),
   //b_lumped(&fes), //u_inflow(&fes), 
   conv_int(velocity), conv_int_remap(v_mesh_coeff), mass_int(), remap(exec_mode_ == 1)
{
    lumpedM.AddDomainIntegrator(new LumpedIntegrator(new MassIntegrator()));
    lumpedM.Assemble();
    lumpedM.Finalize();

    lumpedmassmatrix.SetSize(M.Height());
    lumpedM.SpMat().GetDiag(lumpedmassmatrix);
   //u_inflow.ProjectCoefficient(inflow);

   // distribute the lumped mass matrix entries
   Array<double> lumpedmassmatrix_array(lumpedmassmatrix.GetData(),
                                        lumpedmassmatrix.Size());
   gcomm.Reduce<double>(lumpedmassmatrix_array, GroupCommunicator::Sum);
   gcomm.Bcast(lumpedmassmatrix_array);

   // For bound preservation the boundary condition \hat{u} is enforced
   // via a lumped approximation to < (u_h - u_inflow) * min(v * n, 0 ), w >, i.e.,
   // (u_i - (u_inflow)_i) * \int_F \varphi_i * min(v * n, 0).
   // The integral can be implemented as follows:
   //FunctionCoefficient one_coeff(one);
   //b_lumped.AddBdrFaceIntegrator(
   //   new BoundaryFlowIntegrator(one_coeff, velocity, 1.0));
   //b_lumped.Assemble();

   if(remap)
   {
        conv = &conv_int_remap;
   }
   else
   {
        conv = &conv_int;
   }

   z.SetSize(lumpedmassmatrix.Size());
}

void FE_Evolution::ComputeLOTimeDerivatives(const Vector &u,
                                               Vector &udot) const
{
   udot = 0.0;
   const int nE = fes.GetNE();
   Array<int> dofs;

   for (int e = 0; e < nE; e++)
   {
      auto element = fes.GetFE(e);
      auto eltrans = fes.GetElementTransformation(e);

      // assemble element matrix of convection operator
      conv->AssembleElementMatrix(*element, *eltrans, Ke);

      fes.GetElementDofs(e, dofs);
      ue.SetSize(dofs.Size());
      u.GetSubVector(dofs, ue);
      re.SetSize(dofs.Size());
      re = 0.0;

      for (int i = 0; i < dofs.Size(); i++)
      {
         for (int j = 0; j < i; j++)
         {
            // add low-order stabilization with discrete upwinding
            double dije = max(max(Ke(i,j), Ke(j,i)), 0.0);
            double diffusion = dije * (ue(j) - ue(i));

            re(i) += diffusion;
            re(j) -= diffusion;
         }
      }
      // Add -K_e u_e to obtain (-K_e + D_e) u_e and add element contribution
      // to global vector
      Ke.AddMult(ue, re, -1.0);
      udot.AddElementVector(dofs, re);
   }

   // add boundary condition (u - u_inflow) * b. This is under the assumption that b_lumped has been updated
   //subtract(u, u_inflow, z);
   //z *= b_lumped;
   //udot += z;

   // Distribute
   Array<double> udot_array(udot.GetData(), udot.Size());
   gcomm.Reduce<double>(udot_array, GroupCommunicator::Sum);
   gcomm.Bcast(udot_array);

   // apply inverse lumped mass matrix
   udot /= lumpedmassmatrix;
}

FE_Evolution::~FE_Evolution()
{ }