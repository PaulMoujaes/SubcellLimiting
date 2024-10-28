#include "subcell_clipandscale.hpp"

Subcell_ClipAndScale::Subcell_ClipAndScale(ParFiniteElementSpace &fes_, ParFiniteElementSpace *subcell_fes_,
                              FunctionCoefficient &inflow, VectorCoefficient &velocity,
                              ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction *submesh_vel, int exec_mode_) :
   Subcell_FE_Evolution(fes_, subcell_fes_, inflow, velocity, M, x0_, mesh_vel, submesh_vel, exec_mode_)
{
   umin.SetSize(lumpedmassmatrix.Size());
   umax.SetSize(lumpedmassmatrix.Size());
   udot.SetSize(lumpedmassmatrix.Size());
}

void Subcell_ClipAndScale::ComputeBounds(const Vector &u, Array<double> &u_min,
                                 Array<double> &u_max) const
{
   // iterate over local number of dofs on this processor and compute maximum and minimum over local stencil
   for (int i = 0; i < fes.GetVSize(); i++)
   {
      umin[i] = u(i);
      umax[i] = u(i);

      for (int k = I[i]; k < I[i+1]; k++)
      {
         int j = J[k];
         umin[i] = min(umin[i], u(j));
         umax[i] = max(umax[i], u(j));
      }
   }

   // Distribute min and max to get max and min of local stencil of shared dofs
   gcomm.Reduce<double>(umax, GroupCommunicator::Max);
   gcomm.Bcast(umax);

   gcomm.Reduce<double>(umin, GroupCommunicator::Min);
   gcomm.Bcast(umin);
}

void Subcell_ClipAndScale::Mult(const Vector &x, Vector &y) const
{
   const int dim = fes.GetMesh()->Dimension();
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
   ComputeLOTimeDerivatives(x, udot);
   ComputeBounds(x, umin, umax);

   Array<int> dofs, ddofs;
   for (int e = 0; e < fes.GetNE(); e++)
   {
       
      auto element = fes.GetFE(e);
      auto eltrans = fes.GetElementTransformation(e);
      fes.GetElementDofs(e, dofs);
       
      // assemble element mass and convection matrices
      conv->AssembleElementMatrix(*element, *eltrans, Ke);
      div_int.AssembleElementMatrix2(*element, *element, *eltrans, Ce);
      SparseMatrix Ce_tilde(dofs.Size(), dim * dofs.Size());
      //BuildSubcellDivElementMatrix(e, Ce_tilde);
      BuildSubcellDivElementMatrix2(e, Ce, Ce_tilde);
      mass_int.AssembleElementMatrix(*element, *eltrans, Me);
       
      ue.SetSize(dofs.Size());
      re.SetSize(dofs.Size());
      udote.SetSize(dofs.Size());
      fe.SetSize(dofs.Size());
      fe_star.SetSize(dofs.Size());
      gammae.SetSize(dofs.Size());
      ue_bar.SetSize(dofs.Size());
       
      x.GetSubVector(dofs, ue);
      //MFEM_VERIFY(ue.Min() > -1e-15, "bernstein coefficients negative!")
      udot.GetSubVector(dofs, udote);
       
      v_GFE.ParFESpace()->GetElementVDofs(e, ddofs);
      ve.SetSize(ddofs.Size());
      v_GFE.GetSubVector(ddofs, ve);
      //ve = 1.0;
       
      re = 0.0;
      fe = 0.0;
      gammae = 0.0;
      ue_bar = 0.0;

      auto II = Ce_tilde.GetI();
      auto JJ = Ce_tilde.GetJ();
      auto CC = Ce_tilde.ReadData();

      Vector vi(dim), vj(dim), cij_tilde(dim), cji_tilde(dim), cij(dim), cji(dim);
      for (int i = 0; i < Ce_tilde.Height(); i++)
      {  
         for(int d = 0; d < dim; d++)
         {
            vi(d) = ve(i + d * dofs.Size());
         }

         for(int k = II[i]; k < II[i+1]; k++)
         {
            int j = JJ[k];
            MFEM_VERIFY(abs(CC[k] - Ce_tilde(i,j)) < 1e-15, "index wrong in sparsity pattern" )
            if( j >= i){continue;}

            for(int d = 0; d < dim; d++)
            {
               cij_tilde(d) = Ce_tilde(i, j + d * dofs.Size());
               cji_tilde(d) = Ce_tilde(j, i + d * dofs.Size());
               vj(d) = ve(j + d * dofs.Size());
            }

            double cij_max = max(abs(cij_tilde * vj), abs(cij_tilde * vi));
            double cji_max = max(abs(cji_tilde * vj), abs(cji_tilde * vi));
            double dije_tilde = max(max(cij_tilde * vj, cji_tilde * vi), 0.0);
            dije_tilde = max(cij_max, cji_max);
            double diffusion = dije_tilde * (ue(j) - ue(i));

            re(i) += diffusion;
            re(j) -= diffusion;

            re(i) -= ( (cij_tilde * vj) * ue(j) - (cij_tilde * vi) * ue(i) );
            re(j) -= ( (cji_tilde * vi) * ue(i) - (cji_tilde * vj) * ue(j) );

            // for bounding fluxes
            gammae(i) += dije_tilde;
            gammae(j) += dije_tilde;
             
            // add 2dije * uije and 2djie * ujie
            ue_bar(i) += dije_tilde * (ue(i) + ue(j)) - ( (cij_tilde * vj) * ue(j) - (cij_tilde * vi) * ue(i) );
            ue_bar(j) += dije_tilde * (ue(j) + ue(i)) - ( (cji_tilde * vi) * ue(i) - (cji_tilde * vj) * ue(j) );
             
            // assemble raw antidifussive fluxes f_{i,e} = sum_j m_{ij,e} (udot_i - udot_j) - d_{ij,e} (u_i - u_j)
            // note fije = - fjie
            double fije = - diffusion; // Me(i,j) * (udote(i) - udote(j))
            fe(i) += fije;
            fe(j) -= fije;
         }

         DenseMatrix Cd_tilde;
         Ce_tilde.ToDenseMatrix(Cd_tilde);

         for(int j = 0; j < dofs.Size(); j++)
         {
            for(int d = 0; d < dim; d++)
            {
               cij_tilde(d) = Cd_tilde(i, j + d * dofs.Size());
               cij(d) = Ce(i, j + d * dofs.Size());
               cji(d) = Ce(j, i + d * dofs.Size());
               vj(d) = ve(j + d * dofs.Size());
            }

            fe(i) += Me(i,j) * (udote(i) - udote(j))+ ( (cij_tilde * vj) - (cij * vj) ) * ue(j) - (cji * vj) * ue(j);

         }
      } 

      //*
      //conv->AssembleElementMatrix(*element, *eltrans, Ke);
      Ke.AddMultTranspose(ue, fe);


      gammae *= 2.0;

      double P_plus = 0.0;
      double P_minus = 0.0;

      //Clip
      for (int i = 0; i < dofs.Size(); i++)
      {
         // bounding fluxes to enforce u_i = u_i_min => du/dt >= 0 and vise versa for u_i = u_i_max         
         double fie_max = gammae(i) * (umax[dofs[i]]) - ue_bar(i); //ue_bar(i);
         double fie_min = gammae(i) * (umin[dofs[i]]) - ue_bar(i); //ue_bar(i);
         fie_max = max(0.0, fie_max);
         fie_min = min(0.0, fie_min);
         
         fe_star(i) = min(max(fie_min, fe(i)), fie_max);

         // track positive and negative contributions
         P_plus += max(fe_star(i), 0.0);
         P_minus += min(fe_star(i), 0.0);
      }
      const double P = P_minus + P_plus;

      //and Scale for the sum of fe_star to be 0, i.e., mass conservation
      for (int i = 0; i < dofs.Size(); i++)
      {
         if (fe_star(i) > 0.0 && P > 0.0)
         {
            fe_star(i) *= - P_minus / P_plus;
         }
         else if (fe_star(i) < 0.0 && P < 0.0)
         {
            fe_star(i) *= - P_plus / P_minus;
         }
      }
      // add limited antidiffusive fluxes to element contribution and add to global vector
      re += fe_star;     
      //*/
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
}


Subcell_ClipAndScale::~Subcell_ClipAndScale()
{ }
