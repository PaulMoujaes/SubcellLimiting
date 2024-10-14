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
      udote.SetSize(dofs.Size());
      fe.SetSize(dofs.Size());
      fe_star.SetSize(dofs.Size());
      gammae.SetSize(dofs.Size());
      ue_bar.SetSize(dofs.Size());

      x.GetSubVector(dofs, ue);
      udot.GetSubVector(dofs, udote);
         
      re = 0.0;
      fe = 0.0;
      gammae = 0.0;
      ue_bar = 0.0;

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
            double kije_tilde = KK[k];
            double kjie_tilde = Ke_tilde(j,i);
            double dije_tilde = max(max(kije_tilde, kjie_tilde), 0.0);
            double diffusion = dije_tilde * (ue(j) - ue(i));

            re(i) += diffusion;
            re(j) -= diffusion;


            // for bounding fluxes
            gammae(i) += dije_tilde;
            gammae(j) += dije_tilde;
             
            // add 2dije * uije and 2djie * ujie
            ue_bar(i) += dije_tilde * (ue(i) + ue(j)) - kije_tilde * (ue(j) - ue(i));
            ue_bar(j) += dije_tilde * (ue(j) + ue(i)) - kjie_tilde * (ue(i) - ue(j));
             
            // assemble raw antidifussive fluxes f_{i,e} = sum_j m_{ij,e} (udot_i - udot_j) - d_{ij,e} (u_i - u_j)
            // note fije = - fjie
            double fije = Me(i,j) * (udote(i) - udote(j)) - diffusion;
            fe(i) += fije;
            fe(j) -= fije;
         }
      }
      // add terms for sparsity pattern correction
      Ke_tilde.AddMult(ue, fe, 1.0);
      Ke.AddMult(ue, fe, -1.0);

      // add convective term
      Ke_tilde.AddMult(ue, re, -1.0);

      gammae *= 2.0;

      double P_plus = 0.0;
      double P_minus = 0.0;

      //Clip
      for (int i = 0; i < dofs.Size(); i++)
      {
         
         // bounding fluxes to enforce u_i = u_i_min => du/dt >= 0 and vise versa for u_i = u_i_max         
         double fie_max = gammae(i) * umax[dofs[i]] - ue_bar(i);
         double fie_min = gammae(i) * umin[dofs[i]] - ue_bar(i);
         
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

      //MFEM_VERIFY(abs(fe_star.Sum())< 1e-14, "hmm" );
      /*
      if(abs(fe.Sum()) > 1e-15)
      {
         cout << "Target Scheme nicht masseerhaltend: " << fe.Sum() << endl;
         if(abs(fe_star.Sum()) > 1e-15)
         {
            cout << "C&S nicht masseerhaltend: " << fe.Sum() << endl;
            MFEM_ABORT("C&S nicht masseerhaltend")
         }
         else
         {
            cout << "C&S masseerhaltend " << endl;
         }
      }
      if(fe.Norml2() > 1e-15)
      {
         cout << "vorher " << fe.Norml2() << endl;
         cout << "nachher " << fe_star.Norml2() << "\n\n";
      }
      //*/
     
      y.AddElementVector(dofs, re);
   }
   //MFEM_ABORT("alles ok")

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


Subcell_ClipAndScale::~Subcell_ClipAndScale()
{ }
