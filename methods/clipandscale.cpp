#include "clipandscale.hpp"

ClipAndScale::ClipAndScale(ParFiniteElementSpace &fes_, FunctionCoefficient &inflow,
                           VectorCoefficient &velocity, ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, int exec_mode_):
   FE_Evolution(fes_, inflow, velocity, M, x0_, mesh_vel, exec_mode_)
{
   umin.SetSize(lumpedmassmatrix.Size());
   umax.SetSize(lumpedmassmatrix.Size());
   udot.SetSize(lumpedmassmatrix.Size());
}

void ClipAndScale::ComputeBounds(const Vector &u, Array<double> &u_min,
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

void ClipAndScale::Mult(const Vector &x, Vector &y) const
{
     
    if(remap)
    { 
        const double t = GetTime();
        // since v_gf has is multiplied with -1 for the convection integrator to have the correct direction
        double mt = - t;
        add(x0, mt, v_gf, x_now);

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
    //udot = 0.0;
    ComputeBounds(x, umin, umax);

   Array<int> dofs;
   for (int e = 0; e < fes.GetNE(); e++)
   {
      auto element = fes.GetFE(e);
      auto eltrans = fes.GetElementTransformation(e);

      // assemble element mass and convection matrices
      conv->AssembleElementMatrix(*element, *eltrans, Ke);
      mass_int.AssembleElementMatrix(*element, *eltrans, Me);

      fes.GetElementDofs(e, dofs);
      ue.SetSize(dofs.Size());
      re.SetSize(dofs.Size());
      udote.SetSize(dofs.Size());
      fe.SetSize(dofs.Size());
      fe_star.SetSize(dofs.Size());
      gammae.SetSize(dofs.Size());
         
      x.GetSubVector(dofs, ue);
      udot.GetSubVector(dofs, udote);
         
      re = 0.0;
      fe = 0.0;
      gammae = 0.0;
      for (int i = 0; i < dofs.Size(); i++)
      {
         for (int j = 0; j < i; j++)
         {
            // add low-order diffusion
            // note that dije = djie
            double dije = max(max(Ke(i,j), Ke(j,i)), 0.0);
            double diffusion = dije * (ue(j) - ue(i));

            re(i) += diffusion;
            re(j) -= diffusion;

            // for bounding fluxes
            gammae(i) += dije;
            gammae(j) += dije;
             
            // assemble raw antidifussive fluxes f_{i,e} = sum_j m_{ij,e} (udot_i - udot_j) - d_{ij,e} (u_i - u_j)
            // note fije = - fjie
            double fije = Me(i,j) * (udote(i) - udote(j)) - diffusion;
            fe(i) += fije;
            fe(j) -= fije;
         }
      }

      // add convective term
      Ke.AddMult(ue, re, -1.0);

      gammae *= 2.0;

      double P_plus = 0.0;
      double P_minus = 0.0;

      //Clip
      for (int i = 0; i < dofs.Size(); i++)
      {
         
         // bounding fluxes to enforce u_i = u_i_min => du/dt >= 0 and vise versa for u_i = u_i_max         
         double fie_max = gammae(i) * (umax[dofs[i]] - ue(i));
         double fie_min = gammae(i) * (umin[dofs[i]] - ue(i));
         
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


ClipAndScale::~ClipAndScale()
{ }
