#include "subcell_loworder.hpp"

Subcell_LowOrder::Subcell_LowOrder(ParFiniteElementSpace &fes_, ParFiniteElementSpace &subcell_fes_,
                              FunctionCoefficient &inflow, VectorCoefficient &velocity,
                              ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction &submesh_vel, int exec_mode_) :
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
        add(x0_sub, mt, vsub_gf, xsub_now);

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
      udote.SetSize(dofs.Size());
      fe.SetSize(dofs.Size());
      fe_star.SetSize(dofs.Size());
      gammae.SetSize(dofs.Size());
      ue_bar.SetSize(dofs.Size());
         
      x.GetSubVector(dofs, ue);
         
      re = 0.0;
      fe = 0.0;
      gammae = 0.0;
      ue_bar = 0.0;
      /*
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
             
            // add 2dije * uije and 2djie * ujie
            ue_bar(i) += dije * (ue(i) + ue(j)) - Ke(i,j) * (ue(j) - ue(i));
            ue_bar(j) += dije * (ue(j) + ue(i)) - Ke(j,i) * (ue(i) - ue(j));
             
            // assemble raw antidifussive fluxes f_{i,e} = sum_j m_{ij,e} (udot_i - udot_j) - d_{ij,e} (u_i - u_j)
            // note fije = - fjie
            double fije = Me(i,j) * (udote(i) - udote(j)) - diffusion;
            fe(i) += fije;
            fe(j) -= fije;
         }
      }
      //*/

      auto II = Ke_tilde.GetI();
      auto JJ = Ke_tilde.GetJ();
      auto KK = Ke_tilde.ReadData();

      for(int i = 0; i < Ke_tilde.Height(); i++)
      {
         for(int k = II[i]; k < II[i+1]; k++)
         {
            int j = JJ[k];
            if(j == 0 || j == 2)
            {
               //Ke_tilde(i,j) *= 4.0 / 3.0;
            }
            else if( j>3 && j < 8)
            {
               //Ke_tilde(i,j) *= 2.0 / 3.0;
            }
         }
      }

      Vector rowsums(Ke_tilde.Height());
      Ke_tilde.GetRowSums(rowsums);
      //MFEM_VERIFY(rowsums.Norml2() < 1e-15, "FUCK: " + to_string(rowsums.Norml2()));


      for (int i = 0; i < Ke_tilde.Height(); i++)
      {  
         for(int k = II[i]; k < II[i+1]; k++)
         {
            int j = JJ[k];
            MFEM_VERIFY(abs(KK[k] - Ke_tilde(i,j)) < 1e-15, "index wrong in sparsity pattern" )
            if( j >= i){continue;}
            double dije_tilde = max(max(KK[k], Ke_tilde(j,i)), 0.0);
            double diffusion = dije_tilde * (ue(j) - ue(i));
            re(i) += diffusion;
            re(j) -= diffusion;
         }

      }
      // add convective term
      //re = 0.0;
      Ke_tilde.AddMult(ue, re, -1.0);

      /*

      Vector columnsums(dofs.Size());
      Vector ones = columnsums;
      ones = 1.0;

   

      Ke_tilde.MultTranspose(ones, columnsums);
      Ke.AddMultTranspose(ones, columnsums, -1.0);
      if(columnsums.Norml2() < 1e-14)
      {
         //cout << "There is still hope in this world" << endl;
         //cout << ".";
         //MFEM_ABORT("jo nice")
      }
      else
      {
         //cout << columnsums.Norml2() << " No hope?" << endl;
         if(fes.GetMesh()->Dimension() == 1)
         {
            MFEM_ABORT("Collumnsum not zero: " + to_string(columnsums.Norml2()))
         }
         else if( true)
         {

            cout <<"schade: "<< columnsums.Norml2() <<endl;
            columnsums.Print();
            cout << endl;
            Ke_tilde.MultTranspose(ones, columnsums);

            columnsums.Print();
            cout << endl;
            Ke.MultTranspose(ones, columnsums);
            columnsums.Print();  

            cout <<"-----------" << endl;

         }
         //cout << "/";
      }
      //*/
      //cout << columnsums.Norml2() << endl;
      //MFEM_VERIFY(columnsums.Norml2() < 1e-15, "not mass conserving")
      //Ke_tilde.MultTranspose(ones, columnsums);
      //columnsums.Print();
      //cout << endl;

      //Ke.MultTranspose(ones, columnsums);
      //columnsums.Print();
      

      //cout << "----------------------" << endl;
      /*
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
      //re += fe_star;
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

   //MFEM_ABORT("");
}


Subcell_LowOrder::~Subcell_LowOrder()
{ }
