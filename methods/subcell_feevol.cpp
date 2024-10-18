#include "subcell_feevol.hpp"

Subcell_FE_Evolution::Subcell_FE_Evolution(ParFiniteElementSpace &fes_, ParFiniteElementSpace *subcell_fes_,
                              FunctionCoefficient &inflow, VectorCoefficient &velocity,
                              ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction *submesh_vel, int exec_mode_) :
                              FE_Evolution(fes_, inflow, velocity, M, x0_, mesh_vel, exec_mode_),
                              subcell_fes(subcell_fes_),  vsub_gf(submesh_vel), xsub_now(NULL), v_GFE(&fes_) //, x0_sub(*subcell_fes_->GetMesh()->GetNodes())

{
   if(remap)
   {
      xsub_now = subcell_fes_->GetMesh()->GetNodes();
      x0_sub = *subcell_fes_->GetMesh()->GetNodes();
      v_GFE.ProjectCoefficient(v_mesh_coeff);
   } 
   else
   {
      v_GFE.ProjectCoefficient(velocity);
   }

   subcell_fes->GetMesh()->GetRefinementTransforms().MakeCoarseToFineTable(coarse_to_fine, true);

   Array<int> dofs, subcell_dofs, subcell_el;
   fes.GetElementDofs(0, dofs);
   coarse_to_fine.GetRow(0, subcell_el);
   subcell_fes->GetElementDofs(subcell_el[0], subcell_dofs);
   dofs2subcelldofs.SetSize(subcell_el.Size());

   for(int el = 0; el < subcell_el.Size(); el++)
   {
      int e = subcell_el[el];
      subcell_fes->GetElementDofs(e, subcell_dofs);
      dofs2subcelldofs[el] = new Array<int>(subcell_dofs.Size());

      for(int i = 0; i < subcell_dofs.Size(); i++)
      {  
         dofs2subcelldofs[el]->operator[](i) = -1;
         for(int j = 0; j < dofs.Size(); j++)
         {  
            if(subcell_dofs[i] == dofs[j])
            {
               dofs2subcelldofs[el]->operator[](i) = j;
               break;
            }
         }
         MFEM_VERIFY(dofs2subcelldofs[el]->operator[](i) != -1, "Subcell dof not found in element!");
         //*/
      }
   }
   //MFEM_ABORT("hier alles ok")
   for(int e = 0; e < dofs2subcelldofs.Size(); e++)
   {
      //dofs2subcelldofs[e]->Print();
   }
}

void Subcell_FE_Evolution::ComputeHOTimeDerivatives(const Vector &u,
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

      // Set r = -K_e u_e add element contribution
      // to global vector
      re = 0.0;
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

void Subcell_FE_Evolution::BuildSubcellDivElementMatrix(const int e, SparseMatrix &Ce_tilde) const
{
   int dim = fes.GetMesh()->Dimension();
   Ce_tilde = 0.0;
   Array<int> dofs, subcell_el;
   fes.GetElementDofs(e, dofs);
   Ce_tilde.OverrideSize(dofs.Size(), dim * dofs.Size());
   coarse_to_fine.GetRow(e, subcell_el);

   for(int el = 0; el < subcell_el.Size(); el++)
   {
      int es = subcell_el[el];
      auto element = subcell_fes->GetFE(es);
      auto eltrans = subcell_fes->GetElementTransformation(es);

      sdiv_int.AssembleElementMatrix(*element, *eltrans, Cse);

      Array<int> subcell_dofs = *dofs2subcelldofs[el];
      Array<int> subcell_vdofs(dim * subcell_dofs.Size());

      for(int d = 0; d < dim; d++)
      {
         for(int i = 0; i < subcell_dofs.Size(); i++)
         {
            subcell_vdofs[i + d * subcell_dofs.Size()] = subcell_dofs[i] + d * subcell_dofs.Size();
         }
      }

      Ce_tilde.AddSubMatrix(subcell_dofs, subcell_vdofs, Cse, 0);

   }
   Ce_tilde.Finalize(0);
}

void Subcell_FE_Evolution::BuildSubcellElementMatrix(const int e, SparseMatrix &Ke_tilde) const
{  
   Ke_tilde = 0.0;
   Array<int> dofs, subcell_el;
   fes.GetElementDofs(e, dofs);
   Ke_tilde.OverrideSize(dofs.Size(), dofs.Size());
   coarse_to_fine.GetRow(e, subcell_el);

   for(int el = 0; el < subcell_el.Size(); el++)
   {
      int es = subcell_el[el];
      auto element = subcell_fes->GetFE(es);
      auto eltrans = subcell_fes->GetElementTransformation(es);

      
      // assemble convection matrix on subcell
      conv->AssembleElementMatrix(*element, *eltrans, Kse);
      Array<int> subcell_dofs = *dofs2subcelldofs[el];

      Ke_tilde.AddSubMatrix(subcell_dofs, subcell_dofs, Kse, 0);

      /*
      Array<int> subcell_dofs2;
      subcell_fes->GetElementDofs(es, subcell_dofs2);

      for(int i = 0; i < subcell_dofs.Size(); i++)
      {
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         MFEM_VERIFY(subcell_dofs2[i] == dofs[subcell_dofs[i]], "subcell dofs to dofs wrong in element " + to_string(e) + " on rank " + to_string(rank));
         //cout << "jupp" << endl;
      }
      //*/
   }
   Ke_tilde.Finalize(0);
}

Subcell_FE_Evolution::~Subcell_FE_Evolution()
{ }