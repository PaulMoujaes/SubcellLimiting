#include "subcell_feevol.hpp"

Subcell_FE_Evolution::Subcell_FE_Evolution(ParFiniteElementSpace &fes_, ParFiniteElementSpace &subcell_fes_,
                              FunctionCoefficient &inflow, VectorCoefficient &velocity,
                              ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction &submesh_vel, int exec_mode_) :
                              FE_Evolution(fes_, inflow, velocity, M, x0_, mesh_vel, exec_mode_),
                              xsub_now(*subcell_fes_.GetMesh()->GetNodes()), subcell_fes(subcell_fes_), x0_sub(*subcell_fes_.GetMesh()->GetNodes()),  vsub_gf(submesh_vel)

{
   /*
   FiniteElementSpace dfes(fes.GetMesh(), fes.FEColl(), fes.GetMesh()->Dimension());
   FiniteElementSpace dfes_sub(subcell_fes.GetMesh(), subcell_fes.FEColl(), subcell_fes.GetMesh()->Dimension());
   MixedBilinearForm k_s(&dfes_sub, &subcell_fes);
   MixedBilinearForm k(&dfes, &fes);

   k_s.AddDomainIntegrator(new VectorDivergenceIntegrator ());
   k.AddDomainIntegrator(new VectorDivergenceIntegrator ());

   div_int = new VectorDivergenceIntegrator();

   k_s.Assemble();
   k_s.Finalize(0);

   k.Assemble();
   k.Finalize(0);

   Vector sums(k.Width());
   Vector ones(k.Height());
   ones = 1.0;

   k.MultTranspose(ones, sums);
   k_s.AddMultTranspose(ones, sums, -1.0);

   MFEM_VERIFY(k.Height() == k_s.Height(), "not same hight");
   MFEM_VERIFY(k.Width() == k_s.Width(), "not same Width");

   cout << k.Height() << " x " << k.Width() << endl;
   

   cout << sums.Norml2() << endl;
   //*/

   //MFEM_ABORT("")

   subcell_fes.GetMesh()->GetRefinementTransforms().MakeCoarseToFineTable(coarse_to_fine, true);

   Array<int> dofs, subcell_dofs, subcell_el;
   fes.GetElementDofs(0, dofs);
   coarse_to_fine.GetRow(0, subcell_el);
   subcell_fes.GetElementDofs(subcell_el[0], subcell_dofs);
   dofs2subcelldofs.SetSize(subcell_el.Size());

   for(int el = 0; el < subcell_el.Size(); el++)
   {
      int e = subcell_el[el];
      subcell_fes.GetElementDofs(e, subcell_dofs);
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
      auto element = subcell_fes.GetFE(es);
      auto eltrans = subcell_fes.GetElementTransformation(es);

      
      // assemble convection matrix on subcell
      conv->AssembleElementMatrix(*element, *eltrans, Kse);
      Array<int> subcell_dofs = *dofs2subcelldofs[el];

      Ke_tilde.AddSubMatrix(subcell_dofs, subcell_dofs, Kse, 0);
      Array<int> subcell_dofs2;
      subcell_fes.GetElementDofs(es, subcell_dofs2);

      for(int i = 0; i < subcell_dofs.Size(); i++)
      {
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         MFEM_VERIFY(subcell_dofs2[i] == dofs[subcell_dofs[i]], "subcell dofs to dofs wrong in element " + to_string(e) + " on rank " + to_string(rank));
         //cout << "jupp" << endl;
      }
   }
   Ke_tilde.Finalize(0);
}

Subcell_FE_Evolution::~Subcell_FE_Evolution()
{ }