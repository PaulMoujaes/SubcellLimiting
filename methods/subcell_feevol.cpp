#include "subcell_feevol.hpp"

Subcell_FE_Evolution::Subcell_FE_Evolution(ParFiniteElementSpace &fes_, ParFiniteElementSpace *subcell_fes_,
                              FunctionCoefficient &inflow, VectorCoefficient &velocity,
                              ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction *submesh_vel, int exec_mode_) :
                              FE_Evolution(fes_, inflow, velocity, M, x0_, mesh_vel, exec_mode_),
                              subcell_fes(subcell_fes_),  vsub_gf(submesh_vel), xsub_now(NULL) //, x0_sub(*subcell_fes_->GetMesh()->GetNodes())

{
   if(remap)
   {
      xsub_now = subcell_fes_->GetMesh()->GetNodes();
      x0_sub = *subcell_fes_->GetMesh()->GetNodes();
   }  
   /*
   FiniteElementSpace dfes(fes.GetMesh(), fes.FEColl(), fes.GetMesh()->Dimension());
   FiniteElementSpace dfes_sub(subcell_fes->GetMesh(), subcell_fes->FEColl(), subcell_fes->GetMesh()->Dimension());
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

void Subcell_FE_Evolution::AdjustSubcellElementMatrix(const DenseMatrix &Ke, SparseMatrix &Ke_tilde) const
{
   DenseMatrix K;
   Ke_tilde.ToDenseMatrix(K);
   DenseMatrix Kt = K;
   Kt = 0.0;

   for(int i = 0; i < K.Height(); i++)
   {
      for(int j = 0; j < K.Height(); j++)
      {
         Kt(i,j) = K(j,i);
      }
   }

   Vector ones(Ke.Height());
   ones = 1.0;
   Vector collumnsums(ones.Size()), a(ones.Size()), alpha(ones.Size());

   GMRESSolver bcg;
   bcg.SetOperator(Kt);
   Ke.MultTranspose(ones, collumnsums);
   bcg.SetAbsTol(1e-30);
   //bcg.SetPrintLevel(1);
   bcg.SetMaxIter(100);
   bcg.Mult(collumnsums, a);

   //add(a, 1.0 - a(0), ones, alpha);
   alpha = a;

   Ke_tilde.AddMultTranspose(a, collumnsums, -1.0);
 

   //alpha.Print();
   //cout << " >>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<" <<  endl;
   

   /*
   alpha.Print();
   MFEM_ABORT("hmm")

   for(int i = 0; i < Kt.Width(); i++)
   {
      if(abs(collumnsums(i)) < 1e-15)
      {
         for(int j = 0; j < Kt.Height(); j++)
         {
            Kt(j,i) = 0.0;
         }
         //collumnsums(i) = 1.0;
         Kt(i,i) = 1.0;
         break;
      }
   }
   Kt.Print();
   Kt.Invert();
   Kt.Print();
   MFEM_ABORT("")
   Kt.MultTranspose(collumnsums, alpha);
   MFEM_VERIFY(alpha.Norml2() > 1e-15, "alpha = 0");
   
   Vector aux = alpha;
   Ke_tilde.MultTranspose(alpha, aux);
   aux -= collumnsums;
   cout << "|ktilde alpha| = " <<aux.Norml2() << endl;
   aux.Print();
   cout << "---" << endl;

   Ke.MultTranspose(ones, collumnsums);
   Vector collumnsums2 = collumnsums;
   Ke_tilde.MultTranspose(ones, collumnsums2);
   collumnsums2 -= collumnsums;
   //double init = collumnsums.Norml2();

   //Ke_tilde.AddMultTranspose(ones, collumnsums, 1.0);

   //*/

   auto I = Ke_tilde.GetI();
   auto J = Ke_tilde.GetJ();
   auto KK = Ke_tilde.ReadWriteData();

   for(int i = 0; i < Ke_tilde.Height(); i++)
   {
      for(int k = I[i]; k < I[i+1]; k++)
      {
         KK[k] *= alpha(i);
      }
   }

   //Vector collumnsums3 = collumnsums;
   //Ke_tilde.MultTranspose(ones, collumnsums3);
   //collumnsums3 -= collumnsums;
   //collumnsums2.Print();
   //collumnsums3.Print();

   //MFEM_VERIFY(collumnsums.Norml2() < 1e-13, to_string(collumnsums.Norml2()) );
   //cout << endl;
   //cout << "LETSGOO" << endl;

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
      Array<int> subcell_dofs2;
      subcell_fes->GetElementDofs(es, subcell_dofs2);

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