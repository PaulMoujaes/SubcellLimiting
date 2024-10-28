#include "subcell_feevol.hpp"

Subcell_FE_Evolution::Subcell_FE_Evolution(ParFiniteElementSpace &fes_, ParFiniteElementSpace *subcell_fes_,
                              FunctionCoefficient &inflow, VectorCoefficient &velocity,
                              ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, ParGridFunction *submesh_vel, int exec_mode_) :
                              FE_Evolution(fes_, inflow, velocity, M, x0_, mesh_vel, exec_mode_),
                              subcell_fes(subcell_fes_),  vsub_gf(submesh_vel), xsub_now(NULL), dfes(NULL)//, v_GFE(&fes_) //, x0_sub(*subcell_fes_->GetMesh()->GetNodes())

{
   H1_FECollection fec(fes.GetFE(0)->GetOrder(), fes.GetParMesh()->Dimension(), BasisType::ClosedUniform);
   dfes = new ParFiniteElementSpace(fes.GetParMesh(), &fec, fes.GetParMesh()->Dimension());//fes.FEColl()
   v_GFE.SetSpace(dfes);
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

   //v_GFE.Print();

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

void Subcell_FE_Evolution::AdjustSubcellElementMatrix(const DenseMatrix &Ce, SparseMatrix &Ce_tilde) const
{
   //DenseMatrix C;
   //Ce_tilde.ToDenseMatrix(K);
   //DenseMatrix Kt = K;
   //Kt = 0.0;
   const int dim = fes.GetMesh()->Dimension();
   const int ndofs = Ce.Height();

   SparseMatrix ce_start = Ce_tilde;
   ce_start *= 1.0;

   auto I = Ce_tilde.GetI();
   auto J = Ce_tilde.GetJ();
   auto CC = Ce_tilde.ReadWriteData();

   Vector ones(ndofs);
   ones = 1.0;
   Vector collumnsums(ndofs), alpha(ndofs);

   GMRESSolver bcg;
   bcg.SetAbsTol(1e-30);
   //bcg.SetPrintLevel(1);
   bcg.SetMaxIter(100);

   Vector aux(Ce.Width());
   Ce.MultTranspose(ones, aux);
    

   for(int d = 0; d < dim; d++)
   {  
      collumnsums.SetData(aux.GetData() + d * ndofs);
      SparseMatrix C_kt(ndofs, ndofs);
       
      for(int i = 0; i < ndofs; i++)
      {
          
         for(int k = I[i]; k < I[i+1]; k++)
         {
            int j = J[k];
            if(j >= ndofs) {continue;}

            double cij_k = Ce_tilde(i, j + d * ndofs); 
            C_kt.Set(j,i,cij_k);
            //collumnsums(j) += cij_k;
             
         }
      }
       
      C_kt.Finalize(0);
       
      bcg.SetOperator(C_kt);
      bcg.Mult(collumnsums, alpha);

      //alpha.Print();
      //cout << endl;
       
      for(int i = 0; i < ndofs; i++)
      {
          
         for(int k = I[i]; k < I[i+1]; k++)
         {
            int j = J[k];
            if(j >= ndofs) {continue;}
            Ce_tilde(i, j + d * ndofs) *= alpha(i);
         }
      }
   }
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

      div_int.AssembleElementMatrix2(*element, *element, *eltrans, Cse);

      Array<int> subcell_dofs = *dofs2subcelldofs[el];
      Array<int> subcell_ddofs(dim * subcell_dofs.Size());

      for(int d = 0; d < dim; d++)
      {
         for(int i = 0; i < subcell_dofs.Size(); i++)
         {
            subcell_ddofs[i + d * subcell_dofs.Size()] = subcell_dofs[i] + d * dofs.Size();
         }
      }

      Ce_tilde.AddSubMatrix(subcell_dofs, subcell_ddofs, Cse, 0);

   }
   Ce_tilde.Finalize(0);
}

void Subcell_FE_Evolution::BuildSubcellDivElementMatrix2(const int e, const DenseMatrix &Ce, SparseMatrix &Ce_tilde) const
{
   const int dim = fes.GetMesh()->Dimension();
   const int ndofs = Ce.Height();
   SparseMatrix Me_tilde(ndofs, ndofs);
   BuildSubcellMasstMatrix(e, Me_tilde);

   Ce_tilde.OverrideSize(ndofs, dim * ndofs);

   Vector ml(ndofs);
   Me_tilde.GetRowSums(ml);

   SparseMatrix ML(ml);
   ML.Finalize();
   ML *= -1.0;
   SparseMatrix MC_ML = Me_tilde;
   MC_ML += ML;

   Vector collumnsums(Ce.Width());
   Vector ones(ndofs);
   ones = 1.0;

   Ce.MultTranspose(ones, collumnsums);

   Vector cs_d(ndofs);

   auto I = Me_tilde.GetI();
   auto J = Me_tilde.GetJ();
   
   GMRESSolver gmres;
   gmres.SetAbsTol(1e-30);
   gmres.SetMaxIter(100);
   gmres.SetOperator(MC_ML);

   Vector v_d(ndofs);

   for(int d = 0; d < dim; d++)
   {
      cs_d.SetData(collumnsums.GetData() + d * ndofs);

      gmres.Mult(cs_d, v_d);

      double coeff = - v_d.Sum() / double(ndofs);

      v_d.Add(coeff, ones);

      //if(abs(v_d.Sum()) > 1e-14)
      //{
      //   cout << "v sum = " << v_d.Sum() << endl;
      //}

      for(int i = 0; i < ndofs; i++)
      {
         for(int k = I[i]; k < I[i+1]; k++)
         {
            int j = J[k];
            double cij_d = Me_tilde(i,j) * v_d(i);
            if(i == j) 
            {
               cij_d -= ml(i) * v_d(i); 
            }

            Ce_tilde.Set(i, j + d * ndofs, cij_d);
         }
      }
   }
   Ce_tilde.Finalize(0);

   /*
   Ce_tilde.AddMultTranspose(ones, collumnsums, -1.0);



   if(collumnsums.Norml2() > 1e-14)
   {
      cout << "collumnsums = " << collumnsums.Norml2() << endl;
   }
   //*/
}

void Subcell_FE_Evolution::BuildSubcellMasstMatrix(const int e, SparseMatrix &Me_tilde) const
{
   Me_tilde = 0.0;
   Array<int> dofs, subcell_el;
   fes.GetElementDofs(e, dofs);
   Me_tilde.OverrideSize(dofs.Size(), dofs.Size());
   coarse_to_fine.GetRow(e, subcell_el);
    
   for(int el = 0; el < subcell_el.Size(); el++)
   {
      int es = subcell_el[el];
      auto element = subcell_fes->GetFE(es);
      auto eltrans = subcell_fes->GetElementTransformation(es);

       
      // assemble mass matrix on subcell
      mass_int.AssembleElementMatrix(*element, *eltrans, Kse);
      Array<int> subcell_dofs = *dofs2subcelldofs[el];
       
      Me_tilde.AddSubMatrix(subcell_dofs, subcell_dofs, Kse, 0);
       
   }
   Me_tilde.Finalize(0);
    
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