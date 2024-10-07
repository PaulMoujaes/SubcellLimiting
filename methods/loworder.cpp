#include "loworder.hpp"

// Implementation of Class LowOrderScheme
LowOrderScheme::LowOrderScheme(ParFiniteElementSpace &fes_,
                               FunctionCoefficient &inflow,
                               VectorCoefficient &velocity, ParBilinearForm &M, const Vector &x0_, ParGridFunction &mesh_vel, int exec_mode_):
    FE_Evolution(fes_, inflow, velocity, M, x0_, mesh_vel, exec_mode_)
{ }

void LowOrderScheme::Mult(const Vector &x, Vector &y) const
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

    ComputeLOTimeDerivatives(x, y);
}

LowOrderScheme::~LowOrderScheme()
{ }
