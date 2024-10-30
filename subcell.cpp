
#include <fstream>
#include <iostream>
#include "methods/methodslib.hpp"

using namespace std;
using namespace mfem;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);
void test_velocity_function(const Vector &x, Vector &v)
{
    v = x;
}

// Initial condition
double u0_function(const Vector &x);

// Inflow boundary condition
double inflow_function(const Vector &x);

// Function f = 1 for lumped boundary operator
double one(const Vector &x) {return 1.0;}

int problem, exec_mode;

// Mesh bounding box
Vector bb_min, bb_max;

int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE.
    Mpi::Init();
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // 2. Parse command-line options.
    problem = 0;
    bool subcell = false;
    int mesh_order = 3;
    const char *mesh_file = "../data/periodic-hexagon.mesh";
    int ser_ref_levels = 2;
    int par_ref_levels = 0;
    int order = 3;
    int ode_solver_type = 4;
    int scheme = 1;
    double t_final = 1.0;
    double dt = 0.01;
    bool visualization = true;
    int vis_steps = 50;
    int precision = 8;

    cout.precision(precision);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                    "Mesh file to use.");
    args.AddOption(&problem, "-p", "--problem",
                    "Problem setup to use. See options in velocity_function().");
    //args.AddOption(&exec_mode, "-e", "--execute-mode",
    //                "0 for standard linear advection, 1 for remap mode with moving mesh.");
    args.AddOption(&mesh_order, "-mo", "--mesh-order",
                    "order of the mesh.");
    args.AddOption(&ser_ref_levels, "-r", "--refine-serial",
                    "Number of times to refine the mesh uniformly in serial.");
    args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                    "Number of times to refine the mesh uniformly in parallel.");
    args.AddOption(&order, "-o", "--order",
                    "Order (degree) of the finite elements.");
    args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                    "ODE solver: 1 - Forward Euler,\n\t"
                    "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6,\n\t"
                    "            11 - Backward Euler,\n\t"
                    "            12 - SDIRK23 (L-stable), 13 - SDIRK33,\n\t"
                    "            22 - Implicit Midpoint Method,\n\t"
                    "            23 - SDIRK23 (A-stable), 24 - SDIRK34");
    args.AddOption(&scheme, "-sc", "--scheme",
                    "Finite Element scheme: 1 - Standard Discontinuous Galerkin method,\n\t"
                    "                       11 - Clip and Scale Limiter for continuous Galerkin discretization,\n\t"
                    "                       12 - High-order target schme for continuous Galerkin discretization,\n\t"
                    "                       13 - Low-order schme for continuous Galerkin discretization");
    args.AddOption(&t_final, "-tf", "--t-final",
                    "Final time; start time is 0.");
    args.AddOption(&dt, "-dt", "--time-step",
                    "Time step.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                    "--no-visualization",
                    "Enable or disable GLVis visualization.");
    args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                    "Visualize every n-th timestep.");
    args.Parse();
    if (!args.Good())
    {
        if (Mpi::Root())
        {
            args.PrintUsage(cout);
        }
        return 1;
    }
    if (Mpi::Root())
    {
        args.PrintOptions(cout);
    }

    exec_mode = (int) (problem >= 10);
    subcell = (scheme > 2);

    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    ODESolver *ode_solver = NULL;
    switch (ode_solver_type)
    {
        // Explicit methods
        case 1: ode_solver = new ForwardEulerSolver; break;
        case 2: ode_solver = new RK2Solver(1.0); break;
        case 3: ode_solver = new RK3SSPSolver; break;
        case 4: ode_solver = new RK4Solver; break;
        case 6: ode_solver = new RK6Solver; break;
        default:
        if (Mpi::Root())
        {
           cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
        }
        delete mesh;
        return 5;
    }

    for (int lev = 0; lev < ser_ref_levels; lev++)
    {
        mesh->UniformRefinement();
    }
    if (mesh->NURBSext)
    {
        mesh->SetCurvature(max(order, 1));
    }

    mesh->GetBoundingBox(bb_min, bb_max, max(order, 1));

    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    for (int lev = 0; lev < par_ref_levels; lev++)
    {
        pmesh->UniformRefinement();
    }

    //*
    if(exec_mode == 0)
    {
        mesh_order = 1;
    }

    H1_FECollection mesh_fec(mesh_order, dim, BasisType::GaussLobatto);

    // Current mesh positions.
    ParFiniteElementSpace mesh_pfes(pmesh, &mesh_fec, dim);
    ParGridFunction x(&mesh_pfes);

    if(exec_mode == 1)
    {   
        pmesh->SetCurvature(mesh_order);
        pmesh->SetNodalGridFunction(&x);
    }

    // Store initial mesh positions.
    Vector x0 = x;

    ParGridFunction v_gf(x.ParFESpace());
    //VectorGridFunctionCoefficient v_mesh_coeff(&v_gf);
    if (exec_mode == 1)
    {
        ParGridFunction v(&mesh_pfes);
        VectorFunctionCoefficient v_coeff(dim, velocity_function);
        v.ProjectCoefficient(v_coeff);

        double t = 0.0;
        while (t < t_final)
        {
            // Move the mesh nodes.
            x.Add(min(dt, t_final-t), v);
            t += min(dt, t_final-t);

            // Update the node velocities.
            v.ProjectCoefficient(v_coeff);
        } 

        // Pseudotime velocity.
        add(x0, -1.0, x, v_gf);
        
        // Return the mesh to the initial configuration.
        x = x0;
        t_final = 1.0;
    }
    //*/

    H1_FECollection fec(order, dim, BasisType::Positive);
    ParFiniteElementSpace* pfes = new ParFiniteElementSpace(pmesh, &fec);

    HYPRE_BigInt global_vSize = pfes->GlobalTrueVSize();
    if (Mpi::Root())
    {
        cout << "Number of unknowns: " << global_vSize << endl;
    }

    FunctionCoefficient inflow(inflow_function);
    VectorFunctionCoefficient velocity(dim, velocity_function);
    FunctionCoefficient u0(u0_function);

    ParBilinearForm *mL = new ParBilinearForm(pfes);
    Vector lumpedmassmatrix(mL->Height());
    mL->AddDomainIntegrator(new LumpedIntegrator(new MassIntegrator));
    mL->Assemble();
    mL->Finalize();
    mL->SpMat().GetDiag(lumpedmassmatrix);

    ParBilinearForm *m = new ParBilinearForm(pfes);
    m->AddDomainIntegrator(new MassIntegrator);
    m->Assemble();
    m->Finalize();

    ParGridFunction u(pfes);
    u.ProjectCoefficient(u0);

    double loc_mass = u * lumpedmassmatrix;
    double glob_init_mass = 0.0;
    MPI_Allreduce(&loc_mass, &glob_init_mass, 1, MPI_DOUBLE, MPI_SUM,
              MPI_COMM_WORLD);

    {
        ostringstream mesh_name, sol_name;
        mesh_name << "output/mesh." << setfill('0') << setw(6) << myid;
        sol_name << "output/init." << setfill('0') << setw(6) << myid;
        ofstream omesh(mesh_name.str().c_str());
        omesh.precision(precision);
        pmesh->Print(omesh);
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u.Save(osol);
    }

    socketstream sout;
    if (visualization)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        sout.open(vishost, visport);
        if (!sout)
        {
            if (Mpi::Root())
            {
                cout << "Unable to connect to GLVis server at "
                    << vishost << ':' << visport << endl;
            }
            visualization = false;
            if (Mpi::Root())
            {
                cout << "GLVis visualization disabled.\n";
            }
        }
        else
        {
            sout << "parallel " << num_procs << " " << myid << "\n";
            sout.precision(precision);
            sout << "solution\n" << *pmesh << u;
            sout << "window_title '" << "advection" << "'\n"
              << "window_geometry "
              << 0 << " " << 0 << " " << 1080 << " " << 1080
              << "keys mcjl66666666666666666666666"
              << "66666666666666666666666666666666666666666666666662222222222";
              if(dim ==1)
              {
                sout << "RR";
              }
              else if(dim == 2 && problem != 2)
              {
                if(exec_mode == 1)
                {
                    sout << "PPPPPPPPPPPPPPPPl";
                }
                else
                {
                    sout << "R";
                }
              }
              sout << endl;

            //if (Mpi::Root())
            //{
            //    cout << "GLVis visualization paused."
            //        << " Press space (in the GLVis window) to resume it.\n";
            //}
        }
    }

    ParMesh *subcell_mesh = NULL;
    FiniteElementCollection *fec_sub = NULL;
    ParFiniteElementSpace *pfes_sub = NULL;
    ParFiniteElementSpace *dpfes_sub = NULL;
    ParGridFunction *xsub = NULL;
    ParGridFunction *v_sub_gf;
    VectorGridFunctionCoefficient v_sub_coef;
    Vector x0_sub;

    if(subcell)
    {
        const int btype = BasisType::ClosedUniform;
        subcell_mesh = new ParMesh(ParMesh::MakeRefined(*pmesh, order, btype));
        fec_sub = new H1_FECollection(1, dim, BasisType::ClosedUniform);
        pfes_sub = new ParFiniteElementSpace(subcell_mesh, fec_sub);

        if(exec_mode == 1)
        {
            dpfes_sub = new ParFiniteElementSpace(subcell_mesh, fec_sub, dim);
            xsub = new ParGridFunction(dpfes_sub);
            v_sub_gf = new ParGridFunction(dpfes_sub);
            VectorGridFunctionCoefficient v_sub_coef(&v_gf);
            v_sub_gf->ProjectCoefficient(v_sub_coef);
            subcell_mesh->SetCurvature(1);
            subcell_mesh->SetNodalGridFunction(xsub);
        }
    }

    FE_Evolution *met = NULL;  
    switch (scheme)
    {
        case 0: met = new LowOrderScheme(*pfes, inflow, velocity, *m, x0, v_gf, exec_mode);
            break;
        case 1: met = new ClipAndScale(*pfes, inflow, velocity, *m, x0, v_gf, exec_mode);
            break;  
        case 2: met = new Convex_ClipAndScale(*pfes, inflow, velocity, *m, x0, v_gf, exec_mode);
            break;
        case 3: met = new Subcell_LowOrder(*pfes, pfes_sub, inflow, velocity, *m, x0, v_gf, v_sub_gf, exec_mode);
            break;
        case 4: met = new Subcell_ClipAndScale(*pfes, pfes_sub, inflow, velocity, *m, x0, v_gf, v_sub_gf, exec_mode);
            break;
        default:
            MFEM_ABORT("Unkown scheme!");
    }

    double t = 0.0;
    met->SetTime(t);
    ode_solver->Init(*met);

    bool done = false;
    if(Mpi::Root())
    {
        cout << endl;
        cout << "Preprocessing done! Entering time loop!" << endl;
    }

    tic_toc.Clear();
    tic_toc.Start();
    for (int ti = 0; !done; )
    {
        double dt_real = min(dt, t_final - t);
        ode_solver->Step(u, t, dt_real);
        ti++;

        done = (t >= t_final - 1e-8*dt);

        if (done || ti % vis_steps == 0)
        {
            if (Mpi::Root())
            {
                cout << "time step: " << ti << ", time: " << t << endl;
            }
            if (visualization)
            {
                sout << "parallel " << num_procs << " " << myid << "\n";
                sout << "solution\n" << *pmesh << u << flush;
            }
        }
    }
    
    tic_toc.Stop();
    double min_loc = u.Min();
    double max_loc = u.Max();
    double min_glob, max_glob;

    MPI_Allreduce(&min_loc, &min_glob, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&max_loc, &max_glob, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    mL->BilinearForm::operator=(0.0);
    mL->Assemble();
    mL->SpMat().GetDiag(lumpedmassmatrix);

    loc_mass = u * lumpedmassmatrix;
    double glob_end_mass = 0.0;
    MPI_Allreduce(&loc_mass, &glob_end_mass, 1, MPI_DOUBLE, MPI_SUM,
              MPI_COMM_WORLD);

    Array <double> errors(3);
    errors[0] = u.ComputeL1Error(u0);
    errors[1] = u.ComputeL2Error(u0);
    errors[2] = u.ComputeMaxError(u0);

    if(Mpi::Root())
    {
        cout << "Time stepping loop done in " << tic_toc.RealTime() << " seconds."<< endl;
        cout << endl;

        cout << "L1 error:                           " << errors[0] << '\n';
        cout << "L2 error:                           " << errors[1] << '\n';
        cout << "Linf error                          " << errors[2] << "\n\n";

        cout << "Difference in solution mass: " << abs(glob_init_mass - glob_end_mass) << endl;
        cout << "u in [" << min_glob<< ", " << max_glob<< "]\n\n"; 
    }

    {
        ostringstream sol_name;
        sol_name << "output/final." << setfill('0') << setw(6) << myid;
        ofstream osol(sol_name.str().c_str());
        osol.precision(precision);
        u.Save(osol);
    }

    //delete u;
    delete pfes;
    delete pmesh;
    delete ode_solver;
    delete met;
    delete m;
    delete mL;

    return 0;
}


// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
    int dim = x.Size();

    // map to the reference [-1,1] domain
    Vector X(dim);
    for (int i = 0; i < dim; i++)
    {
        double center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
    }

    switch (problem)
    {
        case 0:
        case 1:
        {
            // Translations in 1D, 2D, and 3D
            switch (dim)
            {
                case 1: v(0) = 1.0; break;
                case 2: v(0) = sqrt(2./3.); v(1) = sqrt(1./3.); break;
                case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.);
                break;
            }
            break;
        }
        case 2:
        {
            v(0) = 2.0 * M_PI * (- X(1));
            v(1) = 2.0 * M_PI * (X(0) ); 
            if(dim > 2)
            {
                v(2) = 0.0;
            }
            break;
        }
        case 3:
        {
            // Clockwise rotation in 2D around the origin
            const real_t w = M_PI/2;
            switch (dim)
            {
                case 1: v(0) = 1.0; break;
                case 2: v(0) = w*X(1); v(1) = -w*X(0); break;
                case 3: v(0) = w*X(1); v(1) = -w*X(0); v(2) = 0.0; break;
            }
         break;
            break;
        }
        case 4:
        {
            // Clockwise twisting rotation in 2D around the origin
            const double w = M_PI/2;
            double d = max((X(0)+1.)*(1.-X(0)),0.) * max((X(1)+1.)*(1.-X(1)),0.);
            d = d*d;
            switch (dim)
            {
                case 1: v(0) = 1.0; break;
                case 2: v(0) = d*w*X(1); v(1) = -d*w*X(0); break;
                case 3: v(0) = d*w*X(1); v(1) = -d*w*X(0); v(2) = 0.0; break;
            }
            break;
        }
        case 12:
        {
            if (dim != 2) { MFEM_ABORT("Not implemented."); }
            for (int d = 0; d < dim; d++) { X(d) = X(d) * 0.5 + 0.5; }
            v(0) =  sin(M_PI*X(0)) * cos(M_PI*X(1));
            v(1) = -cos(M_PI*X(0)) * sin(M_PI*X(1));
            break;
        }
        case 13:
        {
            // Gresho deformation used for mesh motion in remap tests.
            const double r = sqrt(X(0)*X(0) + X(1)*X(1));
            if (r < 0.2)
            {
                v(0) =  5.0 * X(1);
                v(1) = -5.0 * X(0);
            }
            else if (r < 0.4)
            {
                v(0) =  2.0 * X(1) / r - 5.0 * X(1);
                v(1) = -2.0 * X(0) / r + 5.0 * X(0);
            }
            else { v = 0.0; }
            break;
        }
        default:
            MFEM_ABORT("No velocity function for this problem")
    }
}

// Initial condition
double u0_function(const Vector &x)
{
    int dim = x.Size();

    // map to the reference [-1,1] domain
    Vector X(dim);
    for (int i = 0; i < dim; i++)
    {
        double center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
    }

    int problem_num = problem % 10;

    switch (problem_num)
    {
        case 0:
        {
            if(dim != 1)
            {
                MFEM_ABORT("Problem 0 only works in 1D!");
            }
         
            if(X(0) < 0.9 && X(0) > 0.5)
            {
                return exp(10) * exp(1.0 / (0.5 - X(0))) * exp(1.0 / (X(0) - 0.9));
            }
            else if( X(0) < 0.4 && X(0) > 0.2)
            {
                return 1.0;
            }
            else 
            {
                return 0.0;
            }
            break;
        }
        case 1:
        {
            if (dim != 1)
            {
                MFEM_ABORT("Problem 1 only works in 1D.");
            }
            else
            {
                double a = 0.5;
                double z = -0.7;
                double delta = 0.005;
                double alpha = 10.0;
                double beta = log(2) / 36 / delta / delta;
                if(X(0) >= -0.8 && X(0) <= -0.6)
                {
                    double G1 = exp(-beta * (X(0) - (z - delta)) * (X(0) - (z - delta) ));
                    double G2 = exp(-beta * (X(0) - (z + delta)) * (X(0) - (z + delta) ));
                    double G3 = exp(-beta * (X(0) - z) * (X(0) - z ));
                    return 1.0 / 6.0 * ( G1 + G2 + 4.0 * G3);
                }
                else if(X(0) >= -0.4 && X(0) <= -0.2)
                {
                    return 1.0;
                }
                else if(X(0) >= 0.0 && X(0) <= 0.2)
                {
                    return 1.0 - abs(10.0 * (X(0) - 0.1));
                }
                else if(X(0) >= 0.4 && X(0) <= 0.6)
                {
                    double F1 = sqrt( max(1.0 - alpha * alpha * (X(0) - (a - delta)) *  (X(0) - (a - delta)), 0.0));
                    double F2 = sqrt( max(1.0 - alpha * alpha * (X(0) - (a + delta)) *  (X(0) - (a + delta)), 0.0));
                    double F3 = sqrt( max(1.0 - alpha * alpha * (X(0) - a) *  (X(0) - a), 0.0));
                    return 1.0 / 6.0 * ( F1 + F2 + 4.0 * F3);   
                }
                else
                {
                    return 0.0;
                }
            }
            break;
        }
        case 2:
        {
            if (dim != 2) 
            { 
                //MFEM_ABORT("Solid body rotation only works in 2 D."); 
            }
         
            // Initial condition defined on [0,1]^2
            Vector y = X;
            y *= 0.5;
            y += 0.5;
            double s = 0.15;
            double cone = sqrt(pow(y(0) - 0.5, 2.0) + pow(y(1) - 0.25, 2.0));
            double hump = sqrt(pow(y(0) - 0.25, 2.0) + pow(y(1) - 0.5, 2.0));
            return (1.0 - cone / s) * (cone <= s) + 0.25 * (1.0 + cos(M_PI*hump / s)) * (hump <= s) +
                ((sqrt(pow(y(0) - 0.5, 2.0) + pow(y(1) - 0.75, 2.0)) <= s ) && ( abs(y(0) -0.5) >= 0.025 || (y(1) >= 0.85) ) ? 1.0 : 0.0);
            break;
        }
        case 3:
        {
            switch (dim)
            {
                case 1:
                    return exp(-40.*pow(X(0)-0.5,2));
                case 2:
                case 3:
                {
                    double rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
                    if (dim == 3)
                    {
                        const double s = (1. + 0.25*cos(2*M_PI*X(2)));
                        rx *= s;
                        ry *= s;
                    }
                    return ( std::erfc(w*(X(0)-cx-rx))*std::erfc(-w*(X(0)-cx+rx)) *
                        std::erfc(w*(X(1)-cy-ry))*std::erfc(-w*(X(1)-cy+ry)) )/16;
                }
            }
            break;
        }
        case 4:
        {
            double x_ = X(0), y_ = X(1), rho, phi;
            rho = std::hypot(x_, y_);
            phi = atan2(y_, x_);
            return pow(sin(M_PI*rho),2)*sin(3*phi);
            break;
        }
        case 5:
        {
            const double f = M_PI;
            return sin(f*X(0))*sin(f*X(1));
            break;
        }
    }
    return 0.0;
}


/*
// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
    int dim = x.Size();

    // map to the reference [-1,1] domain
    Vector X(dim);
    for (int i = 0; i < dim; i++)
    {
        double center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
    }

    switch (problem)
    {
        case 0:
        {
            // Translations in 1D, 2D, and 3D
            switch (dim)
            {
                case 1:
                    v(0) = 1.0;
                    break;
                case 2:
                    v(0) = sqrt(2./3.);
                    v(1) = sqrt(1./3.);
                    break;
                case 3:
                    v(0) = sqrt(3./6.);
                    v(1) = sqrt(2./6.);
                    v(2) = sqrt(1./6.);
                    break;
            }
            break;
        }
        case 1:
        case 2:
        {
            // Clockwise rotation in 2D around the origin
            const double w = M_PI/2;
            switch (dim)
            {
                case 1:
                    v(0) = 1.0;
                    break;
                case 2:
                    v(0) = w*X(1);
                    v(1) = -w*X(0);
                    break;
                case 3:
                    v(0) = w*X(1);
                    v(1) = -w*X(0);
                    v(2) = 0.0;
                    break;
            }
            break;
        }
        case 3:
        {
            // Clockwise twisting rotation in 2D around the origin
            const double w = M_PI/2;
            double d = max((X(0)+1.)*(1.-X(0)),0.) * max((X(1)+1.)*(1.-X(1)),0.);
            d = d*d;
            switch (dim)
            {
                case 1:
                    v(0) = 1.0;
                    break;
                case 2:
                    v(0) = d*w*X(1);
                    v(1) = -d*w*X(0);
                    break;
                case 3:
                    v(0) = d*w*X(1);
                    v(1) = -d*w*X(0);
                    v(2) = 0.0;
                    break;
            }
            break;
        }
    }
}


// Initial condition
double u0_function(const Vector &x)
{
    int dim = x.Size();

    // map to the reference [-1,1] domain
    Vector X(dim);
    for (int i = 0; i < dim; i++)
    {
        double center = (bb_min[i] + bb_max[i]) * 0.5;
        X(i) = 2 * (x(i) - center) / (bb_max[i] - bb_min[i]);
    }

    switch (problem)
    {
        case 0:
        case 1:
        {
            switch (dim)
            {
                case 1:
                    return exp(-40.*pow(X(0)-0.5,2));
                case 2:
                case 3:
                {
                    double rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
                    if (dim == 3)
                    {
                        const double s = (1. + 0.25*cos(2*M_PI*X(2)));
                        rx *= s;
                        ry *= s;
                    }
                    return ( std::erfc(w*(X(0)-cx-rx))*std::erfc(-w*(X(0)-cx+rx)) *
                            std::erfc(w*(X(1)-cy-ry))*std::erfc(-w*(X(1)-cy+ry)) )/16.0;
                }
            }
        }
        case 2:
        {
            double x_ = X(0), y_ = X(1), rho, phi;
            rho = std::hypot(x_, y_);
            phi = atan2(y_, x_);
            return pow(sin(M_PI*rho),2)*sin(3*phi);
        }
        case 3:
        {
            const double f = M_PI;
            return sin(f*X(0))*sin(f*X(1));
        }
    }
    return 0.0;
} //*/

// Inflow boundary condition (0.0 for the problems considered in this example)
double inflow_function(const Vector &x)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3:
         return 0.0;
   }
   return 0.0;
}
