#ifndef STRUCTDEFS_H_
#define STRUCTDEFS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <string>

using namespace amrex;

struct SettingsStruct
{
	public:
		 int nsteps_max;
		 int plot_int;
		// int max_grid_size;
		 int ndmin;
		 Real startT;
		 Real finalT;
		 bool chatter;
		 std::string testcase;
		 int NCOMP;
		 bool MUSCL;
		 bool rigid_body;
		 std::string rigid_geom;
			 
};

struct ParameterStruct  : public SettingsStruct
{
	public:
		 Vector<Real> prob_lo;
		 Vector<Real> prob_hi;
		 Vector<int>  n_cells;
		 Vector<int> is_periodic; 	//boundary conditions
		 Vector<Real> dx;
		
		 Real coord_sys; 			//cartesian, cylindrical or pherical
		 Vector<Real> PL;
		 Vector<Real> PR;
		 Real x0;
		 int n_vars;
		 Vector<std::string> vars;		
		 std::string IC;
		 int max_level;
		 std::string refinement_prop;
		 std::string refinement_condition;
		 Vector<Real> refinement_grad_fracs;
		 
		 Vector<std::string> problem_variables;
		 
		 //rigid body parameters:
		 Real radius;
		 Vector<Real> centre;
		 

    
    ParameterStruct()
    {
		prob_lo 	= Vector<Real>(AMREX_SPACEDIM);
		prob_hi 	= Vector<Real>(AMREX_SPACEDIM);
		n_cells 	= Vector<int>(AMREX_SPACEDIM);
		is_periodic = Vector<int>(AMREX_SPACEDIM);
		dx 		 	= Vector<Real>(AMREX_SPACEDIM);
		vars		= Vector<std::string>(n_vars);
		PL		 	= Vector<Real>(n_vars);
		PR		 	= Vector<Real>(n_vars);
		
		refinement_grad_fracs = Vector<Real>(max_level);
		
		problem_variables = Vector<std::string>(NCOMP);
		
		centre 		 = Vector<Real>(AMREX_SPACEDIM);
		
	}
};





#endif
