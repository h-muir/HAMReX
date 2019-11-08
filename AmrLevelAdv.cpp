
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include "structdefs.H"
#include "funcdefs.H"

using namespace amrex;
using namespace std;

int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9;
int      AmrLevelAdv::do_reflux       = 1;

int      AmrLevelAdv::NUM_STATE       = 8;  // One variable in the state
int      AmrLevelAdv::NUM_GROW        = 2;  // number of ghost cells

ParameterStruct  AmrLevelAdv::Parameters;
SettingsStruct  AmrLevelAdv::SimSettings;

Vector<BCRec> AmrLevelAdv::bc_rec(NUM_STATE);

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
    flux_reg = 0;
}

//
//The basic constructor.
//
AmrLevelAdv::AmrLevelAdv (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

//
//The destructor.
//
AmrLevelAdv::~AmrLevelAdv () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
void
AmrLevelAdv::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

void 
AmrLevelAdv::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

//
//Write a plotfile to specified directory.
//

void
AmrLevelAdv::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
	//variable types specified in AmrLevelAdv::variableSetUp () below
    AmrLevel::writePlotFile(dir, os, how);
 
}


//
//Define data descriptors.
//
void
AmrLevelAdv::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
			   &cell_cons_interp);

	
    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	lo_bc[i] = hi_bc[i] = BCType::foextrap; //trying to make transmissive: BCType::foextrap  
											//does that require ghost cells on Snew ?
											// periodic boundaries: BCType::int_dir;
    }
    
    BCRec bc(lo_bc, hi_bc);
    
    /* H-----------------------------------------------------
     * SimSettings and Parameters are member variables of 
     * AmrLevelAdv, of struct types: SettingsStruct and 
     * Parameter struct respectively, defined in StructDefs
     * -----------------------------------------------------*/ 
    
    initialiseStructs(SimSettings, Parameters);
    
    if(NUM_STATE != SimSettings.NCOMP){
		amrex::Abort("NUM_STATE != NCOMP in input file");
	}
    
    Vector<std::string> variable_names = Parameters.problem_variables;
    
	for(int comp = 0; comp < NUM_STATE; ++comp){
		amrex::Print() << "variable name " << comp << ": " << variable_names[comp] << "\n";
		desc_lst.setComponent(Phi_Type, comp, variable_names[comp], bc, 
			  StateDescriptor::BndryFunc(phifill));
	}
}

//
//Cleanup data descriptors at end of run.
//
void
AmrLevelAdv::variableCleanUp () 
{
    desc_lst.clear();
}

//
//Initialize grid data at problem start-up.
//
void
AmrLevelAdv::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();
	

    amrex::Print() << "nsteps_max: " << Parameters.nsteps_max << "\n";
    
    
    AccessVariable V(Parameters.problem_variables);
    /*
    cout << "V[rho] : " << V["rho"] << endl;
    
    for(int i = 0; i<NUM_STATE; ++i){
		cout << V[Parameters.problem_variables[i]] << " ";
	}
	cout << endl;
    */
    
    //these cell sizes change with each re-initialisation
    // under refinement
    //Parameters.dx[0] = dx[0];
    //Parameters.dx[1] = dx[1];
    

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        FArrayBox& fab 	   = S_new[mfi];
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        
        Array4<Real> const& prop_arr = fab.array();	
		
		initial_conditions(box, prop_arr, dx, Parameters, V);
		
		   
    }

    if (verbose) {
	amrex::Print() << "Done initializing the level " << level 
                       << " data " << std::endl; 
    }
    
    //TESTING:
    /*
    Vector<Real> v0(3);
	Vector<Real> v1(3,1.0);
	Vector<Real> vR = v0 + v1;
	amrex::Print() << "check addition: \n";
	printVector(vR);
	vR = v0 - v1;
	amrex::Print() << "check subtraction: \n";
	printVector(vR);
	vR = 2.0*v1;
	amrex::Print() << "check multiplication: \n";
	printVector(vR);
	*/
}

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
//
void
AmrLevelAdv::init (AmrLevel &old)
{
    AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{

    MultiFab& S_mm = get_new_data(Phi_Type);
    Real maxval = S_mm.max(0);
    Real minval = S_mm.min(0);
    
    amrex::Print() << "phi max = " << maxval << ", min = " << minval  << std::endl;
    /*
    //Q:not sure if this is needed:
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }
	*/
    MultiFab& S_new = get_new_data(Phi_Type);
    
    MultiFab S_old(S_new.boxArray(), dmap, S_new.nComp(), S_new.nGrow()+1); //H
    FillPatch(*this, S_old, S_new.nGrow()+1, time, Phi_Type, 0, NUM_STATE); //H

    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();


    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
		fine = &getFluxReg(level+1);
		fine->setVal(0.0);
    }

    if (do_reflux && level > 0) {
		current = &getFluxReg(level);
    }

    MultiFab fluxes[BL_SPACEDIM]; //Q: does this make an array of MultiFabs, of length BL_SPACEDIM,
									//like python syntax? 

    if (do_reflux)
    {
	for (int j = 0; j < BL_SPACEDIM; j++)
	{
	    BoxArray ba = S_new.boxArray();
	    ba.surroundingNodes(j);
	    fluxes[j].define(ba, dmap, NUM_STATE, 0);
	}
    }
	
#ifdef _OPENMP
#pragma omp parallel
#endif
    
    
    //H------------------------------------------------------------------------- 
    //create a multifab array to store the fluxes:
	MultiFab flux_arr[BL_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = S_new.boxArray();
        edge_ba.surroundingNodes(dir); //switches box from type:cell to type:node in direction dir ? 
        flux_arr[dir].define(edge_ba, dmap, NUM_STATE, 0); //ngrow=0

    }
    //------------------------------------------------------------------------H/
    
    {
	//FArrayBox flux[BL_SPACEDIM];
	//H-------------------------------------------------------------------------
	Real dxr = dx[0];
	Real dyr = dx[1];
	AccessVariable V(Parameters.problem_variables);
	//------------------------------------------------------------------------H/
	
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& statein 		 = 	 S_old[mfi];	//H
	    FArrayBox& flux_fab_x 	 = 	 flux_arr[0][mfi];	//H
	    //--------------------------------------------------------------------H/		

		//H --------------------------------------------------------------------
		//access the array from the array box:
		
        Array4<Real> const& S_old_arr = statein.array();
        Array4<Real> const& S_new_arr = stateout.array();
        Array4<Real> const& flux_prop_arr_x = flux_fab_x.array();
		
		calc_fluxes_HLLC_x(bx, S_old_arr, flux_prop_arr_x, dt, dxr, NUM_STATE, V, SimSettings.MUSCL);
        cons_update("x-sweep", bx, flux_prop_arr_x, S_old_arr, S_new_arr, dt, dxr, NUM_STATE, V);
        
        //--------------------------------------------------------------------H/
		    		
	    if (do_reflux) {
			rescale_fluxes("x-sweep", bx, flux_prop_arr_x, dt, dyr, NUM_STATE, V); //rescale by dt*dy for Fx
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			fluxes[0][mfi].copy(flux_fab_x,mfi.nodaltilebox(0));
			//fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));
	    }
	    
	} //end x-sweepMFIter	
	
    //copy data from updated S_new --> S_old between dimensional sweeps:
	MultiFab& S_new = get_new_data(Phi_Type);
	MultiFab S_old(S_new.boxArray(), dmap, S_new.nComp(), S_new.nGrow()+1); //H
    FillPatch(*this, S_old, S_new.nGrow()+1, time, Phi_Type, 0, NUM_STATE);
    
	//y-sweep:
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& statein 		 = 	 S_old[mfi];	//H
	    FArrayBox& flux_fab_y 	 =   flux_arr[1][mfi];	//H
	    //H --------------------------------------------------------------------
		
		//H --------------------------------------------------------------------
		//access the array from the array box:
		
        Array4<Real> const& S_old_arr = statein.array();
        Array4<Real> const& S_new_arr = stateout.array();
        Array4<Real> const& flux_prop_arr_y = flux_fab_y.array();
		
		
		calc_fluxes_HLLC_y(bx, S_old_arr, flux_prop_arr_y, dt, dyr, NUM_STATE, V, SimSettings.MUSCL);
        cons_update("y-sweep", bx, flux_prop_arr_y, S_old_arr, S_new_arr, dt, dyr, NUM_STATE, V);
        
        //----------------------------------------------------------------------

		    		
	    if (do_reflux) {
			rescale_fluxes("y-sweep", bx, flux_prop_arr_y, dt, dxr, NUM_STATE, V); //rescale by dt*dx for Fy
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			fluxes[1][mfi].copy(flux_fab_y,mfi.nodaltilebox(1));
			//fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));
	    }
	    
	} //end y-sweep MFIter	
	 		 	 
    } //end local scope
	
    if (do_reflux) {
		
		if (current) {
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			current->FineAdd(fluxes[0],0,0,0,NUM_STATE,1.);
			current->FineAdd(fluxes[1],1,0,0,NUM_STATE,1.);
		}
		
		if (fine) {
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			fine->CrseInit(fluxes[0],0,0,0,NUM_STATE,-1.);
			fine->CrseInit(fluxes[1],1,0,0,NUM_STATE,-1.);
		}
    }
	

    return dt;
}

//
//Estimate time step.
//
Real
AmrLevelAdv::estTimeStep (Real)
{
    // This is just a dummy value to start with 
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real cur_time = state[Phi_Type].curTime();
    MultiFab& S_new = get_new_data(Phi_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	AccessVariable V(Parameters.problem_variables);
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		
		const Box& bx = mfi.tilebox();
	    FArrayBox& stateout     	  =   		S_new[mfi];
	    Array4<Real> const& S_new_arr =   stateout.array();
		
	    for (int i = 0; i < 2; ++i) { //configured for x and y directions only
			Real umax = max_dim_wave_speed(bx, S_new_arr, i, V); //uface[i].norm(0);
			if (umax > 1.e-100) {
				dt_est = std::min(dt_est, dx[i] / umax);
			}
	    }
	}
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    dt_est *= cfl;

    if (verbose) {
	amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
                       << ":  dt_est = " << dt_est << std::endl;
    }
    
    return dt_est;
}

//
//Compute initial time step.
//
Real
AmrLevelAdv::initialTimeStep ()
{
    return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
							   int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Compute new `dt'.
//
void
AmrLevelAdv::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        AmrLevelAdv& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1) 
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    }
    else 
    {
	//
	// Limit dt's by change_max * old dt
	//
	static Real change_max = 1.1;
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	}
    }
    
    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Do work after timestep().
//
void
AmrLevelAdv::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();
}

//
//Do work after regrid().
//
void
AmrLevelAdv::post_regrid (int lbase, int new_finest) {
//particle function (removed from here)
}

//
//Do work after a restart().
//
void
AmrLevelAdv::post_restart() 
{
//particle function (removed from here)
}

//
//Do work after init().
//
void
AmrLevelAdv::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//
void
AmrLevelAdv::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);
    
    //with bordering ghost cells:
    MultiFab S_old(S_new.boxArray(), dmap, S_new.nComp(), S_new.nGrow()+1); //H
    FillPatch(*this, S_old, S_new.nGrow()+1, time, Phi_Type, 0, NUM_STATE); //H
    
    //cout << "level: " << level << endl;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	AccessVariable V(Parameters.problem_variables);
	string	prop = Parameters.refinement_prop;		//refinement based off gradient of property: mass
	string refinement_condition = Parameters.refinement_condition;
	Real grad_frac = Parameters.refinement_grad_fracs[level];
	Real grad_max = 0;
	//determine max rho gradient:
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    const Box&  bx  				= mfi.tilebox();
	    
	    FArrayBox& fab  				= S_old[mfi];
	    
	    Array4<Real> const& S_old_arr 	= fab.array();

        TagBox& tagfab  				= tags[mfi]; //tags is passed to ErrorEst as argument
        
        Array4<char> const& tagarr 		= tagfab.array();
        
        determine_grad_max(bx,S_old_arr,dx, prop, refinement_condition, grad_max, V);

	}
	ParallelDescriptor::ReduceRealMax(grad_max);	
	//cout << "grad_max = " << grad_max << endl;
	//tagging based on propertion of max_rho
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    const Box&  bx  				= mfi.tilebox();
	    
	    FArrayBox& fab  				= S_old[mfi];
	    
	    Array4<Real> const& S_old_arr 	= fab.array();

        TagBox& tagfab  				= tags[mfi]; //creates tagbox with MFI dimensions ?
        
        Array4<char> const& tagarr 		= tagfab.array();
        
        amr_tagging(tagarr,bx,S_old_arr,dx,grad_max, grad_frac, prop, V);
	    
	    
	}
    }
}

//can probably remove this whole function:
void
AmrLevelAdv::read_params ()
{
	
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");   

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);
	
	/*
    // This tutorial code only supports Cartesian coordinates.
    if (! Geometry::IsCartesian()) {
	amrex::Abort("Please set geom.coord_sys = 0");
    }
	
	
    // This tutorial code only supports periodic boundaries.
    if (! Geometry::isAllPeriodic()) {
	amrex::Abort("Please set geom.is_periodic = 1 1 1");
    }
	*/
    //
    // read tagging Parameters from probin file
    //

    std::string probin_file("probin");

    ParmParse ppa("amr");
    ppa.query("probin_file",probin_file);

    int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

    // use a fortran routine to
    // read in tagging parameters from probin file
    get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

}

void
AmrLevelAdv::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);
    
    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = amrex::second() - strt;
	
        ParallelDescriptor::ReduceRealMax(end,IOProc);
	
        amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
                       << " : time = " << end << std::endl;
    }
}

void
AmrLevelAdv::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
AmrLevelAdv::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    AmrLevelAdv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

