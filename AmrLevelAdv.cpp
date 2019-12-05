
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_Level.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_levelset.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_LSCoreBase.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBMultiFabUtil.H>

#include "structdefs.H"
#include "funcdefs.H"

using namespace amrex;
using namespace std;

int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9;
int      AmrLevelAdv::do_reflux       = 1;

int      AmrLevelAdv::NUM_STATE       = 13;  // One variable in the state
int      AmrLevelAdv::NUM_GROW        = 2;  // number of ghost cells

ParameterStruct  AmrLevelAdv::Parameters;
SettingsStruct  AmrLevelAdv::SimSettings;

//string 	 AmrLevelAdv::EoS; //specify EoS as either: "ideal" or "plasma19"
//Plasma19 AmrLevelAdv::AirPlasma("mixture19_cns.txt");


Vector<BCRec> AmrLevelAdv::bc_vec(NUM_STATE);

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
	cout << "begin AmrLevel constructor 1" << endl;
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
/*
void
AmrLevelAdv::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
	//variable types specified in AmrLevelAdv::variableSetUp () below
    AmrLevel::writePlotFile(dir, os, how);
 
}
*/
void
AmrLevelAdv::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
    BL_PROFILE("CNS::writePlotFile()");

//    AmrLevel::writePlotFile(dir, os, how);
    
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++) {
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++) {
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType()) {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
            }
        }
    }

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
	 it != dlist.end();
	 ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
        {
            derive_names.push_back(it->name());
            num_derive++;
	}
    }

    int n_data_items = plot_var_map.size() + num_derive + 1;

    Real cur_time = state[Phi_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            amrex::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (int i = 0; i < plot_var_map.size(); i++)
        {
	    int typ = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for ( std::list<std::string>::iterator it = derive_names.begin();
	      it != derive_names.end(); ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
            os << rec->variableName(0) << '\n';
        }

        // volfrac
        os << "vfrac\n";

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (int i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string sLevel = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
        FullPath += '/';
    FullPath += sLevel;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(FullPath, 0755))
            amrex::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = sLevel;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }

        // volfrac threshhold for amrvis
        if (level == parent->finestLevel()) {
            for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                os << "1.0e-6\n";
            }
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow, MFInfo(), Factory());
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (int i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
#ifdef BL_TESTING
        // to avoid fcompare failure
        if (typ == Cost_Type) {
            plotMF.setVal(0.0, cnt, 1, nGrow);
        }
#endif
	cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::iterator it = derive_names.begin();
	     it != derive_names.end(); ++it)
	{
            auto derive_dat = derive(*it,cur_time,nGrow);
            MultiFab::Copy(plotMF,*derive_dat,0,cnt,1,nGrow);
	    cnt++;
	}
    }

    plotMF.setVal(0.0, cnt, 1, nGrow);

    //MultiFab::Copy(plotMF,volFrac(),0,cnt,1,nGrow);

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
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
    
    /* H-----------------------------------------------------
     * SimSettings and Parameters are member variables of 
     * AmrLevelAdv, of struct types: SettingsStruct and 
     * Parameter struct respectively, defined in StructDefs
     * -----------------------------------------------------*/ 
    
    initialiseStructs(SimSettings, Parameters);
    
    AccessVariable V(Parameters.problem_variables);
    
    setBoundaryConditions(bc_vec, Parameters, NUM_STATE, V);
    
    if(Parameters.IC == "mach_theta"){
		if(Parameters.vars[0] == "rho" && Parameters.vars[1] == "u" && \
		   Parameters.vars[2] == "v"   && Parameters.vars[3] == "p"){
			Vector<Real> u_vec{Parameters.PR[1], Parameters.PR[2], 0};
			Prim WR(Parameters.PR[0], u_vec, Parameters.PR[3]);
			Vector<Real> normal_vec{cos(Parameters.theta),-sin(Parameters.theta),0};
			Prim WL = rankine_hugoniot(WR, Parameters.mach, normal_vec);
			Parameters.PL[0] = WL.rho;
			Parameters.PL[1] = WL.u_vec[0]; 
			Parameters.PL[2] = WL.u_vec[1];
			Parameters.PL[3] = WL.p;
		}else{
			amrex::Abort("P vars input properties not input as assumed");
		}
	}
    
    if(NUM_STATE != SimSettings.NCOMP){
		amrex::Abort("NUM_STATE != NCOMP in input file");
	}
    
    Vector<std::string> variable_names = Parameters.problem_variables;
    
    StateDescriptor::BndryFunc 	func;
    func = StateDescriptor::BndryFunc(phifill);
    
    /*
    if(SimSettings.testcase == "inclined_place"){
		func = StateDescriptor::BndryFunc(mixedfill);
	}else{
		func = StateDescriptor::BndryFunc(phifill);
	}
    */
    
	for(int comp = 0; comp < NUM_STATE; ++comp){
		amrex::Print() << "variable name " << comp << ": " << variable_names[comp] << "\n";
		desc_lst.setComponent(Phi_Type, comp, variable_names[comp], bc_vec[comp], 
			  func);
	}
	
	cout << "end of variableSetUp()" << endl;
	
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
    
    //create new multifab with ngrow (for finite difference routine to operate on):
	MultiFab S_old(S_new.boxArray(), dmap, S_new.nComp(), NUM_GROW); //H
    FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE); //H
	

    amrex::Print() << "nsteps_max: " << Parameters.nsteps_max << "\n";
    
    AccessVariable V(Parameters.problem_variables);

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }
	
	//EB:
	const Vector<int> ngrow_vec{NUM_GROW,NUM_GROW,NUM_GROW};
	std::unique_ptr<amrex::EBFArrayBoxFactory> factory_ptr = makeEBFabFactory(geom, S_new.boxArray(), dmap, ngrow_vec, EBSupport::full);
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(S_new.Factory()));  
    const FabArray<EBCellFlagFab>& flags = factory->getMultiEBCellFlagFab();
    MultiCutFab const& centroid = factory->getCentroid();
	
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        FArrayBox& fab 	   = S_new[mfi];
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        
        Array4<Real> const& S_new_arr = fab.array();
        
        FabType EB_fab_type		 	= 	flags[mfi].getType(box); //EB
	    auto const& flags_bx		= 	flags[mfi];
		Array4<const EBCellFlag> flags_arr	= flags_bx.array();
		
		initial_conditions(NUM_STATE, box, S_new_arr, dx, Parameters, V);
		
		if(EB_fab_type == FabType::covered) {
			//cout << "box is covered \n";
			S_new[mfi].setVal(-1.0, box, V["EB"], 1);
		}else if(EB_fab_type == FabType::regular){
			//cout << "box is regular \n";
			S_new[mfi].setVal(1.0, box, V["EB"], 1);
		}else if(EB_fab_type == FabType::singlevalued){
			//cout << "box has cut cells: \n";
			find_embedded_boundary(box, S_new_arr, flags_arr, NUM_STATE, V);
			//able to get cutFab:
			const auto& centroid_fab = centroid[mfi];
		}
		  
    }
	
	initialise_levelset_geometry(level, geom, S_new.boxArray(), S_new.DistributionMap(), SimSettings, Parameters, V, S_new);
	
	for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     		= mfi.validbox();
        FArrayBox& fab_new 	  	= S_new[mfi];
	    FArrayBox& fab_old 		= S_old[mfi];	
		
        Array4<Real> const& S_old_arr = fab_new.array();
        Array4<Real> const& S_new_arr = fab_old.array();

		initialise_levelset_normals(box, S_old_arr, S_new_arr, dx, V);

	}//end: MFIter
	

    if (verbose) {
	amrex::Print() << "Done initializing the level " << level 
                       << " data " << std::endl; 
    }

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
    
    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);
    
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    
    AccessVariable V(Parameters.problem_variables);

    MultiFab& S_new = get_new_data(Phi_Type);
    
    MultiFab S_old(S_new.boxArray(), dmap, S_new.nComp(), NUM_GROW); //H
    FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE); //H
    const Box& prob_domain = geom.Domain();
    special_boundary(prob_domain, S_old, SimSettings, dx, V);
    
    
    //FillDomainBoundary(S_old, geom, bc_vec);
    //S_old.FillBoundary(geom.periodicity());

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
    
    //------------------------EMBEDDED BOUNDARY--------------------------------
	
	const Vector<int> ngrow_vec{NUM_GROW,NUM_GROW,NUM_GROW};
	std::unique_ptr<amrex::EBFArrayBoxFactory> factory_ptr = makeEBFabFactory(geom, S_new.boxArray(), dmap, ngrow_vec, EBSupport::full);
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(S_new.Factory()));  
    const FabArray<EBCellFlagFab>& flags = factory->getMultiEBCellFlagFab();
    MultiCutFab const& centroid = factory->getCentroid();
    
    /*
    if(Parameters.EBgeom != "none"){
		//set up level set multifabs
		int ls_pad = 1;
		MultiFab LSMF(S_new.boxArray(), dmap, 1, ls_pad);			//ncomp = 1
		iMultiFab ivalid_LSMF(S_new.boxArray(), dmap, 1, ls_pad);
		//fills level set multifabs: LSMF & ivalid_LSMF:
		constructLevelSet(level, geom, S_new.boxArray(), dmap, SimSettings, Parameters, LSMF, ivalid_LSMF);
		//0th element is the level set data in LSMF -- > 1 component copied to level_set index of S_new:
		MultiFab::Copy(S_new, LSMF, 0, V["level_set"], 1, S_new.nGrow());	
		//FillPatch(*this, S_new, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
		
    }
    */
    
    //-------------------------------------------------------------------------
    
    {
	//FArrayBox flux[BL_SPACEDIM];
	//H-------------------------------------------------------------------------
	Real dxr = dx[0];
	Real dyr = dx[1];
	
	//HALF TIME STEP RADIAL SOURCE TERM INTEGRATION:------------------------------------------------------
	if(Parameters.coord_sys == 1){
		for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
		{
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& statein 		 = 	 S_old[mfi];	//H
	    //--------------------------------------------------------------------H/		

		//H --------------------------------------------------------------------
		//access the array from the array box:
		
        Array4<Real> const& S_old_arr = statein.array();
        Array4<Real> const& S_new_arr = stateout.array();
		
		radial_sourceterm_integration(bx, S_old_arr, S_new_arr, dt/2.0, dxr, V);
		}
		
		MultiFab::Copy(S_old, S_new, 0, 0, NUM_STATE, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V);
	}
	
	//----------------------------------------------------------------------------------------------------
	//const EBCellFlagFab& flag_fab_x = flags[0];
	//EBFArrayBox flag_fab_arr  	= 	EBFArrayBox(flag_fab_x, bx, NUM_STATE);		//EB
	//x-sweep:--------------------------------------------------------------------------
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& statein 		 = 	 S_old[mfi];	
	    FArrayBox& flux_fab_x 	 = 	 flux_arr[0][mfi];	
	    
	    FabType EB_fab_type		 	= 	flags[mfi].getType(bx); //EB
	    auto const& flags_bx		= 	flags[mfi];
	    
	    //--------------------------------------------------------------------H/		

		//H --------------------------------------------------------------------
		//access the array from the array box:
		
        Array4<Real> const& S_old_arr 		= statein.array();
        Array4<Real> const& S_new_arr 		= stateout.array();
        Array4<Real> const& flux_prop_arr_x = flux_fab_x.array();
        
        Array4<const EBCellFlag> flags_arr	= flags_bx.array();
		
		if(EB_fab_type == FabType::covered) {
			//cout << "box is covered \n";
			S_new[mfi].setVal(-1.0, bx, V["EB"], 1);
		}else if(EB_fab_type == FabType::regular){
			//cout << "box is regular \n";
			S_new[mfi].setVal(1.0, bx, V["EB"], 1);
		}else if(EB_fab_type == FabType::singlevalued){
			//cout << "box has cut cells: \n";
			find_embedded_boundary(bx, S_new_arr, flags_arr, NUM_STATE, V); //EB val set in here
			//able to get cutFab:
			const auto& centroid_fab = centroid[mfi];
		}
		
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
	MultiFab::Copy(S_old, S_new, 0, 0, NUM_STATE, S_new.nGrow());
    FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    special_boundary(prob_domain, S_old, SimSettings, dx, V);
    //note: S_new.nGrow() = 0
    
	//y-sweep:--------------------------------------------------------------------------------------
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
	MultiFab::Copy(S_old, S_new, 0, 0, NUM_STATE, S_new.nGrow());
    FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    special_boundary(prob_domain, S_old, SimSettings, dx, V);
	
	//HALF TIME STEP RADIAL SOURCE TERM INTEGRATION:------------------------------------------------------
	if(Parameters.coord_sys == 1){
		for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
		{
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& statein 		 = 	 S_old[mfi];	//H
	    //--------------------------------------------------------------------H/		

		//H --------------------------------------------------------------------
		//access the array from the array box:
		
        Array4<Real> const& S_old_arr = statein.array();
        Array4<Real> const& S_new_arr = stateout.array();
		
		radial_sourceterm_integration(bx, S_old_arr, S_new_arr, dt/2.0, dxr, V);
		}

	}
	 	 
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
        
        amr_tagging(tagarr,bx,S_old_arr,dx,grad_max, grad_frac, prop, refinement_condition, V);
	    
	    
	}
    }//end local scope
    MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow()); //copy compute grad(rho) to S_new 
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
    
    //cout << "end of read_params() function" << endl;

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

