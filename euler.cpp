#include "structdefs.H"
#include "funcdefs.H"
#include "classdefs.H"
#include <cmath>
#include <algorithm>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_TagBox.H>
#include <AmrLevelAdv.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>


using namespace amrex;
using namespace std;

/* ---------------------GLOBAL VARIABLES: ----------------------------*/
//#define gamma 1.4   		//heat capacity ratio for standard air
//int nT; 					//final number of simulation time steps
string EoS;		//system Equation of State --> set in initialiseStructs function.
Plasma19 AirPlasma("mixture19_cns.txt");
string scheme;
/*--------------------------------------------------------------------*/

/*----------------------------VECTOR DEFINITIONS--------------------------------------------------
 * 		
 * 		names:			   State vector (S):	     Flux vector (F):    	 Flux vector (G):
 * 
 * 	{	mass	}			{	 rho	}			{	  rho*u		}		{	  rho*v		}
 *  {	mom_x	} 			{	rho*u	}			{   rho*u^2 +p	}		{    rho*v*u	}
 *  {	mom_y	}			{	rho*v	}			{	 rho*u*v	}		{   rho*v^2 +p	}
 *  {	 en		}			{	  E		}			{	u*[E + p]	}		{	v*[E + p]	}
 *  {	rho		}			{	 rho	}			{  	    0		}		{  	    0		}
 *  {	 u		}			{	  u		}			{  	    0		}		{  	    0		}
 *  {	 v		}			{	  v		}			{	    0		}		{  	    0		}
 *  {	 p 		}			{	  p 	}			{	    0 		}		{  	    0		}
 * 
 *   			 			  where: E = rho*e + 1/2*rho*|(u,v)|^2
 * 						  				and: e = EoS(rho,p)
 -----------------------------------------------------------------------------------------------*/

//the following variable name sare assumed for the current system:
Vector<string> assumed_variables{"mass", "mom_x", "mom_y", "en", "rho", "u", "v", "p", "EoS19", "EB", "level_set", "nx", "ny"};

void check_assumed_variables(Vector<string> variable_names){
	//check that the variable names which are called in this system 
	// all exist in the inputted variables list
	
	for(int i = 0; i<assumed_variables.size(); ++i){
		if(std::find(variable_names.begin(), variable_names.end(), assumed_variables[i]) == variable_names.end())
		{
			cout << assumed_variables[i] << endl;
			amrex::Abort("above assumed variable was not found in the inputted list of variable names");
		}else{
			//assumed variable found
		}
	}
	cout << "all assumed variables have been inputted correctly from amr.derive_plot_vars" << endl;
}


void initialiseStructs(SettingsStruct& SimSettings, ParameterStruct& Parameters)
{
	/* -----------------------------------------------------
     * Global simulation parameters.
     * -----------------------------------------------------*/
     
    ParmParse pp;
	pp.get("nsteps_max", SimSettings.nsteps_max);
	pp.get("plot_int", SimSettings.plot_int);
    pp.get("ndmin", SimSettings.ndmin);
    pp.get("startT", SimSettings.startT);
    pp.get("chatter", SimSettings.chatter);
    pp.get("testcase", SimSettings.testcase);
    pp.get("NCOMP", SimSettings.NCOMP);
    pp.get("NGROW", SimSettings.NGROW);
    pp.get("MUSCL", SimSettings.MUSCL);
    pp.get("scheme", SimSettings.scheme);
    pp.get("EoS", SimSettings.EoS);
    SimSettings.EBgeom = "none";
    pp.query("EBgeom", SimSettings.EBgeom);

    
    /* -----------------------------------------------------
     * Test Case specific parameters.
     * -----------------------------------------------------*/   
    
    ParmParse ppGeom("geometry");
    ppGeom.getarr("prob_lo", Parameters.prob_lo); 
    ppGeom.getarr("prob_hi", Parameters.prob_hi);
    ppGeom.get("coord_sys", Parameters.coord_sys);
    ppGeom.getarr("is_periodic", Parameters.is_periodic);
     
    ParmParse ppTest(SimSettings.testcase);
    ppTest.get("finalT", SimSettings.finalT);
    ppTest.get("n_vars", Parameters.n_vars);
    ppTest.getarr("vars", Parameters.vars);
    ppTest.getarr("PL", Parameters.PL);
    ppTest.getarr("PR", Parameters.PR);
    ppTest.getarr("n_cells", Parameters.n_cells);
    ppTest.get("IC", Parameters.IC);
    ppTest.getarr("bc_conds", Parameters.bc_conds);
    if(Parameters.IC == "x0" || Parameters.IC == "y0"){ 
		ppTest.get("x0", Parameters.x0);
	}else if(Parameters.IC == "mach_theta"){
		ppTest.get("x0", Parameters.x0);
		ppTest.get("mach", Parameters.mach);
		ppTest.get("theta",Parameters.theta);	
	}
	
	ParmParse ppEB(SimSettings.EBgeom);
	if(SimSettings.EBgeom == "sphere"){
		ppEB.get("R", Parameters.EB_R);
		ppEB.getarr("centre", Parameters.EB_centre);
	}else if(SimSettings.EBgeom == "plane"){
		ppEB.getarr("point", Parameters.EB_point);
		ppEB.getarr("normal", Parameters.EB_normal);
	}else if(SimSettings.EBgeom == "none"){
		//no additional settings
	}else{
		amrex::Abort("invalid EBgeom passed from inputs file");
	}
		
    
    ParmParse ppAmr("amr");
    ppAmr.getarr("derive_plot_vars", Parameters.problem_variables);
    ppAmr.get("max_level", Parameters.max_level);
    ppAmr.get("refinement_prop", Parameters.refinement_prop);
    ppAmr.get("refinement_condition", Parameters.refinement_condition);
    ppAmr.getarr("refinement_grad_fracs", Parameters.refinement_grad_fracs);
    
    //pp.getarr("problem_variables", Parameters.problem_variables);
      
    if( Parameters.problem_variables.size() != SimSettings.NCOMP){ 
		cout << "Parameters.problem_variables.size() : " << Parameters.problem_variables.size() << endl;
		cout << "SimSettings.NCOMP : " << SimSettings.NCOMP << endl;
		amrex::Abort("problem variables list must be of length NCOMP");
	}
	
	check_assumed_variables(Parameters.problem_variables);
    
    //set the global variable
    EoS = SimSettings.EoS;
    if(EoS != "ideal" && EoS != "plasma19"){
		amrex::Abort("invalid EoS input in test input file");
	}
    cout << "global variable EoS set to: " << EoS << endl;
    
    scheme = SimSettings.scheme;
    if(scheme != "SuperBee" && scheme != "MinBee" && scheme != "UltraBee" && scheme != "VanLeer"){
		amrex::Abort("invalid MUSCL limiter scheme in test input file");
	}
		
}

void setBoundaryConditions(Vector<BCRec>& bc_vec, const ParameterStruct& p, const int& nc, AccessVariable& V){
	
	/*----------------------------------------------
					  hi_bc[1]
				 ________________
				|				 |
				|				 |
	  lo_bc[0]  |				 |  hi_bc[0]
				|				 |
				|________________|
				
					  lo_bc[1]
	
	
	order of bc_conds = [left, right, bottom, top]
	
	-----------------------------------------------*/
	
	int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    
    //check BC's are all avild:
    for(int s = 0; s<4; ++s){
		if( !(p.bc_conds[s] == "transmissive") && !(p.bc_conds[s] == "reflective")){
			amrex::Abort("bc_conds from input file are invalid - must be 4 x transmissive/reflective");
		}
	}
    
    if(BL_SPACEDIM > 2){
		amrex::Abort("setBoundaryConditions not configured for BL_SPACEDIM > 2" );
	}
    
	for(int n=0; n<nc; ++n){
		for (int i = 0; i < BL_SPACEDIM; ++i) {
				
			if(i == 0){
				//lo_bc --> left:
				if(p.bc_conds[0] == "transmissive"){
					lo_bc[i] =  BCType::foextrap; //transmissive 
				}else if(p.bc_conds[0] == "reflective"){
					if(n == V["u"] || n == V["mom_x"]){
						lo_bc[i] = BCType::reflect_odd;  //reflect x velocities
					}else{
					lo_bc[i] = BCType::reflect_even; //transmit property
					}
				}
				//hi_bc --> right:
				if(p.bc_conds[1] == "transmissive"){
					hi_bc[i] = BCType::foextrap; //transmissive
				}else if(p.bc_conds[1] == "reflective"){
					if(n == V["u"] || n == V["mom_x"]){
						hi_bc[i] = BCType::reflect_odd;  //reflect x velocities
					}else{
						hi_bc[i] = BCType::reflect_even; //transmit property
					}
				}
			}
			if(i == 1){
				//lo_bc --> bottom:
				 if(p.bc_conds[2] == "transmissive"){
					lo_bc[i] =  BCType::foextrap; //transmissive
				}else if(p.bc_conds[2] == "reflective"){
					if(n == V["v"] || n == V["mom_y"]){
						lo_bc[i] = BCType::reflect_odd;  //reflect y velocities
					}else{
						lo_bc[i] = BCType::reflect_even; //transmit property
					}
				}
				//hi_bc --> top:
				if(p.bc_conds[3] == "transmissive"){
					hi_bc[i] =  BCType::foextrap; //transmissive
				}else if(p.bc_conds[3] == "reflective"){
					if(n == V["v"] || n == V["mom_y"]){
						hi_bc[i] = BCType::reflect_odd;  //reflect y velocities
					}else{
						hi_bc[i] = BCType::reflect_even; //transmit property
					}
				}
			}
			
		}//end i<dim 
		
		//store BC in bc_vec:
		BCRec bc(lo_bc, hi_bc);
		bc_vec[n] = bc;
		
	}//end n<ncomp		
	
}

//EoS function:
int EoS_properties(std::string eos, const double& rho, const double& p, const double& e, double& out_val, std::string prop){
	//EoS read from case: "ideal" or "EoS19"
	//properties: rho, p and e passed by reference
	//calculated property determined by strong prop: "e", "p", "a", "sigma_e"
	bool  err_flag = 0; 
	int EoS19;
	double gamma = 1.4;
	//double out_val;
	double rho_min_plasma19 = 1.001*1.0e-3; //plus 1% tolerance
	double p_min_plasma19 = 1.001*1013.25;
	double temp_p;
	
	//check if p or rho are below tabulated range:
	//if so, revert to ideal EoS:
	
	if(eos == "plasma19"){
		if(prop == "p"){
			//NB: this probably isn't rhobust
			temp_p = rho*(gamma - 1)*e;
			if(rho < rho_min_plasma19 || temp_p < p_min_plasma19){
				//switch to ideal
				eos = "ideal";
			}
		}else{
			if(rho < rho_min_plasma19 || p < p_min_plasma19){
				//switch to ideal
				eos = "ideal";
			}
		}
	}
	
	//cout << "1. (rho, p, e, prop) = " << rho << p << e << prop << endl;
	
	if(eos == "plasma19"){
		EoS19 = 1;
		if(prop == "e"){
			out_val = interp2D(AirPlasma.m_InternalEnergies, 
								AirPlasma.v_Pressures, AirPlasma.v_Densities, 
								rho, p, e, err_flag);
		}else if(prop == "a"){
			out_val = interp2D(AirPlasma.m_SoundSpeeds, 
								AirPlasma.v_Pressures, AirPlasma.v_Densities, 
								rho, p, e, err_flag);
		}else if(prop == "p"){
			out_val = interp2D_p(AirPlasma.m_InternalEnergies, 
								AirPlasma.v_Pressures, AirPlasma.v_Densities, 
								rho, p, e, err_flag);
		}else if(prop == "sigma_e"){
			out_val = interp2D(AirPlasma.m_SigmaElectrical, 
								AirPlasma.v_Pressures, AirPlasma.v_Densities, 
								rho, p, e, err_flag);
		}else if(prop == "T"){
			out_val = interp2D(AirPlasma.m_Temperatures, 
								AirPlasma.v_Pressures, AirPlasma.v_Densities, 
								rho, p, e, err_flag);
		}else if(prop == "gamma"){
			//requires new function
			out_val = gamma;
		}else{
			amrex::Abort("property function not defined");
		}
		if(err_flag){
			//catch remainaing out of range plasma properties
			//NOTE: this formulation may not be robust:
					//be sure to default to ideal for cases exceding low range and not high range
			//cout << "caught out of range plasma property" << endl;
			
			eos = "ideal"; //switch to ideal evaulation
			//recursive function call now using ideal EoS:
			//cout << "2. (rho, p, e, prop) = " << rho << p << e << prop << endl;
			//amrex::Abort("Plasma19 OOB: EoS_properties function abort");
			//EoS_properties("ideal", rho, p, e, prop);
		}
	}
	if(eos == "ideal"){
		EoS19 = 0;
		if(prop == "e"){
			out_val = p/(rho*(gamma-1));
		}else if(prop == "p"){
			out_val = rho*(gamma - 1)*e;
		}else if(prop == "a"){
			out_val = pow(gamma*p/rho,0.5);
		}else if(prop == "sigma_e"){
			//"can't produce conductivity vals from ideal EoS"
			amrex::Abort("can't produce conductivity vals from ideal EoS");
		}else if(prop == "T"){
			out_val = p/(rho*287);
		}else if(prop == "gamma"){
			out_val = gamma;
		}else{
			amrex::Abort("property function not defined");
		}
	}
	
	
	return EoS19;
}

//conversion functions:

void S_conWtoU(Box const& box, Array4<Real> const& prop_arr, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double e;
    
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				prop_arr(i,j,k,V["mass"]) = prop_arr(i,j,k,V["rho"]);
				prop_arr(i,j,k,V["mom_x"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["u"]);
				prop_arr(i,j,k,V["mom_y"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["v"]);
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), 0.0, e, "e"); //compute e
				prop_arr(i,j,k,V["en"]) = prop_arr(i,j,k,V["rho"])*e + 0.5*prop_arr(i,j,k,V["rho"])*\
											(pow(prop_arr(i,j,k,V["u"]),2) + pow(prop_arr(i,j,k,V["v"]),2));
			}
		}
	}
	
}

void S_conUtoW(Box const& box, Array4<Real> const& prop_arr, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double e, p;
    
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {	
				prop_arr(i,j,k,V["rho"]) = prop_arr(i,j,k,V["mass"]);
				prop_arr(i,j,k,V["u"])	 = prop_arr(i,j,k,V["mom_x"])/prop_arr(i,j,k,V["rho"]);
				prop_arr(i,j,k,V["v"])	 = prop_arr(i,j,k,V["mom_y"])/prop_arr(i,j,k,V["rho"]);
				e = 1.0/prop_arr(i,j,k,V["rho"])*(prop_arr(i,j,k,V["en"]) - \
						0.5*prop_arr(i,j,k,V["rho"])*(pow(prop_arr(i,j,k,V["u"]),2) + pow(prop_arr(i,j,k,V["v"]),2))); 
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), 0.0, e, p, "p"); //compute p
				prop_arr(i,j,k,V["p"]) = p;
			}
		}
	}
   
}	

void F_conWtoF(Box const& box, Array4<Real> const& prop_arr, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double e, E;
    
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {	
				prop_arr(i,j,k,V["mass"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["u"]);
				prop_arr(i,j,k,V["mom_x"]) = prop_arr(i,j,k,V["rho"])*pow(prop_arr(i,j,k,V["u"]),2)+prop_arr(i,j,k,V["p"]);
				prop_arr(i,j,k,V["mom_y"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["u"])*prop_arr(i,j,k,V["v"]);
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), 0.0, e, "e");
				E = prop_arr(i,j,k,V["rho"])*e + 0.5*prop_arr(i,j,k,V["rho"])*\
				    (pow(prop_arr(i,j,k,V["u"]),2) + pow(prop_arr(i,j,k,V["v"]),2));
				prop_arr(i,j,k,V["en"]) = prop_arr(i,j,k,V["u"])*(E+prop_arr(i,j,k,V["p"]));
			}
		}
	}
   
}			
	

//Conversion functions between classes:

ConsF conWtoF(Prim const& W){
    double mass = W.rho*W.u_vec[0];
    Vector<Real> mom_vec(3);
    mom_vec[0] = W.rho*pow(W.u_vec[0],2)+W.p;
    mom_vec[1] = W.rho*W.u_vec[0]*W.u_vec[1];
    mom_vec[2] = 0;
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E =  W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    double en = W.u_vec[0]*(E + W.p );
    ConsF F(mass, mom_vec, en);
    return F;
}

Prim conUtoW(ConsU const& U){
    double rho = U.rho;
    Vector<Real> u_vec(3);
    u_vec[0] = U.rhou_vec[0]/rho;
    u_vec[1] = U.rhou_vec[1]/rho;
    u_vec[2] = 0;
    double e = 1/rho*(U.E - 0.5*rho*pow(mag(u_vec),2));
    double p;
    int EoS19 = EoS_properties(EoS, rho, 0.0, e, p, "p");
    Prim W(rho, u_vec, p);
    W.r = U.r;
    return W;
}

ConsU conWtoU(Prim const& W){
    double rho = W.rho;
    Vector<Real> rhou_vec(3);
    rhou_vec[0] = W.rho*W.u_vec[0];
    rhou_vec[1] = W.rho*W.u_vec[1];
    rhou_vec[2] = W.rho*W.u_vec[2];
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E = W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    ConsU U(rho, rhou_vec, E);
    U.r = W.r;
    return U;
}


ConsF conWtoG(Prim const& W){
    double mass = W.rho*W.u_vec[1];
    Vector<Real> mom_vec(3);
    mom_vec[0] = W.rho*W.u_vec[1]*W.u_vec[0];
    mom_vec[1] = W.rho*W.u_vec[1]*W.u_vec[1]+W.p;
    mom_vec[2] = 0;
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E =  W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    double en = W.u_vec[1]*(E + W.p);
    ConsF G(mass, mom_vec, en);
    return G;
}

ConsF conWtoSE(Prim const& W){
    double mass = ((W.r > 0)? -1/W.r * W.rho*W.u_vec[0] : 0);
    Vector<Real> mom_vec(3);
    mom_vec[0] = ((W.r > 0)? -1/W.r * W.rho*pow(W.u_vec[0],2) : 0);
    mom_vec[1] = 0;
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E =  W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    double en = ((W.r > 0)? -1/W.r * W.u_vec[0]*(E+W.p) : 0);
    ConsF SE(mass, mom_vec, en);
    return SE;
}

//rankine-hugoniot calculation:
Prim rankine_hugoniot(Prim const& W, double const& MS, Vector<Real> const& normal_vec){
	//rankine hugoniot relation derived for ideal gas with gamma = 1.4
	//assume intitial conditions are in ideal regime when this function 
	//is called
	//also assumes a normal shock along normalised_vec:
	double gamma = 1.4;
	double eR;
	int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, eR, "e");
	if(EoS19){ 
		cout << "WARNING: rankine hugoniot calc assumes ideal gas regime \
					but plasma validity range redetected" << endl;
		EoS_properties("ideal", W.rho, W.p, 0.0, eR, "e");
	}
	double aR;
	EoS19 = EoS_properties("ideal", W.rho, W.p, eR, aR, "a");
	Vector<Real> MR_vec{W.u_vec[0]/aR, W.u_vec[1]/aR, 0};
	double MR 	= mag(MR_vec);
	double rho 	= (((gamma+1)*pow(MR-MS,2))/((gamma-1)*pow(MR-MS,2)+2))*W.rho;
	double p   	= ((2*gamma*pow(MR-MS,2)-(gamma-1))/(gamma+1))*W.p;
	double S3 	= MS*aR;
	Vector<Real> S3_vec = S3*normal_vec;
	Vector<Real> u_vec = (1-W.rho/rho)*S3_vec + (W.rho/rho)*W.u_vec;
	
	cout << "computed properties from analytic rankine hugoniot relations: \n";
	cout << "gamma = " << gamma << endl;
	cout << "aR = " << aR << endl;
	cout << "MR = " << MR << endl;
	cout << "rho = " << rho << endl;
	cout << "p = " << p << endl;
	cout << "S3_vec : \n";
	printVector(S3_vec);
	cout << "u_vec : \n";
	printVector(u_vec);

	Prim WL(rho, u_vec, p);
    return WL;
	
}

//problem initialisation:
void initial_conditions(int const& nc, Box const& box, Array4<Real> const& prop_arr, const Real* dx, ParameterStruct const& p, AccessVariable& V){

    const auto lo = lbound(box);
    const auto hi = ubound(box);
	
	//amrex::Print() << "lo : " << lo << ", hi : " << hi << "\n";

	//will need changing to *dx when cell sizes change in AMR
    int int_x0 = (p.x0/p.prob_hi[0])*p.n_cells[0];
    int int_y0 = (p.x0/p.prob_hi[1])*p.n_cells[1];
	
	//AccessVariable V(p.problem_variables);
	
	//cout << "PL and PR size() " << p.PL.size() << " : " << p.PR.size() << endl;
	
	Real R = 0.2;
	Real h;
	
	//set all intitial data to zero - this seems to make visit happy:
	for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				for(int n=0; n<nc; ++n){
					prop_arr(i,j,k,n) = 0.0;
				}
			}
		}
	}
	
	
	if(p.IC == "source"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					
					h = sqrt(pow((i+0.5)*dx[0]-0.5,2) + pow((j+0.5)*dx[1]-0.5,2));
					if(h <= R){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n];
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}
					
				}
			}
		}
	}		
	}else if(p.IC == "x0"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(i*dx[0] <= p.x0){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n]; 
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}
					
				}
			}
		}
		}
	}else if(p.IC == "y0"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(j*dx[1] <= p.x0){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n]; 
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}
					
				}
			}
		}
		}
	}else if(p.IC == "mach_theta"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(i*dx[0] <= p.x0 + (j*dx[1])*tan(p.theta)){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n]; 
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}	
				}
			}
		}
		}
	}else{
		amrex::Abort("invalid intitial condition - check inputs file for adv.IC");
	}
	
	S_conWtoU(box, prop_arr, V);
	
}




void special_boundary(Box const& dom, MultiFab& S_old, SettingsStruct const& sim, const Real* dx, AccessVariable& V){
    
    //S_new box passed here. box is whole physical domain as a box
    const auto prob_lo = lbound(dom); 
    const auto prob_hi = ubound(dom);
    
    //cout << "prob_hi.x = " << prob_hi.x << endl;
    
    ParmParse ppTest(sim.testcase);
    string specialBC = "none";
    ppTest.query("specialBC", specialBC);
    if(specialBC == "partial_right" || specialBC == "partial_lower"){
		Real dist = 0.2;
		ppTest.query("transmissive_dist", dist);
		
		for (MFIter mfi(S_old, true); mfi.isValid(); ++mfi){
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    Array4<Real> const& S_old_arr = S_old[mfi].array();	//H
	    //--------------------------------------------------------------------H/		
		
		const auto lo = lbound(bx); //lo bound of real domain (excluding ghost cells)
		const auto hi = ubound(bx);
		
		if(specialBC == "partial_lower" && lo.y == prob_lo.y && lo.x*dx[0] <= dist){	
		//overwrite bottom reflective portion with transmissive
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = prob_lo.y-S_old.nGrow(); j < prob_lo.y; ++j) {
					for(int i = lo.x-S_old.nGrow(); i*dx[0] <= dist; ++i) {
						S_old_arr(i,j,k,V["mass"]) 	= S_old_arr(i,prob_lo.y,k,V["mass"]);
						S_old_arr(i,j,k,V["mom"]) 	= S_old_arr(i,prob_lo.y,k,V["mom"]);
						S_old_arr(i,j,k,V["en"]) 	= S_old_arr(i,prob_lo.y,k,V["en"]);
						S_old_arr(i,j,k,V["rho"]) 	= S_old_arr(i,prob_lo.y,k,V["rho"]);
						S_old_arr(i,j,k,V["u"]) 	= S_old_arr(i,prob_lo.y,k,V["u"]);
						S_old_arr(i,j,k,V["v"]) 	= S_old_arr(i,prob_lo.y,k,V["v"]);
						S_old_arr(i,j,k,V["p"]) 	= S_old_arr(i,prob_lo.y,k,V["p"]);
					}
				}
			}//end i,j,k loop
		
		}//end partial_lower altered portion update
		
		if(specialBC == "partial_right" && hi.x == prob_hi.x && lo.y*dx[1] <= dist){
		//overwrite bottom reflective portion with transmissive
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = lo.y; j*dx[1] <= dist; ++j) {
					for(int i = prob_hi.x; i < prob_hi.x + S_old.nGrow(); ++i) {
						S_old_arr(i,j,k,V["mass"]) 	= S_old_arr(prob_hi.x,j,k,V["mass"]);
						S_old_arr(i,j,k,V["mom"]) 	= S_old_arr(prob_hi.x,j,k,V["mom"]);
						S_old_arr(i,j,k,V["en"]) 	= S_old_arr(prob_hi.x,j,k,V["en"]);
						S_old_arr(i,j,k,V["rho"]) 	= S_old_arr(prob_hi.x,j,k,V["rho"]);
						S_old_arr(i,j,k,V["u"]) 	= S_old_arr(prob_hi.x,j,k,V["u"]);
						S_old_arr(i,j,k,V["v"]) 	= S_old_arr(prob_hi.x,j,k,V["v"]);
						S_old_arr(i,j,k,V["p"]) 	= S_old_arr(prob_hi.x,j,k,V["p"]);
					}
				}
			}//end i,j,k loop
		
		}//end partial_right altered portion update
		
		}//end MFIter
	}else{
		//no additional boundary treatements necessary
	}
}
 

double r_calc(double n_val, double l_val){
    //alternative conventions for dividing by zero RHS slope:
    // -3 encoded to represent -ve infinity
    // +3 encoded to represent +ve infinity 
    double r = l_val;
    if(n_val == 0){
        if(l_val== 0){
            r = 0;           //both left and right slopes are zero
        }else if(l_val < 0){
            r = -3;          //negative slope on left, zero slope on right
        }else{
            r = 3;           //positive slope on left, zero slope on right
        }
    }else{ r = r/n_val;}     //permissible division by delta_i_n

    return r;
}

double delta_bar_calc(string& scheme, double& r, double& delta_i, double& sigma_R, double& sigma_L){
    double delta_bar;
    if(scheme == "SuperBee"){
        if(r < 0){
            delta_bar = 0;
        }else if(0 <= r && r < 0.5){
            delta_bar = 2*r*delta_i;
        }else if(0.5 <= r && r < 1){
            delta_bar = delta_i;
        }else if(1 <= r && r < 2){
            delta_bar = ((r < sigma_R)? 
                    r*delta_i : sigma_R*delta_i);
        }else{
            delta_bar = ((sigma_R < 2)? 
                    sigma_R*delta_i : 2*delta_i);
        }
    }else if(scheme == "MinBee"){

        if(r <= 0){
            delta_bar = 0;
        }else if(0 <= r && r <= 1.0){
            delta_bar = r*delta_i;
        }else{
            delta_bar = ((1 < sigma_R)? delta_i : sigma_R*delta_i);
        }
    }else if(scheme == "VanLeer"){
        if(r < 0){
            delta_bar = 0;
        }else{
            delta_bar = ((2*r/(1+r) < sigma_R)? 
                    (2*r/(1+r))*delta_i : sigma_R*delta_i);
        }
    }else if(scheme == "UltraBee"){
        if(r <= 0){
            delta_bar = 0;
        }else{
            delta_bar = ((sigma_L < sigma_R)? 
                    sigma_L*delta_i : sigma_R*delta_i);
        }
    }else{
		cout << "no scheme name given to delta_bar_calc" << endl;
        return delta_i;
    }

    return delta_bar;
}

Prim delta(double& omega, Prim& Wl, Prim& Wi, Prim& Wr, string& scheme){
    //omega E [-1, 1]
    //Wl = W_{i-1}, Wi = Wi, Wr = W_{i+1}
    Prim delta_i_l(Wi.rho - Wl.rho, Wi.u_vec - Wi.u_vec, 
                    Wi.p - Wl.p);
    Prim delta_i_n(Wr.rho - Wi.rho, Wr.u_vec - Wi.u_vec, 
                    Wr.p - Wi.p);
                    
    //r calc for each W property
    double r_rho = r_calc(delta_i_n.rho, delta_i_l.rho);
    vector<double> r_u_vec = delta_i_l.u_vec;
    for(int i = 0; i<3; ++i){
		r_u_vec[i] = r_calc(delta_i_n.u_vec[i], delta_i_l.u_vec[i]);
	}
    double r_p = r_calc(delta_i_n.p, delta_i_l.p);
	
	//delta calc:
    Prim delta_i(0.5*(1+omega)*delta_i_l.rho+0.5*(1-omega)*delta_i_n.rho,
                  0.5*(1+omega)*delta_i_l.u_vec+0.5*(1-omega)*delta_i_n.u_vec,
                  0.5*(1+omega)*delta_i_l.p+0.5*(1-omega)*delta_i_n.p);
	
	//sigma calcs:
	Vector<Real> sigma_u_vec(3);
	
	for(int i=0; i<3; ++i){
		sigma_u_vec[i] = 2.0/(1-omega+(1+omega)*r_u_vec[i]);
	}
	
    Prim sigma_R(2/(1-omega+(1+omega)*r_rho), sigma_u_vec, 
                    2/(1-omega+(1+omega)*r_p));
                    
    for(int i=0; i<3; ++i){
		sigma_u_vec[i] = 2.0*r_u_vec[i]/(1-omega+(1+omega)*r_u_vec[i]);
	}

    Prim sigma_L(2*r_rho/(1-omega+(1+omega)*r_rho), sigma_u_vec, 
                 2*r_p/(1-omega+(1+omega)*r_p));

    //adjustments for infinity encoded cases:
    if(r_rho == -3 || r_rho == 3){
        sigma_R.rho = 0;
    }
    for(int i=0; i<3; ++i){
		if(r_u_vec[i] == -3 || r_u_vec[i] == 3){
			sigma_R.u_vec[i] = 0;
		}
	}
    if(r_p == -3 || r_p == 3){
        sigma_R.p = 0;
    }
    
    //delta bar calc for each W property
    Vector<Real> v0(3); //zero vec for initialisation 
    Prim delta_bar(0,v0,0);
    
    delta_bar.rho = delta_bar_calc(scheme, r_rho, delta_i.rho, sigma_R.rho, sigma_L.rho);
    for(int i=0; i<3; ++i){
		delta_bar.u_vec[i] = delta_bar_calc(scheme, r_u_vec[i], delta_i.u_vec[i], sigma_R.u_vec[i], sigma_L.u_vec[i]);
	}
	delta_bar.p = delta_bar_calc(scheme, r_p, delta_i.p, sigma_R.p, sigma_L.p);
    
    return delta_bar;

}


void calc_fluxes_HLLC_x(Box const& box, Array4<Real> const& S_old, Array4<Real> const& F, 
					Real const dt, Real const dx, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  

	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF FL(0,v0,0), FR(0,v0,0);
	ConsF FLstar(0,v0,0), FRstar(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Fi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
    
    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x+1; ++i) { 
				LL = i-2;
				L = i-1;
				I = i;
				R = i+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(LL,j,k,V["u"]);
				v0[1] = S_old(LL,j,k,V["v"]);
				Wll = Prim(S_old(LL,j,k,V["rho"]), v0, S_old(LL,j,k,V["p"]));
				v0[0] = S_old(R,j,k,V["u"]);
				v0[1] = S_old(R,j,k,V["v"]);
				Wr = Prim(S_old(R,j,k,V["rho"]), v0, S_old(R,j,k,V["p"]));
				}
				
				v0[0] = S_old(L,j,k,V["u"]);
				v0[1] = S_old(L,j,k,V["v"]);
				Wl = Prim(S_old(L,j,k,V["rho"]), v0, S_old(L,j,k,V["p"]));
				v0[0] = S_old(I,j,k,V["u"]);
				v0[1] = S_old(I,j,k,V["v"]);
				Wi = Prim(S_old(I,j,k,V["rho"]), v0, S_old(I,j,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.rho - dt/dx*Wl.rho*delta_l.u_vec[0]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[0] - dt/dx*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[1]);
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[2]);
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.p - dt/dx*(Wl.rho*pow(aL,2))*delta_l.u_vec[0]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.rho + dt/dx*Wi.rho*delta_i.u_vec[0]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[0] + dt/dx*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[1]);
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[2]);
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.p + dt/dx*(Wi.rho*pow(aL,2))*delta_i.u_vec[0]);
				
				WL = WbarL; WR = WbarR;
				}else{
					WL = Wl; WR = Wi;
				}
				//--------------------------------------MUSCL half time step----------------------------------------------
				
				//direct wave speed estimates:
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
				SL = ((WL.u_vec[0] - aL < WR.u_vec[0] - aR)? WL.u_vec[0] - aL : WR.u_vec[0] - aR);
				SR = ((WL.u_vec[0] + aL > WR.u_vec[0] + aR)? WL.u_vec[0] + aL : WR.u_vec[0] + aR);
				Sstar = (WR.u_vec[0]*WR.rho*(SR-WR.u_vec[0]) - WL.u_vec[0]*WL.rho*(SL-WL.u_vec[0]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[0])-WL.rho*(SL-WL.u_vec[0]));
				Splus = ((fabs(WL.u_vec[0])+aL > fabs(WR.u_vec[0])+aR)? \
						  fabs(WL.u_vec[0])+aL : fabs(WR.u_vec[0])+aR);
						  
						  
				//W-->U and W-->F conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				FL = conWtoF(WL);
				FR = conWtoF(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Fi = FL;
				}else if(SL <= 0 && 0 < Sstar){
					//cout << "left star state \n";
					//left star state terms:
					ULstar.rho = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar));
					ULstar.rhou_vec[0] = WL.rho*Sstar*((SL-WL.u_vec[0])/(SL-Sstar));
					ULstar.rhou_vec[1] = WL.rho*WL.u_vec[1]*((SL-WL.u_vec[0])/(SL-Sstar)); 
					ULstar.E = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar))*(UL.E/WL.rho+(Sstar-WL.u_vec[0])\
								*(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[0]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						amrex::Abort("aborted in MUSCL_HLLC_x");
						break;
					}
					//left star state flux:
					FLstar.mass = FL.mass + SL*(ULstar.rho - UL.rho);
					FLstar.mom_vec = FL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
					FLstar.en = FL.en + SL*(ULstar.E - UL.E);
					Fi = FLstar;
				}else if(Sstar <= 0 && 0 < SR){
					//cout << "right star state \n";
					//right star state terms:
					URstar.rho = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar));
					URstar.rhou_vec[0] = WR.rho*Sstar*((SR-WR.u_vec[0])/(SR-Sstar));
					URstar.rhou_vec[1] = WR.rho*WR.u_vec[1]*((SR-WR.u_vec[0])/(SR-Sstar)); 
					URstar.E = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar))*(UR.E/WR.rho+(Sstar-WR.u_vec[0])\
								*(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[0]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						amrex::Abort("aborted in MUSCL_HLLC_x");
						break;
					}
					FRstar.mass = FR.mass + SR*(URstar.rho - UR.rho);
					FRstar.mom_vec = FR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
					FRstar.en = FR.en + SR*(URstar.E - UR.E);
					Fi = FRstar;
				}else if(SR <= 0){
					Fi = FR; 
				}
				
				//-------------------convert F to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					F(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite calculated fluxes:
				F(i,j,k,V["mass"]) = Fi.mass;
				F(i,j,k,V["mom_x"]) = Fi.mom_vec[0];
				F(i,j,k,V["mom_y"]) = Fi.mom_vec[1];
				F(i,j,k,V["en"]) = Fi.en;
				//----------------------------------------------------------------										
			
			}
		}
	}//end i,j,k
	
}


void calc_fluxes_HLLC_y(Box const& box, Array4<Real> const& S_old, Array4<Real> const& G, 
					Real const dt, Real const dy, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  
    
  
	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF GL(0,v0,0), GR(0,v0,0);
	ConsF GLstar(0,v0,0), GRstar(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Gi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
    
    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y+1; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) { 
				LL = j-2;
				L = j-1;
				I = j;
				R = j+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(i,LL,k,V["u"]);
				v0[1] = S_old(i,LL,k,V["v"]);
				Wll = Prim(S_old(i,LL,k,V["rho"]), v0, S_old(i,LL,k,V["p"]));
				v0[0] = S_old(i,R,k,V["u"]);
				v0[1] = S_old(i,R,k,V["v"]);
				Wr = Prim(S_old(i,R,k,V["rho"]), v0, S_old(i,R,k,V["p"]));
				}
				
				v0[0] = S_old(i,L,k,V["u"]);
				v0[1] = S_old(i,L,k,V["v"]);
				Wl = Prim(S_old(i,L,k,V["rho"]), v0, S_old(i,L,k,V["p"]));
				v0[0] = S_old(i,I,k,V["u"]);
				v0[1] = S_old(i,I,k,V["v"]);
				Wi = Prim(S_old(i,I,k,V["rho"]), v0, S_old(i,I,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.rho - dt/dy*Wl.rho*delta_l.u_vec[1]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[0]);
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[1] - dt/dy*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[2]);
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.p - dt/dy*(Wl.rho*pow(aL,2))*delta_l.u_vec[1]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.rho + dt/dy*Wi.rho*delta_i.u_vec[1]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[0]);
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[1] + dt/dy*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[2]);
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.p + dt/dy*(Wi.rho*pow(aL,2))*delta_i.u_vec[1]);
				
				WL = WbarL; WR = WbarR;
				}else{
				WL = Wl; WR = Wi;
				}
				//--------------------------------------MUSCL half time step----------------------------------------------
				//direct wave speed estimates:
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
				SL = ((WL.u_vec[1] - aL < WR.u_vec[1] - aR)? WL.u_vec[1] - aL : WR.u_vec[1] - aR);
				SR = ((WL.u_vec[1] + aL > WR.u_vec[1] + aR)? WL.u_vec[1] + aL : WR.u_vec[1] + aR);
				Sstar = (WR.u_vec[1]*WR.rho*(SR-WR.u_vec[1]) - WL.u_vec[1]*WL.rho*(SL-WL.u_vec[1]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[1])-WL.rho*(SL-WL.u_vec[1]));
				Splus = ((fabs(WL.u_vec[1])+aL > fabs(WR.u_vec[1])+aR)? \
						  fabs(WL.u_vec[1])+aL : fabs(WR.u_vec[1])+aR);
						  
						  
				//W-->U and W-->G conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				GL = conWtoG(WL);
				GR = conWtoG(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Gi = GL;
				}else if(SL <= 0 && 0 < Sstar){
					//cout << "left star state \n";
					//left star state terms:
					ULstar.rho = WL.rho*((SL-WL.u_vec[1])/(SL-Sstar));
					ULstar.rhou_vec[0] = WL.rho*WL.u_vec[0]*((SL-WL.u_vec[1])/(SL-Sstar));
					ULstar.rhou_vec[1] = ULstar.rho*Sstar; 
					ULstar.E = ((SL-WL.u_vec[1])/(SL-Sstar))*(UL.E+WL.rho*(Sstar-WL.u_vec[1])*\
									(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[1]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						break;
					}
					//left star state flux:
					GLstar.mass = GL.mass + SL*(ULstar.rho - UL.rho);
					GLstar.mom_vec = GL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
					GLstar.en = GL.en + SL*(ULstar.E - UL.E);
					Gi = GLstar;
				}else if(Sstar <= 0 && 0 < SR){
					//cout << "right star state \n";
					//right star state terms:
					URstar.rho = WR.rho*((SR-WR.u_vec[1])/(SR-Sstar));
					URstar.rhou_vec[0] = WR.rho*WR.u_vec[0]*((SR-WR.u_vec[1])/(SR-Sstar));
					URstar.rhou_vec[1] = URstar.rho*Sstar; 
					URstar.E = ((SR-WR.u_vec[1])/(SR-Sstar))*(UR.E + WR.rho*(Sstar-WR.u_vec[1])*\
									(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[1]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						break;
					}
					GRstar.mass = GR.mass + SR*(URstar.rho - UR.rho);
					GRstar.mom_vec = GR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
					GRstar.en = GR.en + SR*(URstar.E - UR.E);
					Gi = GRstar;
				}else if(SR <= 0){
					Gi = GR; 
				}
				
				//-------------------convert G to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					G(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite computed fluxes:
				G(i,j,k,V["mass"]) = Gi.mass;
				G(i,j,k,V["mom_x"]) = Gi.mom_vec[0];
				G(i,j,k,V["mom_y"]) = Gi.mom_vec[1];
				G(i,j,k,V["en"]) = Gi.en;
				//----------------------------------------------------------------		
				
				
			}
		}
	}//end i,j,k
	
}//end HLLC_y
			
	
void cons_update(const std::string sweep, Box const& box, Array4<Real> const& F, Array4<Real> const& S_old, 
        Array4<Real> const& S_new, Real const dt, Real const dx, int const nc, AccessVariable& V){
    
    const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);
    Prim Wn(0,v0,0);
    Prim Wnew(0,v0,0);
    ConsU Un(0,v0,0);
    ConsU Unew(0,v0,0);
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
				
				//-------------------convert P to U---------------------------
				v0[0] = S_old(i,j,k,V["u"]);
				v0[1] = S_old(i,j,k,V["v"]);
				Wn = Prim(S_old(i,j,k,V["rho"]), v0, S_old(i,j,k,V["p"]));
				Un = conWtoU(Wn);
				
				//-------------conservative update on Un----------------------
				//where flux_prop_arr_x(rho*u, rho*u^2+p, rho*u*v, u(E+p)
                if(sweep == "x-sweep"){
					Unew.rho = Un.rho + dt/dx * (F(i,j,k,V["mass"]) - F(i+1,j,k,V["mass"])); 
					Unew.rhou_vec[0] = Un.rhou_vec[0] + dt/dx * (F(i,j,k,V["mom_x"]) - F(i+1,j,k,V["mom_x"])); 
					Unew.rhou_vec[1] = Un.rhou_vec[1] + dt/dx * (F(i,j,k,V["mom_y"]) - F(i+1,j,k,V["mom_y"])); 
					Unew.E = Un.E + dt/dx * (F(i,j,k,V["en"]) - F(i+1,j,k,V["en"])); 
                    
                }else if(sweep == "y-sweep"){
                    Unew.rho = Un.rho + dt/dx * (F(i,j,k,V["mass"]) - F(i,j+1,k,V["mass"])); 
					Unew.rhou_vec[0] = Un.rhou_vec[0] + dt/dx * (F(i,j,k,V["mom_x"]) - F(i,j+1,k,V["mom_x"])); 
					Unew.rhou_vec[1] = Un.rhou_vec[1] + dt/dx * (F(i,j,k,V["mom_y"]) - F(i,j+1,k,V["mom_y"])); 
					Unew.E = Un.E + dt/dx * (F(i,j,k,V["en"]) - F(i,j+1,k,V["en"])); 
                }
				
				//--------------convert Unew to S_new-----------------------
				//Wnew = conUtoW(Unew);
				S_new(i,j,k,V["mass"]) = Unew.rho;
				S_new(i,j,k,V["mom_x"]) = Unew.rhou_vec[0];
				S_new(i,j,k,V["mom_y"]) = Unew.rhou_vec[1];
				S_new(i,j,k,V["en"]) = Unew.E;
								
            }
        }
    }//end i,j,k loop
    
    S_conUtoW(box, S_new, V);
    
}

void radial_sourceterm_integration(Box const& box, Array4<Real> const& S_old, Array4<Real> const& S_new, 
										Real const dt, Real const dx, AccessVariable& V){
    
    //NOTE for half time step integration: assumes dt is passed as: dt/2
    const auto lo = lbound(box);
    const auto hi = ubound(box);
    
	Real R;
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
				
				//-------------------radial position---------------------------
				R = (i+0.5)*dx; //cell centre: should avoid div by R=0
				//-------------radial source term integration----------------------
				S_new(i,j,k,V["mass"]) = S_old(i,j,k,V["mass"]) + dt * (-1.0/R * S_old(i,j,k,V["rho"]) * S_old(i,j,k,V["u"]));
				S_new(i,j,k,V["mom_x"]) = S_old(i,j,k,V["mom_x"]) + dt * (-1.0/R * S_old(i,j,k,V["rho"]) * pow(S_old(i,j,k,V["u"]), 2));
				S_new(i,j,k,V["en"]) = S_old(i,j,k,V["en"]) + dt * (-1.0/R * S_old(i,j,k,V["u"]) * (S_old(i,j,k,V["en"]) + S_old(i,j,k,V["p"])) );								
            }
        }
    }//end i,j,k loop
    
    //update primitive variables from new conservative variables: 
    S_conUtoW(box, S_new, V);
    
}

Real max_dim_wave_speed(Box const& bx, Array4<Real> const& S_new, int const i, AccessVariable& V){
	
	//convert from i dimension to n_index
	
	//this calculation should be sufficient to track EoS19 variable (?)
	
	string u_dim;
	if(i == 0){
		u_dim = "u";
	}else if(i==1){
		u_dim = "v";
	}
	
	
	Real max_wave_speed = 0.0;
	Real sound_speed, wave_speed;
	
	const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
				
				S_new(i,j,k,V["EoS19"]) = fabs(EoS_properties(EoS, S_new(i,j,k,V["rho"]), S_new(i,j,k,V["p"]), 0.0, sound_speed, "a")); 
				wave_speed = sound_speed + fabs(S_new(i,j,k,V[u_dim]));
				
				if(wave_speed > max_wave_speed){
					max_wave_speed = wave_speed;
				}
				
			}
		}
	}//end i,j,k loop
	
	//cout << "max_wave_speed = " << max_wave_speed << endl;
	return max_wave_speed;
	
}



void rescale_fluxes(const std::string sweep, Box const& box, Array4<Real> const& F, 
					Real const dt, Real const dx, int const nc, AccessVariable& V){
	//nodal box is passed in 
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);

	int ynode;
	int xnode;
	if(sweep == "x-sweep"){
		ynode = 0;
		xnode = 1;
	}else if(sweep == "y-sweep"){
		ynode = 1;
		xnode = 0;
	}
	
	//int const nc = 1;
	for(int n = 0; n<nc; ++n){
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y + ynode; ++j) {
            for(int i = lo.x; i <= hi.x + xnode; ++i) {
				
                if(sweep == "x-sweep"){
                    F(i,j,k,n) = F(i,j,k,n) * (dt * dx); //dx = dy for Fx
                }else if(sweep == "y-sweep"){
                    F(i,j,k,n) = F(i,j,k,n) * (dt * dx); //dx = dx for Fy
                }
				

            }
        }
    }
	}
    
}

void determine_grad_max(const Box& box, Array4<Real> const& S_old_arr, const Real* dx, string prop, 
							string const refinement_condition, Real& grad_max, AccessVariable& V)
{
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real local_grad, phi_x, phi_y;
    
    Real alpha = 1.0; 	//log scaling parameter
    
    //find max gradient
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				phi_x = fabs((S_old_arr(i+1,j,k,V[prop])-S_old_arr(i-1,j,k,V[prop]))/(2*dx[0]));
				phi_y = fabs((S_old_arr(i,j+1,k,V[prop])-S_old_arr(i,j-1,k,V[prop]))/(2*dx[1]));
				if(refinement_condition == "grad"){
					local_grad = sqrt(pow(phi_x,2) + pow(phi_y,2));
					//local_grad = (phi_x > phi_y ? phi_x : phi_y);
					if(prop == "mass"){
						local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
					}															//for rho only as must be non-zero everywhere
				}else if(refinement_condition == "log_grad"){
					local_grad = (phi_x > phi_y ? phi_x : phi_y); 
					if(prop == "mass"){
						local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
					}															//for rho only as must be non-zero everywhere
					local_grad = log(alpha*local_grad + 1);
				}else{
					amrex::Abort("refinement_condition invalid - check inputs file");
				}
				if(local_grad > grad_max){
					grad_max = local_grad;
				}
			}
		}
	}//end i,j,k
	
	//cout << "grad_max = " << grad_max << endl;
	
}
	

void amr_tagging(Array4<char> const& tagarr, const Box& box, Array4<Real> const& S_old_arr, const Real* dx, Real grad_max, Real grad_frac,
					string prop, string const refinement_condition, AccessVariable& V)
{
	//tagging function for refinement
	//based on gradient of rho
	//assumes rho is ncomp element 0
	//box has S_new dimensions, but S_old_arr has at least +1 ngrow
	//therefore, we can use central different function for gradient
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real local_grad, phi_x, phi_y;
    //Real prop_max = 0.5;
        
    Real alpha = 1.0; 	//log scaling parameter
    
	//set tags based on gradient relative to max gradient: 
	//ie. prop_max is proportion of max gradient
	for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				phi_x = fabs((S_old_arr(i+1,j,k,V[prop])-S_old_arr(i-1,j,k,V[prop]))/(2*dx[0]));
				phi_y = fabs((S_old_arr(i,j+1,k,V[prop])-S_old_arr(i,j-1,k,V[prop]))/(2*dx[1]));
				if(refinement_condition == "grad"){
					local_grad = sqrt(pow(phi_x,2) + pow(phi_y,2));
					//local_grad = (phi_x > phi_y ? phi_x : phi_y);
					if(prop == "mass"){
						local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
					}															//for rho only as must be non-zero everywhere
				}else if(refinement_condition == "log_grad"){
					local_grad = sqrt(pow(phi_x,2) + pow(phi_y,2));
					//local_grad = (phi_x > phi_y ? phi_x : phi_y); 
					if(prop == "mass"){
						local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
					}															//for rho only as must be non-zero everywhere
					local_grad = log(alpha*local_grad + 1);
				}else{
					amrex::Abort("refinement_condition invalid - check inputs file");
				}
				/*
				if(prop == "mass"){
					S_old_arr(i,j,k,V["grad_rho"]) = local_grad;
				}
				*/
				if(local_grad > grad_frac*grad_max){
					tagarr(i,j,k,V[prop]) = TagBox::SET;
					//cout << "tagged" << endl;
				}else{
					tagarr(i,j,k,V[prop]) = TagBox::CLEAR;
				}
			}
		}
	}//end i,j,k
		
}


void calc_fluxes_EB_HLLC_x(Box const& box, Array4<Real> const& S_old, Array4<Real> const& F, 
					Real const dt, Real const dx, int const nc, AccessVariable& V, bool MUSCL, Array4<const EBCellFlag>& flags_arr)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
     Vector<Real> v0(3);  

	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF FL(0,v0,0), FR(0,v0,0);
	ConsF FLstar(0,v0,0), FRstar(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Fi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
    
    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x+1; ++i) { 
				
				for(int n = 0; n<nc; ++n){
						F(i,j,k,n) = 0; //default set all fluxes to zero
					}
				
				if(flags_arr(i,j,k).isCovered()){
					//fluxes remain zero
				}else if(flags_arr(i,j,k).isSingleValued()){
					//fluxes remain zero
				}else if(flags_arr(i,j,k).isRegular()){
						
					LL = i-2;
					L = i-1;
					I = i;
					R = i+1;
					//-------------------convert S to W---------------------------
					if(MUSCL){
					v0[0] = S_old(LL,j,k,V["u"]);
					v0[1] = S_old(LL,j,k,V["v"]);
					Wll = Prim(S_old(LL,j,k,V["rho"]), v0, S_old(LL,j,k,V["p"]));
					v0[0] = S_old(R,j,k,V["u"]);
					v0[1] = S_old(R,j,k,V["v"]);
					Wr = Prim(S_old(R,j,k,V["rho"]), v0, S_old(R,j,k,V["p"]));
					}
					
					v0[0] = S_old(L,j,k,V["u"]);
					v0[1] = S_old(L,j,k,V["v"]);
					Wl = Prim(S_old(L,j,k,V["rho"]), v0, S_old(L,j,k,V["p"]));
					v0[0] = S_old(I,j,k,V["u"]);
					v0[1] = S_old(I,j,k,V["v"]);
					Wi = Prim(S_old(I,j,k,V["rho"]), v0, S_old(I,j,k,V["p"]));
					//-------------------convert S to W---------------------------
					
					
					//--------------------------------------MUSCL half time step----------------------------------------------
					if(MUSCL){
					delta_l = delta(omega, Wll, Wl, Wi, scheme);
					delta_i = delta(omega, Wl, Wi, Wr, scheme);
					
					//WbarL:
					WbarL.rho = Wl.rho + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.rho - dt/dx*Wl.rho*delta_l.u_vec[0]);
					WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[0] - dt/dx*((1/Wl.rho)*delta_l.p));
					WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[1]);
					WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[2]);
					S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
					//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
					WbarL.p = Wl.p + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.p - dt/dx*(Wl.rho*pow(aL,2))*delta_l.u_vec[0]);					
					
					//WbarR:
					WbarR.rho = Wi.rho - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.rho + dt/dx*Wi.rho*delta_i.u_vec[0]);
					WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[0] + dt/dx*((1/Wi.rho)*delta_i.p));
					WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[1]);
					WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[2]);
					S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
					//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
					WbarR.p = Wi.p - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.p + dt/dx*(Wi.rho*pow(aL,2))*delta_i.u_vec[0]);
					
					WL = WbarL; WR = WbarR;
					}else{
						WL = Wl; WR = Wi;
					}
					//--------------------------------------MUSCL half time step----------------------------------------------
					
					//direct wave speed estimates:
					S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
					//cout << "aL = " << aL << endl;
					S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
					SL = ((WL.u_vec[0] - aL < WR.u_vec[0] - aR)? WL.u_vec[0] - aL : WR.u_vec[0] - aR);
					SR = ((WL.u_vec[0] + aL > WR.u_vec[0] + aR)? WL.u_vec[0] + aL : WR.u_vec[0] + aR);
					Sstar = (WR.u_vec[0]*WR.rho*(SR-WR.u_vec[0]) - WL.u_vec[0]*WL.rho*(SL-WL.u_vec[0]) + (WL.p-WR.p))\
							/(WR.rho*(SR-WR.u_vec[0])-WL.rho*(SL-WL.u_vec[0]));
					Splus = ((fabs(WL.u_vec[0])+aL > fabs(WR.u_vec[0])+aR)? \
							  fabs(WL.u_vec[0])+aL : fabs(WR.u_vec[0])+aR);
							  
							  
					//W-->U and W-->F conversions:
					UL = conWtoU(WL);
					UR = conWtoU(WR);
					FL = conWtoF(WL);
					FR = conWtoF(WR);
					
					//flux evaluation at (x/t)=0:	
					if( 0 < SL){
						Fi = FL;
					}else if(SL <= 0 && 0 < Sstar){
						//cout << "left star state \n";
						//left star state terms:
						ULstar.rho = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar));
						ULstar.rhou_vec[0] = WL.rho*Sstar*((SL-WL.u_vec[0])/(SL-Sstar));
						ULstar.rhou_vec[1] = WL.rho*WL.u_vec[1]*((SL-WL.u_vec[0])/(SL-Sstar)); 
						ULstar.E = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar))*(UL.E/WL.rho+(Sstar-WL.u_vec[0])\
									*(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[0]))));
						if(ULstar.E < 0){
							cout << "negative energy computed in cell" << endl;
							amrex::Abort("aborted in MUSCL_HLLC_x");
							break;
						}
						//left star state flux:
						FLstar.mass = FL.mass + SL*(ULstar.rho - UL.rho);
						FLstar.mom_vec = FL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
						FLstar.en = FL.en + SL*(ULstar.E - UL.E);
						Fi = FLstar;
					}else if(Sstar <= 0 && 0 < SR){
						//cout << "right star state \n";
						//right star state terms:
						URstar.rho = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar));
						URstar.rhou_vec[0] = WR.rho*Sstar*((SR-WR.u_vec[0])/(SR-Sstar));
						URstar.rhou_vec[1] = WR.rho*WR.u_vec[1]*((SR-WR.u_vec[0])/(SR-Sstar)); 
						URstar.E = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar))*(UR.E/WR.rho+(Sstar-WR.u_vec[0])\
									*(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[0]))));
						if(ULstar.E < 0){
							cout << "negative energy computed in cell" << endl;
							amrex::Abort("aborted in MUSCL_HLLC_x");
							break;
						}
						FRstar.mass = FR.mass + SR*(URstar.rho - UR.rho);
						FRstar.mom_vec = FR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
						FRstar.en = FR.en + SR*(URstar.E - UR.E);
						Fi = FRstar;
					}else if(SR <= 0){
						Fi = FR; 
					}
					
					//-------------------convert F to flux_arr------------------------
					//overwrite calculated fluxes:
					F(i,j,k,V["mass"]) = Fi.mass;
					F(i,j,k,V["mom_x"]) = Fi.mom_vec[0];
					F(i,j,k,V["mom_y"]) = Fi.mom_vec[1];
					F(i,j,k,V["en"]) = Fi.en;
					//----------------------------------------------------------------		
				}								
			
			}
		}
	}//end i,j,k
      
}




