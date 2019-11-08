#include "structdefs.H"
#include "funcdefs.H"
#include "classdefs.H"
#include <cmath>
#include <algorithm>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_TagBox.H>

using namespace amrex;
using namespace std;

/* ---------------------GLOBAL VARIABLES: ----------------------------*/
#define gamma 1.4   		//heat capacity ratio for standard air
int nT; 					//final number of simulation time steps
std::string EoS = "ideal";		//system Equation of State
int EoS19 = 0;				//whether EoS19 is included
/*--------------------------------------------------------------------*/

//the following variable name sare assumed for the current system:
Vector<string> assumed_variables{"mass", "mom_x", "mom_y", "en", "rho", "u", "v", "p"};

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
    pp.get("MUSCL", SimSettings.MUSCL);
    pp.get("rigid_body", SimSettings.rigid_body);
    if(SimSettings.rigid_body){
		pp.get("rigid_geom", SimSettings.rigid_geom);
		ParmParse ppRigid(SimSettings.rigid_geom);
		if(SimSettings.rigid_geom == "sphere"){
			ppRigid.get("radius", Parameters.radius);
			ppRigid.getarr("centre", Parameters.centre);
		}
	}
    
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
    if(Parameters.IC == "x0" || Parameters.IC == "y0"){ 
		ppTest.get("x0", Parameters.x0);
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
    
}


//problem initialisation:
void initial_conditions(Box const& box, Array4<Real> const& prop_arr, const Real* dx, ParameterStruct const& p, AccessVariable& V){

    const auto lo = lbound(box);
    const auto hi = ubound(box);
	
	//amrex::Print() << "lo : " << lo << ", hi : " << hi << "\n";

	//will need changing to *dx when cell sizes change in AMR
    int int_x0 = (p.x0/p.prob_hi[0])*p.n_cells[0];
    int int_y0 = (p.x0/p.prob_hi[1])*p.n_cells[1];
	
	int nc = p.NCOMP;
	
	//AccessVariable V(p.problem_variables);
	
	//cout << "PL and PR size() " << p.PL.size() << " : " << p.PR.size() << endl;
	
	Real R = 0.2;
	Real h;
	
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
	}else{
		amrex::Abort("invalid intitial condition - check inputs file for adv.IC");
	}
	
	S_conWtoU(box, prop_arr, V);
	
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
	string scheme = "UltraBee";
	
	
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
				aL = EoS_properties(EoS, Wl.rho, Wl.p, 0, "a", EoS19);
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.p - dt/dx*(Wl.rho*pow(aL,2))*delta_l.u_vec[0]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.rho + dt/dx*Wi.rho*delta_i.u_vec[0]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[0] + dt/dx*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[1]);
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[2]);
				aR = EoS_properties(EoS, Wi.rho, Wi.p, 0, "a", EoS19);
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.p + dt/dx*(Wi.rho*pow(aL,2))*delta_i.u_vec[0]);
				
				WL = WbarL; WR = WbarR;
				}else{
					WL = Wl; WR = Wi;
				}
				//--------------------------------------MUSCL half time step----------------------------------------------
				
				//direct wave speed estimates:
				aL = EoS_properties(EoS, WL.rho, WL.p, 0, "a", EoS19);
				//cout << "aL = " << aL << endl;
				aR = EoS_properties(EoS, WR.rho, WR.p, 0, "a", EoS19);
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
				F(i,j,k,V["mass"]) = Fi.mass;
				F(i,j,k,V["mom_x"]) = Fi.mom_vec[0];
				F(i,j,k,V["mom_y"]) = Fi.mom_vec[1];
				F(i,j,k,V["en"]) = Fi.en;
				F(i,j,k,V["rho"]) = 0.0;
				F(i,j,k,V["u"]) = 0.0;
				F(i,j,k,V["v"]) = 0.0;
				F(i,j,k,V["p"]) = 0.0;
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
	string scheme = "UltraBee";
	
	
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
				aL = EoS_properties(EoS, Wl.rho, Wl.p, 0, "a", EoS19);
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.p - dt/dy*(Wl.rho*pow(aL,2))*delta_l.u_vec[1]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.rho + dt/dy*Wi.rho*delta_i.u_vec[1]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[0]);
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[1] + dt/dy*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[2]);
				aR = EoS_properties(EoS, Wi.rho, Wi.p, 0, "a", EoS19);
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.p + dt/dy*(Wi.rho*pow(aL,2))*delta_i.u_vec[1]);
				
				WL = WbarL; WR = WbarR;
				}else{
				WL = Wl; WR = Wi;
				}
				//--------------------------------------MUSCL half time step----------------------------------------------
				//direct wave speed estimates:
				aL = EoS_properties(EoS, WL.rho, WL.p, 0, "a", EoS19);
				//cout << "aL = " << aL << endl;
				aR = EoS_properties(EoS, WR.rho, WR.p, 0, "a", EoS19);
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
				G(i,j,k,V["mass"]) = Gi.mass;
				G(i,j,k,V["mom_x"]) = Gi.mom_vec[0];
				G(i,j,k,V["mom_y"]) = Gi.mom_vec[1];
				G(i,j,k,V["en"]) = Gi.en;
				G(i,j,k,V["rho"]) = 0.0;
				G(i,j,k,V["u"]) = 0.0;
				G(i,j,k,V["v"]) = 0.0;
				G(i,j,k,V["p"]) = 0.0;
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

Real max_dim_wave_speed(Box const& bx, Array4<Real> const& S_new, int const i, AccessVariable& V){
	
	//convert from i dimension to n_index
	
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
				
				sound_speed = fabs(EoS_properties(EoS, S_new(i,j,k,V["rho"]), S_new(i,j,k,V["p"]), 0, "a", EoS19)); 
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
					string prop, AccessVariable& V)
{
	//tagging function for refinement
	//based on gradient of rho
	//assumes rho is ncomp element 0
	//box has S_new dimensions, but S_old_arr has at least +1 ngrow
	//therefore, we can use central different function for gradient
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real grad_rho, phi_x, phi_y;
    //Real prop_max = 0.5;
    
	//set tags based on gradient relative to max gradient: 
	//ie. prop_max is proportion of max gradient
	for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				phi_x = fabs((S_old_arr(i+1,j,k,V[prop])-S_old_arr(i-1,j,k,V[prop]))/(2*dx[0]));
				phi_y = fabs((S_old_arr(i,j+1,k,V[prop])-S_old_arr(i,j-1,k,V[prop]))/(2*dx[1]));
				grad_rho = sqrt(pow(phi_x,2) + pow(phi_y,2));
				if(grad_rho > grad_frac*grad_max){
					tagarr(i,j,k,V[prop]) = TagBox::SET;
					//cout << "tagged" << endl;
				}else{
					tagarr(i,j,k,V[prop]) = TagBox::CLEAR;
				}
			}
		}
	}//end i,j,k
		
}
