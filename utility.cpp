#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB_levelset.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_LSCoreBase.H>

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <assert.h>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <limits>

#include "structdefs.H"
#include "classdefs.H"
#include "funcdefs.H"

using namespace amrex;
using namespace std;


//typedef vector<double> Vector;

//template<typename T>
void printVector(Vector<Real> vec){
    int len = vec.size();
    for(int i; i<len; ++i){
        cout << vec[i] << ' ';
    }
    cout << endl;
}//Vector-array print function 

//operator overloading for addition, subtraction and multiplication of vectors

//template<typename T>
Vector<Real> operator+(Vector<Real> const& first, Vector<Real> const& second)
{
    assert(first.size() == second.size());

    Vector<Real> result;
    result.reserve(first.size());

    transform(first.begin(), first.end(), second.begin(), back_inserter(result), plus<Real>());
    return result;
}

//template<typename T>
Vector<Real> operator-(Vector<Real> const& first, Vector<Real> const& second)
{
    assert(first.size() == second.size());

    Vector<Real> result;
    result.reserve(first.size());

    transform(first.begin(), first.end(), second.begin(), back_inserter(result), minus<Real>());
    return result;
}

//template<typename T>
Vector<Real> operator*(double const& scal, Vector<Real> const& first)
{
    
    vector<double> second(first.size(), scal);
    assert(first.size() == second.size());

    Vector<Real> result;
    result.reserve(first.size());

    transform(first.begin(), first.end(), second.begin(), back_inserter(result), multiplies<Real>());
    return result;
}



//vector functions:
double mag(Vector<Real> v){
    assert(v.size() == 3);
    double mag_v = pow( pow(v[0],2) + pow(v[1],2) + pow(v[2],2), 0.5);
    return mag_v;
}

double dot(Vector<Real> v1, Vector<Real> v2){
    assert(v1.size() == 3 && v2.size() == 3);
    double dot_v = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    return dot_v;
}

//accessvariable wrapper:

AccessVariable::AccessVariable(Vector<string> const &problem_variables){
	
	int num_vars = problem_variables.size();
	
	for(int n =0; n<num_vars; ++n){
		variables.insert(std::pair<string,int>(problem_variables[n],n));
	}
	
}

int& AccessVariable::operator[](string const &var)
{
	return variables[var];
}


//Embedded boundary stuff:


void makeEmbeddedBoundary(Geometry& geom, const int& max_level){
	
	ParmParse pp;
	string EBgeom;
	pp.get("EBgeom", EBgeom);
	ParmParse ppEB(EBgeom);
	int ls_pad = 1;		//number of cells away from EB in which level set is computed
	if(EBgeom == "sphere"){
		
		Real radius;
		Array<Real,AMREX_SPACEDIM> centre{0.5, 0.5}; //default Center of the sphere
		Vector<Real> EB_centre;
		EB_centre = Vector<Real>(BL_SPACEDIM);
		ppEB.get("R", radius);
		ppEB.getarr("centre", EB_centre);
		for(int i = 0; i<AMREX_SPACEDIM; ++i){
			centre[i] = EB_centre[i];
		}
		
		bool inside = false;  // Is the fluid inside the sphere --> false
		EB2::SphereIF sphere(radius, centre, inside);

		auto shop = EB2::makeShop(sphere);

		EB2::Build(shop, geom, max_level, max_level, ls_pad); 
		//arguments: (gshop, geom, required_coarsening_level, max_coarsening_level, ngrow)
	
	}else if(EBgeom == "plane"){
		Vector<Real> point, normal;
		ppEB.getarr("point", point);
		ppEB.getarr("normal", normal);
		
		EB2::PlaneIF plane({AMREX_D_DECL(point[0], point[1], 0.)},
                             {AMREX_D_DECL(normal[0], normal[1], 0.)}, false);
        
        auto shop = EB2::makeShop(plane);

		EB2::Build(shop, geom, max_level, max_level, ls_pad); 
		//arguments: (gshop, geom, required_coarsening_level, max_coarsening_level, ngrow)                     
                             
		
	}else if(EBgeom == "none"){
		//no additional settings
		cout << "EBgeom is set to \"none\" \n";
		//EB2 build still required for valid intitialisation under EB = TRUE in makefile
		EB2::Build(geom, max_level, max_level, ls_pad);
		
	}else{
		amrex::Abort("invalid EBgeom passed from inputs file");
	}
	
}

void initialise_levelset_geometry(const int& level, Geometry& geom, const BoxArray& ba, 
						const DistributionMapping& dm, SettingsStruct const& sim, ParameterStruct const& p, 
						  AccessVariable& V, MultiFab& S_new)
{
	int ls_pad = 4; //number of cells away from EB in which level set is computed
	int eb_pad = 2; //sim.NGROW;
	int ls_ref = 1; //pow(2, level);
	int eb_ref = 1;
	int ebt_size = 32;
	
	if(sim.EBgeom == "sphere"){
		
	}else if(sim.EBgeom == "plane"){
		//Define EB:
		EB2::PlaneIF plane({AMREX_D_DECL(p.EB_point[0], p.EB_point[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal[0], p.EB_normal[1], 0.)}, false);
        EB2::GeometryShop<EB2::PlaneIF> plane_gshop(plane);
        
        // Build level-set factory
		LSFactory level_set(level, ls_ref, eb_ref, ls_pad, eb_pad, ba, geom, dm);

		// Build EB
		const Geometry& eb_geom = level_set.get_eb_geom();
		EB2::Build(plane_gshop, eb_geom, p.max_level, p.max_level);
		//parameters: gshop, eb_geom, required_coarsening_level, max_coarsening_level
		
		const EB2::IndexSpace& plane_ebis = EB2::IndexSpace::top();
		const EB2::Level&      plane_lev  = plane_ebis.getLevel(eb_geom);

		// Build EB factory
		EBFArrayBoxFactory eb_factory(plane_lev, eb_geom, ba, dm, {eb_pad, eb_pad, eb_pad}, EBSupport::full);
		//const Vector<int> ngrow_vec{eb_pad,eb_pad,eb_pad};
		//std::unique_ptr<amrex::EBFArrayBoxFactory> factory_ptr = makeEBFabFactory(geom, S_new.boxArray(), dm, ngrow_vec, EBSupport::full);
		//auto eb_factory = dynamic_cast<EBFArrayBoxFactory const*>(&(S_new.Factory()));  
		
		// Fill level-set (factory)
		GShopLSFactory<EB2::PlaneIF> plane_lsgs(plane_gshop, level_set);
		std::unique_ptr<MultiFab> plane_mf_impfunc = plane_lsgs.fill_impfunc();

		MultiFab::Copy(S_new, *plane_mf_impfunc, 0, V["level_set"], 1, S_new.nGrow());
			
	}else{
		amrex::Abort("invalid EBgeom setting");
	}
	
}

void initialise_levelset_normals(Box const& box, Array4<Real> const& S_old, Array4<Real> const& S_new, 
									const Real* dx, AccessVariable& V)
{
    
    const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real nx, ny;
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
				nx = (S_old(i+1,j,k,V["level_set"]) - S_old(i-1,j,k,V["level_set"]))/(2*dx[0]);
				ny = (S_old(i,j+1,k,V["level_set"]) - S_old(i,j-1,k,V["level_set"]))/(2*dx[1]);
				S_new(i,j,k,V["nx"]) = nx/pow(pow(nx,2) + pow(ny,2),0.5);
				S_new(i,j,k,V["ny"]) = ny/pow(pow(nx,2) + pow(ny,2),0.5);
			}
		}
	} // end i,j,k loop
	
}



void constructLevelSet(const int& level, Geometry& geom, const BoxArray& ba, 
						const DistributionMapping& dm, SettingsStruct const& sim, ParameterStruct const& p, 
						  MultiFab& LSMF, iMultiFab& ivalid_LSMF){
	
	//update LSMF and ivalid_LSMF by reference (ncomp = 1)
	
	int ls_pad = 1; //number of cells away from EB in which level set is computed
	int eb_pad = 1; //sim.NGROW;
	int ls_ref = 1;
	int eb_ref = 1;
	//const IntVect ebt_size(AMREX_D_DECL(64, 64, 64));
	//const IntVect ebt_size(AMREX_D_DECL(eb_pad, eb_pad, eb_pad)); //EB tile size  ???? check this
	int ebt_size = 32;
	
	
	if(sim.EBgeom == "sphere"){
		/*
		Array<Real,AMREX_SPACEDIM> centre{p.EB_centre[0], p.EB_centre[1]};
		bool inside = false;  // Is the fluid inside the sphere --> false
		EB2::SphereIF sphere(p.EB_R, centre, inside);
		EB2::GeometryShop<EB2::SphereIF> sphere_gshop(sphere);
		//EB2::Build(sphere_gshop, geom, p.max_level, p.max_level, ls_pad); 
		//arguments: (gshop, geom, required_coarsening_level, max_coarsening_level, ngrow)
	
		//const EB2::IndexSpace& EB_indexspace = EB2::IndexSpace::top();
		//const EB2::Level&      EB_level  = EB_indexspace.getLevel(geom);

		//build eb_factor:
		//EBFArrayBoxFactory eb_factory(EB_level, geom, ba, dm, {eb_pad, eb_pad, eb_pad}, EBSupport::full);
		// Fill implicit function
		GShopLSFactory<EB2::SphereIF> EB_levelset(sphere_gshop, geom, ba, dm, ls_pad);
		std::unique_ptr<MultiFab> EB_impfuncMF = EB_levelset.fill_impfunc();
		LSFactory(level, 1, 1, ls_pad, eb_pad, ba, geom, dm, ebt_size);
		// Fill level-set:
		//LSFactory::fill_data(LSMF, ivalid_LSMF, eb_factory, *EB_impfuncMF,
			//			 ebt_size, 1, 1, geom, geom);  //1,1 are level set and EB refinement factors: 
													   // ie. no additional refinement here
		//alternate function call:
		//LSFactory::fill_data(LSMF, ivalid_LSMF, *EB_impfuncMF, eb_pad, geom);  
		
		//copy data from ls multipfabs created in function --> multifabs passed into function:
		MultiFab::Copy(LSMF, *EB_impfuncMF, 0, 0, 1, ls_pad);
		
		//copy data from ls multipfabs created in function --> multifabs passed into function:
		//MultiFab::Copy(LSMF, ls_grid, 0, 0, 1, ls_pad); //check ls_pad as last parameter here
		*/
                     
    }else if(sim.EBgeom == "plane"){
		//Define EB:
		EB2::PlaneIF plane({AMREX_D_DECL(p.EB_point[0], p.EB_point[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal[0], p.EB_normal[1], 0.)}, false);
        EB2::GeometryShop<EB2::PlaneIF> plane_gshop(plane);
        
        // Build level-set factory
		LSFactory level_set(level, ls_ref, eb_ref, ls_pad, eb_pad, ba, geom, dm);
		//std::unique_ptr<iMultiFab> ilevel_set_ptr = level_set.copy_valid(dm);
		//iMultiFab ilevel_set(ba, dm, 1, ls_pad);
		//iMultiFab::Copy(ilevel_set, *ilevel_set_ptr, 0, 0, 1, ls_pad);

		// Build EB
		const Geometry& eb_geom = level_set.get_eb_geom();
		EB2::Build(plane_gshop, eb_geom, p.max_level, p.max_level);
		
		const EB2::IndexSpace& plane_ebis = EB2::IndexSpace::top();
		const EB2::Level&      plane_lev  = plane_ebis.getLevel(eb_geom);

		// Build EB factory
		EBFArrayBoxFactory eb_factory(plane_lev, eb_geom, ba, dm, {eb_pad, eb_pad, eb_pad}, EBSupport::full);
		
		// Fill level-set (factory)
		GShopLSFactory<EB2::PlaneIF> plane_lsgs(plane_gshop, level_set);
		std::unique_ptr<MultiFab> plane_mf_impfunc = plane_lsgs.fill_impfunc();
		//MultiFab plane_mf(ba, dm, 1, ls_pad);  
		//cout << "FLAG 1 \n";
		//MultiFab::Copy(plane_mf, *plane_mf_impfunc, 0, 0, 1, ls_pad);
		//amrex::Print() << "eb_factory.nComp, eb_factory.nGrow : " << eb_factory.nComp() << ", " << eb_factory.nGrow << "\n";
		//amrex::Print() << "plane_mf_impfunc.nComp, plane_mf_impfunc.nGrow : " << plane_mf.nComp() << ", " << plane_mf.nGrow() << "\n";
		//level_set.Fill(eb_factory, *plane_mf_impfunc);  //segfault here
		//eb_factory.fillLevelSet(LSMF, eb_geom);
		
		//cout << "FLAG 2 \n";
		//MultiFab::Copy(LSMF, *level_set.get_data(), 0, 0, 1, ls_pad);
		MultiFab::Copy(LSMF, *plane_mf_impfunc, 0, 0, 1, ls_pad);
		//cout << "FLAG 3 \n";
		//LSFactory::fill_data(LSMF, ivalid_LSMF, eb_factory, *plane_mf_impfunc, ebt_size, ls_ref, eb_ref, geom, geom);  //segfault
		//LSFactory::fill_data(LSMF, ivalid_LSMF, *plane_mf_impfunc, eb_pad, geom); //segfault
			
	}else{
		amrex::Abort("invalid EBgeom setting");
	}
	
}
	
void find_embedded_boundary(const Box& box, Array4<Real> const& S_new_arr, Array4<const EBCellFlag>& flags_arr, const int& nc, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    
    //find cut cells, covered cells and regular cells:
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				
				if(flags_arr(i,j,k).isCovered()){
					//cout << "cell is covered \n";
					//for(int n = 0; n<nc; ++n){
					//	S_new_arr(i,j,k,n) = 0;
					//}
					S_new_arr(i,j,k, V["EB"]) = -1;
				}else if(flags_arr(i,j,k).isSingleValued()){
					//cout << "cell is cut \n";
					S_new_arr(i,j,k, V["EB"]) = 0;
				}else if(flags_arr(i,j,k).isRegular()){
					//cout << "cell is regular \n";
					S_new_arr(i,j,k, V["EB"]) = 1;
				}

			}
		}
	}//end i,j,k loop
	
}

//interpolation functions for Plasma19 EoS:

int find_index(const vector<double>& vec_ref, double &val1, double &val2, double val, bool &err_flag){
	int index;
	int vec_length = vec_ref.size();
	
	//cout << "vec_length: " << vec_length << endl;
	
	if( (val < vec_ref[0]) || (val > vec_ref[vec_length-1]) ){
		err_flag = true;
		//amrex::Abort("Plasma19 OOB: value passed to find_index function is out of range (low)");
		return 0;
		//cout << "throwing interp error" << endl;
		//throw "interp_error";
	}else{
		//assumes vec_ref in ascending value order:
		//assert(val > vec_ref[0] \
				&& "EoS input property is not within valid tabulated data range");
		//assert(val < vec_ref[vec_length-1] \
				&& "EoS input property is not within valid tabulated data range");
		
		index = 0;
		while(val >= vec_ref[index]){
			++index;
			if(index >= vec_length-1){
				cout << "value not in range \n";
				//amrex::Abort("Plasma19 OOB: value passed to find_index function is out of range (high)");
				break;
			}
		}
		--index;
		val1 = vec_ref[index];
		val2 = vec_ref[index+1];
		
		return index;
	}
}

double interp_1Dvec(const vector<double>& vec_ref, const vector<double>& vec_prop, double val, bool &err_flag){
		
		double frac;
		if(vec_ref.size() != vec_prop.size()){
			amrex::Abort("Plasma19 OOB: interp_1Dvec (utility.cpp) inconsistent vector lengths");
		}
		
		double val1;
		double val2;
		int index = find_index(vec_ref, val1, val2, val, err_flag);
		frac = (val - val1)/(val2 - val1);	
		
		double out_val = vec_prop[index] + frac*(vec_prop[index+1]-vec_prop[index]);
		//output by variable passed by reference
		
		return out_val;	
}

double interp_finalval(double ref1, double ref2, double val1, double val2, double ref_val, bool &err_flag){
	
	if( !((ref1 <= ref_val && ref_val <= ref2) || (ref2 <= ref_val && ref_val <= ref1)) ){
		err_flag = true;
		amrex::Abort("Plasma19 OOB: value passed to interp_finalval function is out of range");
		return 0;
		//cout << "throwing interp error" << endl;
		//throw "interp_error";
	}else{	
		//cout << "vals: " << ref1 << ' ' << ref2 << ' ' << ref_val << endl;
		//assert((ref1 <= ref_val && ref_val <= ref2) || (ref2 <= ref_val && ref_val <= ref1)\
				&& "EoS input property is not within valid tabulated data range");
		double out_val = val1 + (ref_val-ref1)/(ref2-ref1) * (val2-val1);
		
		return out_val;
	}
}

double interp2D(const vector< vector<double> >& m_OutputProperty, 
				const vector<double>& v_Pressures, const vector<double>& v_Densities, 
				const double rho_val, const double p_val, double e_val, bool &err_flag){
	
	//e,a=f(rho, p)
	//get p1 p2 val and index:
	double p1_val;
	double p2_val;

	int p1 = find_index(v_Pressures, p1_val, p2_val, p_val, err_flag);
	int p2 = p1+1;
	
	//check if OOB value has been encountered
	if(err_flag){
		return 0;
	}
	
	//get corresponding output vals from rho interp of p1 row and p2 row
	double out_val1 = interp_1Dvec(v_Densities, m_OutputProperty[p1], rho_val, err_flag);
	double out_val2 = interp_1Dvec(v_Densities, m_OutputProperty[p2], rho_val, err_flag);
						
	//obtain final 2D interp value
	double final_val = interp_finalval(p1_val, p2_val, out_val1, out_val2, p_val, err_flag);
	
	return final_val;
}

double interp2D_p(const vector< vector<double> >& m_InternalEnergies, 
				const vector<double>& v_Pressures, const vector<double>& v_Densities, 
				const double rho_val, double p_val, const double e_val, bool &err_flag){
	
	//p=f(rho, e)
	//can generalise later (obtain p or rho from another matrix input property)
	double rho1_val;
	double rho2_val;
	int rho1 = find_index(v_Densities, rho1_val, rho2_val, rho_val, err_flag);
	int rho2 = rho1+1;
	
	//contruct row vectors from output property columns
	//corresponding to rho values
	int vec_length = v_Densities.size();
	vector<double> eVector1(vec_length);
	vector<double> eVector2(vec_length);
	for(int i = 0; i<vec_length; ++i){
		eVector1[i] = m_InternalEnergies[i][rho1];
		eVector2[i] = m_InternalEnergies[i][rho2];
	}
		
	//get corresponding output vals from rho interp of p1 row and p2 row
	double out_val1 = interp_1Dvec(eVector1, v_Pressures, e_val, err_flag);
	double out_val2 = interp_1Dvec(eVector2, v_Pressures, e_val, err_flag);
						
	//obtain final 2D interp value
	double final_val = interp_finalval(rho1_val, rho2_val, out_val1, out_val2, rho_val, err_flag);
	
	return final_val;
}








