#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
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
	

//EoS function:
double EoS_properties(std::string EoS, double rho, double p, double e, std::string prop, int EoS19){
	//EoS read from case: "ideal" or "EoS19"
	//properties: rho, p and e passed by reference
	//calculated property determined by strong prop: "e", "p", "a", "sigma_e"
	bool  err_flag = 0; 
	double gamma = 1.4;
	
	if(EoS == "ideal"){
		EoS19 = 0; //EOS19 equation of state not activated
		if(prop == "e"){
			e = p/(rho*(gamma-1));
			return e;
		}else if(prop == "p"){
			p = rho*(gamma - 1)*e;
			return p;
		}else if(prop == "a"){
			double a = pow(gamma*p/rho,0.5);
			return a;
		}else if(prop == "sigma_e"){
			//"can't produce conductivity vals from ideal EoS"
			return 0.0; 
		}else if(prop == "T"){
			return p/(rho*287);
		}else{
			assert(0 && "property function not defined" );
			return 0;
		}
	}else if(EoS == "EoS19"){
		EoS19 = 1; //EOS19 equation of state activated
		double out_val;
		return out_val;
		
	}else{
		assert(false && "EoS not specified or invalid");
        return 0;
	}
	
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
				e = EoS_properties("ideal", prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), 0, "e", 0);
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
				p = EoS_properties("ideal", prop_arr(i,j,k,V["rho"]), 0, e, "p", 0);
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
				e = EoS_properties("ideal", prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), 0, "e", 0);
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
    double e = EoS_properties("ideal", W.rho, W.p, 0, "e", 0);
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
    double p = EoS_properties("ideal", rho, 0, e, "p", 0);
    //double p = (gamma-1)*(U.E - 0.5*rho*pow(mag(u_vec),2));
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
    double e = EoS_properties("ideal", W.rho, W.p, 0, "e", 0);
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
    double e = EoS_properties("ideal", W.rho, W.p, 0, "e", 0);
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
    double e = EoS_properties("ideal", W.rho, W.p, 0, "e", 0);
    double E =  W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    double en = ((W.r > 0)? -1/W.r * W.u_vec[0]*(E+W.p) : 0);
    ConsF SE(mass, mom_vec, en);
    return SE;
}




#include "TestCases.H"

//some global variables and functions:
ofstream ufile, rhofile, pfile, efile, Bfile, Tfile, EoSfile, elecfile, namefile, MSfile;
Vector<Real> openfiles(TestCase, Prim);
void writetimestep(TestCase, const vector<Prim>&, Vector<Real>&);
void writeconductivities(TestCase, const vector<Prim>&);
void write_velocities(TestCase, const vector<Prim>&);
void writegeometry(TestCase);
void writecurrent(vector<double>&, vector<double>&);
void closefiles(TestCase);
Vector<Real> domain(10);






