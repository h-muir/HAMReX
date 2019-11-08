#ifndef FUNCDEFS_H_
#define FUNCDEFS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <string>

#include "classdefs.H"

void main_main ();

// in utility.cpp:
//operator functions 
void printVector(Vector<Real>);
Vector<Real> operator+(Vector<Real> const&, Vector<Real> const&);
Vector<Real> operator-(Vector<Real> const&, Vector<Real> const&);
Vector<Real> operator*(double const&, Vector<Real> const&);

double mag(Vector<Real>); 
double dot(Vector<Real>, Vector<Real>);

double EoS_properties(std::string, double, double, double, std::string, int);

//conversion functions:
void S_conWtoU(Box const&, Array4<Real> const&, AccessVariable&);
void S_conUtoW(Box const&, Array4<Real> const&, AccessVariable&);
void F_conWtoF(Box const&, Array4<Real> const&, AccessVariable&);

//class conversion functions:
ConsU conWtoU(Prim const&); ConsF conWtoF(Prim const&); ConsF conWtoG(Prim const&); 
ConsF conWtoSE(Prim const&); Prim conUtoW(ConsU const&);

//euler.cpp functions:
void initialiseStructs(SettingsStruct&, ParameterStruct&);
void initial_conditions(Box const&, Array4<Real> const&, const Real* dx, ParameterStruct const&, AccessVariable&);
void check_assumed_variables(Vector<std::string>);
void determine_grad_max(const Box&, Array4<Real> const&, const Real*, std::string, std::string const, Real&, AccessVariable&);
void amr_tagging(Array4<char> const&, const Box&, Array4<Real> const&, const Real*, 
					Real, Real, std::string, AccessVariable&);
double delta_bar_calc(std::string&, double&, double&, double&, double&);
Prim delta(double&, Prim&, Prim&, Prim&, std::string&);
void calc_fluxes_HLLC_x(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void calc_fluxes_HLLC_y(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void cons_update(const std::string, Box const&, Array4<Real> const&, Array4<Real> const&, 
        Array4<Real> const&, Real const, Real const, int const, AccessVariable&);
Real max_dim_wave_speed(Box const&, Array4<Real> const&, int const, AccessVariable&);
void rescale_fluxes(const std::string, Box const&, Array4<Real> const&, Real const, 
						Real const, int const, AccessVariable&);


#endif 
