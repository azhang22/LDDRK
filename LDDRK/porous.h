//++++++++++++++++++++++++++++++filename: porous.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/**************************************************************************************/ 
/* eff_density is effective density
/* Kp is gas compressibility
/* IMAX,JMAX the max number of points for y and z direction
/* for one dimensition, set IMAX=1;therefore v_n is null;w_n and p_n is velocity 
/* and pressure for 1-D respectively
/***************************************************************************************/ 
#ifndef _POROUS_H
#define _POROUS_H
#include "mygeneral.h"
#include "calculation_2D.h"

class porous:public calculation_2D 
{
	protected:
		double eff_density,porosity,resistivity,surf_porosity,var_porosity,Kp;//kg/m3;mks.Rayls/m(or Pa.s/m2),m.s2/kg
		double Beta,Gama;
	public:
		porous(){};
		porous(scheme_list,DifferenceStep,char *coordi,MovingFrame MF,const double,PoreStruct PorePara,
				AirStruct AirPara,int mpi_rank1,int mpi_size1,int mpi_yarea1,int mpi_zarea1,int mpi_porous1);
		~porous();
		virtual void cal_fvw(double **temp_v,double **temp_w,
			double **temp_p);
		virtual void cal_fp(double **temp_v,double **temp_w,
			double **temp_p);
};
#endif






