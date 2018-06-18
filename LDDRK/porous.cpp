//+++++++++++++++++++++++++++filename: porous.cpp ++++++++++++++++++++++++++++++//

//-----------------------porous media acoustic wave equation---------------------//
/*******************************************************************************/ 
/* the program just apply Navier-Stokes equation to calculate acoustic pressure 
/* when sound propagation in the ground (porous media) with Eulerian time-domain model 
/*******************************************************************************/ 

//-----------------------scheme from Victor W.Sparrow--------------------------//
/*******************************************************************************/ 
/* Victor W.Sparrow
/* calculate pressure first, and then calculate velocity
/* apply non-staggered grid, velocity and pressure are at the same grid point
/* point pieces at y direction(air): 0,1,...IMAX-1;
/* point pieces at z direction(air): 0,1,...JMAX-1;
/********************************************************************************/

//---------------------------scheme from Erik M.Salomons-----------------------//
/*******************************************************************************/ 
/* Eulerian Time-Domain Model for sound propagation over 
/* a finite impedance ground surface. 
/* Comparison with frequency-domain models. 
/* Author: Erik M.Salomons.
/* ACTA Vol.88(2002) 483-492
/* calculate velocity first, and then calculate pressure
/* apply staggered grid, velocity and pressure are at the different grid point
/* point pieces at y direction(air): 0,1,...IMAX-1;
/* point pieces at z direction(air): 0,1,...JMAX-1;
/* east and west boundary point of v is v[1][j] and 
/* north and south boundary point of w is w[i][1] and w[i][JMAX-1]
/* the four boundary point v[1][j], v[IMAX-1][j], w[i][1], w[i][JMAX-1]
/* are calculated through equation, not from interpolation. This is different from
/* above colocated scheme
/*******************************************************************************/ 
#include <stdio.h>
#include <math.h>
#include "porous.h"

porous::porous(scheme_list scheme1,DifferenceStep Step,char *coordi,MovingFrame MF,
			   const double gauss_width,PoreStruct PorePara,AirStruct AirPara,
									int mpi_rank1,int mpi_size1,int mpi_yarea1,int mpi_zarea1,int mpi_porous1)
	:calculation_2D(scheme1,Step,coordi,MF,gauss_width,AirPara,mpi_rank1,mpi_size1,mpi_yarea1,mpi_zarea1,mpi_porous1)
{
	porosity=PorePara.porosity;
	eff_density=PorePara.eff_density*porosity;
	
	resistivity=PorePara.resistivity;
	Kp=PorePara.Kp;
	Beta=resistivity*diff_t/eff_density;
	Gama=diff_t/Kp;
	surf_porosity=porosity;
}

porous::~porous()
{
}

//----------------------------calculate velocity and pressure------------------//
void porous::cal_fvw(double **temp_v,double **temp_w,
			double **temp_p)
{
	int i,j;
	for(j=1;j<(mpi_JMAX[mpi_jindex]-1);j++){
		for(i=1;i<(mpi_IMAX[mpi_iindex]-1);i++){
			var_porosity=1.0*surf_porosity;

			v_nn[i][j]=-temp_v[i][j]-1.0/(1.0+Beta*var_porosity*Aav[j])*diff_t/eff_density/Jacobian[i][j]*Aav[j]*var_porosity
				*(Z_over_z[i][j]/diff_y*(temp_p[i][j]-temp_p[i-1][j])-
				Z_over_y[i][j]/diff_z*(temp_p[i][j]-temp_p[i][j-1]))
				+1.0/(1.0+Beta*Aav[j]*var_porosity)*temp_v[i][j];

			w_nn[i][j]=-temp_w[i][j]-1.0/(1.0+Beta*var_porosity*Aav[j])*diff_t/eff_density/Jacobian[i][j]*Aav[j]*var_porosity
				*(Y_over_y[i][j]/diff_z*(temp_p[i][j]-temp_p[i][j-1])-
				Y_over_z[i][j]/diff_y*(temp_p[i][j]-temp_p[i-1][j]))
				+1.0/(1.0+Beta*Aav[j]*var_porosity)*temp_w[i][j];
	}
}
}

void porous::cal_fp(double **temp_v,double **temp_w,
			double **temp_p)
{
	int i,j;

	for(j=1;j<(mpi_JMAX[mpi_jindex]-1);j++){
		for(i=1;i<(mpi_IMAX[mpi_iindex]-1);i++){
			var_porosity=1.0*surf_porosity;
			p_nn[i][j]=-(Gama/var_porosity)/Jacobian[i][j]*
				(Z_over_z[i][j]/diff_y*(temp_v[i+1][j]-temp_v[i][j])-
				Z_over_y[i][j]/diff_z*(temp_v[i][j+1]-temp_v[i][j])
				+Y_over_y[i][j]/diff_z*(temp_w[i][j+1]-temp_w[i][j])-
				Y_over_z[i][j]/diff_y*(temp_w[i+1][j]-temp_w[i][j]));
		}
	}
}
