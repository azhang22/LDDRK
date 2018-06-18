//++++++++++++++++++++++++++filename: air.cpp +++++++++++++++++++++++++++++++++//

//-----------------------air acoustic wave equation--------------------------//
/*******************************************************************************/ 
/* the program just apply Navier-Stokes equation to calculate acoustic pressure 
/* when sound propagation with Eulerian time-domain model 
/*******************************************************************************/ 

//-----------------------scheme from Victor W.Sparrow--------------------------//
/*******************************************************************************/ 
/* Victor W.Sparrow
/* calculate pressure first, and then calculate velocity
/* apply non-staggered grid(colocated), velocity and pressure are at the same grid point
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
/* east and west boundary point of v is v[1][j] and v[IMAX-1][j]
/* north and south boundary point of w is w[i][1] and w[i][JMAX-1]
/* the four boundary point v[1][j], v[IMAX-1][j], w[i][1], w[i][JMAX-1]
/* are calculated through equation, not from interpolation. This is different from
/* above colocated scheme
/*******************************************************************************/ 
#include <stdio.h>
#include <math.h>
#include "air.h"
#include "mygeneral.h"

SpeedProfile::SpeedProfile(double velocity_coef,int VM) //double velocity_coef,double sound_speed,int VM
{
	
	alpha=1.256431;
	b=velocity_coef;
	velocity_method=VM;
	/*
	c=sound_speed;velocity_method=VM;
	switch(velocity_method){
	case 2:
		{
			y0=0;z0=0;L=0.03;
			vorticity=b*(1+2*alpha)*c/L;
			circulation=vorticity*L*L*PI/alpha;
		}
		break;
	case 3:
		{
			y0=0.66;z0=0;L=0.03;
			circulation=b;
		}
		break;
	case 6:
		{
			y0=3.5;z0=4;L=1;
			vorticity=b*(1+2*alpha)*c/L;
			circulation=vorticity*L*L*PI/alpha;
		}
		break;
	default:
		;
	}
	*/
  }

double SpeedProfile::Vav(double x,double y,double z)
{
	switch(velocity_method){
	case 1://"Eulerian time-domain ...",Acta Acustica united with Acustica,Vol.88(2002) 483-492
		{
			//////only for H=4////
			if(z>=8.0){
				return b*8.0;
			}else{
				return b*z;
			}
		}
		break;
	case 2://"study of the sound-vortex interaction ...", Eur.Phys.J.B 32,237-242(2003)
		{
			/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.000001) r=1e-20;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return (z-z0)/r*Vr;
			*/
		  }
		break;
	case 3://"transmission of sound through a single vortex", Eur.Phys.J.B 37,229-239(2004)
		{
			/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if(r<=L){//Vr=circulation/(2*PI*L)*r/L;Vav=(z-z0)/r*Vr 
				return circulation/(2*PI*L)/L*(z-z0);
			}
			Vr=circulation/(2*PI)/r;
			return (z-z0)/r*Vr;
			*/	
		 }
		break;
	case 4://reverse coordinate from case 1
		{
			return 0;
		}
		break;
	case 5://book of Solomans: terrain
		{
			return b*log(fabs(z)/0.1+1);
		}
		break;
	case 6:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.00001) return 0;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return (z-z0)/r*Vr;
		*/
		}
		break;
	default:
		return 0;
	}
}

double SpeedProfile::Wav(double x,double y,double z)
{
	switch(velocity_method){
	case 1:
		{
			return 0;
		}
		break;
	case 2:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.000001) r=1e-10;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return -(y-y0)/r*Vr;
		*/		
		 }
		break;
	case 3:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if(r<=L){//Vr=circulation/(2*PI*L)*r/L;Wav=-(y-y0)/r*Vr;
				return circulation/(2*PI*L)/L*(-(y-y0));
			}
			Vr=circulation/(2*PI)/r;
			return -(y-y0)/r*Vr;
		*/
		}
		break;
	case 4:
		{
			return 0;
		}
		break;
	case 5:
		{
			return b*y;
		}
		break;
	case 6:
		{
		/*
			double r,Vr;
			r=sqrt(pow((y-y0),2)+pow((z-z0),2));
			if((r/L)<=0.00001) return 0;
			Vr=circulation/(2*PI*r)*(1-exp(-alpha*r*r/L/L));
			return -(y-y0)/r*Vr;
		*/
		}
		break;
	default:
		return 0;//b*log(fabs(27-y)/0.1+1)
	}
}
SpeedProfile::~SpeedProfile(){}
//---------------------------------air member function---------------------------------------//

air::air(scheme_list scheme1,DifferenceStep Step,char *coordi,MovingFrame MF,const double gauss_width,
		 double velocity_coef,int velocity_method1,AirStruct AirPara,double AbsorbZmax,PoreStruct PorePara,
		 hillpore hill1,int mpi_rank1,int mpi_size1,int mpi_yarea1,int mpi_zarea1,int mpi_porous1)
   :calculation_2D(scheme1,Step,coordi,MF,gauss_width,AirPara,mpi_rank1,mpi_size1,mpi_yarea1,mpi_zarea1,mpi_porous1)
{
	int i,j;
	adiabatic_coef=AirPara.adiabatic_coef;
	Pav=AirPara.Pav; //Aav=AirPara.Aav,sound_speed=AirPara.sound_speed;
	coef_tao=AirPara.coef_tao;
	coef_eta=AirPara.coef_eta;
	coef_y=AirPara.coef_y;
	velocity_method=velocity_method1;
	SP=new SpeedProfile(velocity_coef,velocity_method); //velocity_coef,AirPara.sound_speed,velocity_method
	//calculate partial differential of wind speed
	Vav=new double* [mpi_IMAX[mpi_iindex]];Wav=new double* [mpi_IMAX[mpi_iindex]];
	Vav_over_y=new double* [mpi_IMAX[mpi_iindex]];Vav_over_z=new double* [mpi_IMAX[mpi_iindex]];
	Wav_over_y=new double* [mpi_IMAX[mpi_iindex]];Wav_over_z=new double* [mpi_IMAX[mpi_iindex]];
	coef1=new double* [mpi_IMAX[mpi_iindex]];coef2=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		Vav[i]=new double [mpi_JMAX[mpi_jindex]];Wav[i]=new double [mpi_JMAX[mpi_jindex]];
		Vav_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Vav_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
		Wav_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Wav_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
		coef1[i]=new double [mpi_JMAX[mpi_jindex]];coef2[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	coef_weighty=new double *[2];
	coef_weightz=new double *[2];
	for (i=0;i<2;i++){
		coef_weighty[i]=new double [mpi_IMAX[mpi_iindex]];
		coef_weightz[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	for (i=0;i<mpi_IMAX[mpi_iindex];i++){
		if (i==0){
			coef_weighty[0][i]=1.0;
			coef_weighty[1][i]=1.0;
		}else{
			coef_weighty[0][i]=cal_weight(coef_y,i);
			coef_weighty[1][i]=cal_weight(coef_y-1.0,i);
		}
	}
	for (j=0;j<mpi_JMAX[mpi_jindex];j++){
		if (j==0){
			coef_weightz[0][j]=1.0;
			coef_weightz[1][j]=1.0;
		}else{
			coef_weightz[0][j]=cal_weight(coef_y,j);
			coef_weightz[1][j]=cal_weight(coef_y-1.0,j);
		}
	}
	pow_y1=1.0/pow(diff_y,coef_y);
	pow_z1=1.0/pow(diff_z,coef_y);
	pow_y2=1.0/pow(diff_y,coef_y-1.0);
	pow_z2=1.0/pow(diff_z,coef_y-1.0);
	var_miu1=0.5/cos(3.1415926*coef_y*0.5);
	var_miu2=0.5/cos(3.1415926*(coef_y-1.0)*0.5);

	Qv=new double*[mpi_IMAX[mpi_iindex]];Qw=new double *[mpi_IMAX[mpi_iindex]];Qp=new double *[mpi_IMAX[mpi_iindex]];
	for (i=0;i<mpi_IMAX[mpi_iindex];i++){
		Qv[i]=new double[mpi_JMAX[mpi_jindex]];
		Qw[i]=new double[mpi_JMAX[mpi_jindex]];
		Qp[i]=new double[mpi_JMAX[mpi_jindex]];
	}
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			Qv[i][j]=0;Qw[i][j]=0;Qp[i][j]=0;
		}
	}
	whole_Qv=new double*[IMAX];whole_Qw=new double*[IMAX];whole_Qp=new double*[IMAX];
	for (i=0;i<IMAX;i++){
		whole_Qv[i]=new double[mpi_JMAX[mpi_jindex]];
		whole_Qw[i]=new double[mpi_JMAX[mpi_jindex]];
		whole_Qp[i]=new double[mpi_JMAX[mpi_jindex]];
	}
	for(i=0;i<IMAX;i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			whole_Qv[i][j]=0;whole_Qw[i][j]=0;whole_Qp[i][j]=0;
		}
	}

	PML_AbsorbZmax=AbsorbZmax;PML_alpha=2.0;
	if(JMAX1!=0){
		PML_Height1=fabs(leftbound_Z[JMAX1-1]-leftbound_Z[0]);
		PML_z1=leftbound_Z[JMAX1-1];
	}
	if(JMAX2!=0){
		PML_Height2=fabs(leftbound_Z[JMAX-1]-leftbound_Z[JMAX-JMAX2]);
		PML_z2=leftbound_Z[JMAX-JMAX2];
	}
	// PML for y direction
	PML_AbsorbYmax = AbsorbZmax; PML_alphay = PML_alpha;
	if (IMAX2 != 0){ ////////////////// PML RH
		PML_y1 = bottombound_Y[IMAX2 - 1];//bottombound_Y[0];
		PML_y2 = bottombound_Y[IMAX - IMAX2];
		PML_Width1 = fabs(bottombound_Y[IMAX2 - 1] - bottombound_Y[0]);// fabs(bottombound_Y[0] - bottombound_Y[0]);
		PML_Width2 = fabs(bottombound_Y[IMAX - 1] - bottombound_Y[IMAX - IMAX2]);
		//PML_AbsorbYmax = 1e3;
		//PML_alphay = PML_alpha;
		if (mpi_rank == 0) cout << "PML_y1 = " << PML_y1 << ", PML_y2 = " << PML_y2 << ", PML_Width1 = " << PML_Width1 << ", PML_Width2 = " << PML_Width2 << endl;
	}
	// if cr_judge is 1, porous media hill exists there
	cr_judge=hill1.curve_judge;
	ystar=hill1.curve_ystar;
	ystop=hill1.curve_ystop;
	zstar=hill1.curve_zstar;
	zstop=hill1.curve_zstop;
	cr_a=hill1.curve_coefa;
	cr_b=hill1.curve_coefb;
	cr_c=hill1.curve_coefc;
	cr_d=hill1.curve_coefd;
	cr_e=hill1.curve_coefe;
	cr_f=hill1.curve_coeff;
	cr_g=hill1.curve_coefg;
	cr_h=hill1.curve_coefh;
	cr_i=hill1.curve_coefi;
	cr_j=hill1.curve_coefj;
	cr_k=hill1.curve_coefk;
	cr_l=hill1.curve_coefl;
	cr_points=hill1.curve_points;

	crpo_eff_density=PorePara.eff_density;
	crpo_porosity=PorePara.porosity;
	crpo_resistivity=PorePara.resistivity;
	crpo_Kp=PorePara.Kp;
// crpo_reistivity should be a variable to be chose	
	crpo_Beta=crpo_resistivity*diff_t/crpo_eff_density; 
	crpo_Gama=diff_t/crpo_Kp/crpo_porosity;
	
	cr_y=new double [cr_points]; cr_z=new double [cr_points]; 
	cr_uc=new int*[IMAX];
	for (i=0;i<IMAX;i++) cr_uc[i]=new int[JMAX];
	
	for(i=0;i<IMAX;i++){
		for (j=0;j<JMAX;j++) cr_uc[i][j] = 0.0;
	}
	
	Ny = new int[JMAX];	Nz = new int[IMAX]; Ny1 = new int[JMAX]; Nz1 = new int[IMAX];
	Sy = new int[JMAX];	Sz = new int[IMAX]; Sy1 = new int[JMAX]; Sz1 = new int[IMAX];
	zone_y = new int[JMAX]; zone_z = new int[IMAX];
	for (i = 0; i < IMAX; i++) { Nz[i] = 0; Sz[i] = 0; Nz1[i] = 0; Sz1[i] = 0; zone_z[i] = 0; }
	for (j = 0; j < JMAX; j++) { Ny[j] = 0; Sy[j] = 0; Ny1[j] = 0; Sy1[j] = 0; zone_y[j] = 0; }
	
	istar = 0;	istop = 0;	jstar = 0;	jstop = 0;

	for (i = 0; i < IMAX; i++){
		if (fabs(bottombound_Y[i] - ystar) <= diff_y / 8.0) istar = i;
		if (fabs(bottombound_Y[i] - ystop) <= diff_y / 8.0){
			istop = i;
			break;
		}
	}
	for (j = 0; j < JMAX; j++){
		if (fabs(leftbound_Z[j] - zstar) <= diff_z / 8.0) jstar = j;
		if (fabs(leftbound_Z[j] - zstop) <= diff_z / 8.0){
			jstop = j;
			break;
		}
	}
	if (mpi_rank ==0 ) cout <<"istar = " << istar <<", istop = " << istop << ", jstar = " << jstar << ", jstop = " << jstop << endl;   

	// input the parameters of the bone
	speed_air1=2000;
	Aav1=1.0/1500.0;
	coef_tao1=-1.957321621e-9;
	coef_eta1=-6.200185788e-7;
	int Max_z = 0;
	int Max_y = 0;

	// Immersed Boundary searching
	switch (cr_judge){
	case 1: 
		for (i = 1; i < IMAX - 1; i++) {
			for (j = 1; j < JMAX - 1; j++) {
				cr_uc[i][j] = 0;
				if (bottombound_Y[i] >= 0.00) {
					//double cr_length=sqrt(pow((bottombound_Y[i]-0.03),2)+pow((leftbound_Z[j]-0.03),2));
					//if (cr_length<=0.015){
					cr_uc[i][j] = 1;
					if (Ny[j] == 0) Sy[j] = i;
					if (Nz[i] == 0) Sz[i] = j;
					Ny[j]++;
					Nz[i]++;
				}

				if (cr_uc[i - 1][j] == 0 && cr_uc[i][j] == 1) zone_y[j]++;
				if (cr_uc[i][j - 1] == 0 && cr_uc[i][j] == 1) zone_z[i]++;
			}
		}

		for (i = 1; i < IMAX - 1; i++) if (Nz[i] >= Max_z) Max_z = Nz[i];
		for (j = 1; j < JMAX - 1; j++) if (Ny[j] >= Max_y) Max_y = Ny[j];

		break;
	case 2:
		for (i = 1; i < IMAX - 1; i++) {
			for (j = 1; j < JMAX - 1; j++) {
				cr_uc[i][j] = 0;
				double cr_length = sqrt(pow((bottombound_Y[i] - 0.015), 2) + pow((leftbound_Z[j] - 0.01), 2));
				if (cr_length <= 0.005) {
					//double cr_length=sqrt(pow((bottombound_Y[i]-0.03),2)+pow((leftbound_Z[j]-0.03),2));
					//if (cr_length<=0.015){
					cr_uc[i][j] = 1;
					if (Ny[j] == 0) Sy[j] = i;
					if (Nz[i] == 0) Sz[i] = j;
					Ny[j]++;
					Nz[i]++;
				}
				if (cr_uc[i - 1][j] == 0 && cr_uc[i][j] == 1) zone_y[j]++;
				if (cr_uc[i][j - 1] == 0 && cr_uc[i][j] == 1) zone_z[i]++;
			}
		}

		for (i = 1; i < IMAX - 1; i++) if (Nz[i] >= Max_z) Max_z = Nz[i];
		for (j = 1; j < JMAX - 1; j++) if (Ny[j] >= Max_y) Max_y = Ny[j];

		break;
	case 3: // two strips
		for (i = 1; i < IMAX - 1; i++) {
			for (j = 1; j < JMAX - 1; j++) {
				cr_uc[i][j] = 0;
				if (bottombound_Y[i] >= 0.005 && bottombound_Y[i] <= 0.01 || bottombound_Y[i] >= 0.015 && bottombound_Y[i] <= 0.02) {
					//double cr_length=sqrt(pow((bottombound_Y[i]-0.03),2)+pow((leftbound_Z[j]-0.03),2));
					//if (cr_length<=0.015){
					cr_uc[i][j] = 1;
				}
				if (cr_uc[i - 1][j] == 0 && cr_uc[i][j] == 1) zone_y[j]++;
				if (cr_uc[i][j - 1] == 0 && cr_uc[i][j] == 1) zone_z[i]++;
				
				if (cr_uc[i][j] == 1) {
					if (zone_y[j] == 1 && Ny[j] == 0) Sy[j] = i;
					if (zone_y[j] == 1) Ny[j]++;
					if (zone_y[j] == 2 && Ny1[j] == 0) Sy1[j] = i;
					if (zone_y[j] == 2) Ny1[j]++;

					if (zone_z[i] == 1 && Nz[i] == 0) Sz[i] = j;
					if (zone_z[i] == 1) Nz[i]++;
					if (zone_z[i] == 2 && Nz1[i] == 0) Sz1[i] = j;
					if (zone_z[i] == 2) Nz1[i]++;
				}
			}
		}

		//for (i = 1; i < IMAX - 1; i++) if (Nz[i] >= Max_z) Max_z = Nz[i];
		//for (j = 1; j < JMAX - 1; j++) if (Ny[j] >= Max_y) Max_y = Ny[j];
		Max_z = JMAX;
		Max_y = IMAX;

		break; 
	case 4: // two strips
			for (i = 1; i < IMAX - 1; i++) {
				for (j = 1; j < JMAX - 1; j++) {
					cr_uc[i][j] = 0;
					double cr_length = sqrt(pow((bottombound_Y[i] - 0.015), 2) + pow((leftbound_Z[j] - 0.01), 2));
					if (cr_length <= 0.005 && cr_length >= 0.002) {
						cr_uc[i][j] = 1;
					}
					if (cr_uc[i - 1][j] == 0 && cr_uc[i][j] == 1) zone_y[j]++;
					if (cr_uc[i][j - 1] == 0 && cr_uc[i][j] == 1) zone_z[i]++;

					if (cr_uc[i][j] == 1) {
						if (zone_y[j] == 1 && Ny[j] == 0) Sy[j] = i;
						if (zone_y[j] == 1) Ny[j]++;
						if (zone_y[j] == 2 && Ny1[j] == 0) Sy1[j] = i;
						if (zone_y[j] == 2) Ny1[j]++;

						if (zone_z[i] == 1 && Nz[i] == 0) Sz[i] = j;
						if (zone_z[i] == 1) Nz[i]++;
						if (zone_z[i] == 2 && Nz1[i] == 0) Sz1[i] = j;
						if (zone_z[i] == 2) Nz1[i]++;
					}
				}
			}

			//for (i = 1; i < IMAX - 1; i++) if (Nz[i] >= Max_z) Max_z = Nz[i];
			//for (j = 1; j < JMAX - 1; j++) if (Ny[j] >= Max_y) Max_y = Ny[j];
			Max_z = JMAX;
			Max_y = IMAX;

			break;

	}
	

	//Get_geo();
	MPI_Barrier(MPI_COMM_WORLD);
	
	gy=new double *[mpi_IMAX[mpi_iindex]]; gz=new double *[mpi_IMAX[mpi_iindex]];
	wy=new double *[mpi_IMAX[mpi_iindex]]; wz=new double *[mpi_IMAX[mpi_iindex]];
	for (i=0;i<mpi_IMAX[mpi_iindex];i++){
		gy[i]=new double[mpi_JMAX[mpi_jindex]];	gz[i]=new double[mpi_JMAX[mpi_jindex]];
		wy[i]=new double[mpi_JMAX[mpi_jindex]];	wz[i]=new double[mpi_JMAX[mpi_jindex]];
	}
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			gy[i][j]=0.0; gz[i][j]=0.0;
			wy[i][j]=0.0; wz[i][j]=0.0;
		}
	}

	for (i=0;i<mpi_IMAX[mpi_iindex];i++){
		for (j = 0; j<mpi_JMAX[mpi_jindex];j++){
			int in=mpi_i1[mpi_iindex]+i;
			int jn=mpi_j1[mpi_jindex]+j;
			if (cr_uc[in][jn] != 0){
				gy[i][j] = cal_weight(coef_y,i,Sy[j]);
				gz[i][j] = cal_weight(coef_y,j,Sz[i]);
				wy[i][j] = cal_weight(coef_y - 1.0,i,Sy[j]);
				wz[i][j] = cal_weight(coef_y - 1.0,j,Sz[i]);				
			}
		}
	}
	//Get_g(0);
	gp1 = new double[Max_y + 1]; gp2 = new double[Max_y + 1];
	gq1 = new double[Max_z + 1]; gq2 = new double[Max_z + 1];
	for (i = 0; i < Max_y + 1; i++) { gp1[i] = 0.0; gp2[i] = 0.0; }
	for (j = 0; j < Max_z + 1; j++) { gq1[j] = 0.0; gq2[j] = 0.0; }

	cal_weight(coef_y, Max_y + 1, gp1);
	cal_weight(coef_y, Max_z + 1, gq1);
	cal_weight(coef_y - 1.0, Max_y + 1, gp2);
	cal_weight(coef_y - 1.0, Max_z + 1, gq2);

}

void air::Get_g(int Rank)
{
	int i, j;
	if (mpi_rank == Rank){
		char fileg[100] = "./Weight_check_rank", temp[10];
		int ff;
		ff = sprintf_s(temp,"_%d.dat",mpi_rank);
		strcat_s(fileg,temp);
		ofstream ofg(fileg, ios::out|ios::binary);
		ofg << "VARIABLES = \"Y\", \"Z\", \"gy\", \"gz\", \"wy\", \"wz\"" << endl;
		ofg << "ZONE T = \"Weight coefficients Zone\"" << endl;
		ofg << "I=" << mpi_IMAX[mpi_iindex] << ",J=" << mpi_JMAX[mpi_jindex] << ",F=POINT" << endl;				
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){				
			int in=mpi_i1[mpi_iindex]+i;
			int jn=mpi_j1[mpi_jindex]+j;
				//ofg << Y[i][j] << " " << Z[i][j] << " " << gy[i][j] << " " << gz[j] << " " << wy[i] << " " << wz[j] << endl; 
				ofg << bottombound_Y[in] << " " << leftbound_Z[jn] << " " << gy[i][j] << " " << gz[i][j] << " " << wy[i][j] << " " << wz[i][j] << endl;
			}
		}
		ofg.close();
	}
	if(mpi_rank == Rank) cout << "Weight output done!" << endl;	
}

void air::Get_geo(){
	int Rank = 0;
	if (mpi_rank == Rank){
		char filens[100] = "./Geo_check_rank", temp[10];
		int ff;
		ff = sprintf_s(temp,"_%d.dat",mpi_rank);
		strcat_s(filens,temp);
		ofstream ofNS(filens, ios::out|ios::binary);
		ofNS << "VARIABLES = \"Y\", \"Z\", \"C\", \"Ny\", \"Nz\", \"Sy\", \"Sz\"" << endl;
		ofNS << "ZONE T = \"Geo Zone\",";
		ofNS << "I=" << IMAX << ",J=" << JMAX << ",F=POINT" << endl;				
		for (int j=0;j<JMAX;j++){
			for(int i=0;i<IMAX;i++){
				ofNS << bottombound_Y[i] << " " << leftbound_Z[j] << " " << cr_uc[i][j] << " " << Ny[j] << " " << Nz[i] << " " << Sy[j] << " " << Sz[i] << endl; 
			}
		}
		ofNS.close();
	}
	if(mpi_rank == Rank) cout << "Geo output done!" << endl;
}

air::~air()
{
	for(int i=0;i<mpi_IMAX[mpi_iindex];i++){
		delete[] Vav[i];delete[] Wav[i];
		delete[] Vav_over_y[i];delete[] Vav_over_z[i];
		delete[] Wav_over_y[i];delete[] Wav_over_z[i];
		delete[] coef1[i];delete[] coef2[i];
		delete[] Qv[i];delete[] Qw[i];delete[] Qp[i];
		delete[] cr_uc[i]; 
		delete[] gy[i];delete[] gz[i];delete[] wy[i];delete[] wz[i];
	}
	for (int i=0;i<2;i++){
		delete[] coef_weighty[i];delete[] coef_weightz[i];
	}
	for (int i=0;i<IMAX;i++){
		delete[] whole_Qv[i];
		delete[] whole_Qw[i];
		delete[] whole_Qp[i];
	}
	delete[] whole_Qv;delete whole_Qw; delete whole_Qp;
	delete[] Vav;delete[] Wav;
	delete[] Vav_over_y;delete[] Vav_over_z;
	delete[] Wav_over_y;delete[] Wav_over_z;
	delete[] coef1;delete[] coef2;
	delete[] Qv; delete[] Qw;delete[] Qp;delete[] cr_uc;
	delete[] cr_y;delete[] cr_z; 
	delete[] coef_weighty;delete[] coef_weightz;
	delete[] Ny; delete[] Nz; delete[] Sy; delete[] Sz;
	delete[] gy; delete[] gz; delete[] wy; delete[] wz;
	delete SP;
}

void air::SetWindProfile(int N)
{
	int i,j;
	//calculate Z0 on local ground level
	for(j=0;j<mpi_JMAX[mpi_jindex];j++){
		for(i=0;i<mpi_IMAX[mpi_iindex];i++){
			double Z0;
			Z0=0;
			Vav[i][j]=SP->Vav(0,Y[i][j],Z[i][j]-Z0);
			Wav[i][j]=SP->Wav(0,Y[i][j],Z[i][j]-Z0);
		}
	}
	for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
		for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){

			coef1[i][j]=(Vav[i][j]*Z_over_z[i][j]-Wav[i][j]*Y_over_z[i][j])
			                                        *diff_t/diff_y/Jacobian[i][j];
			coef2[i][j]=(Wav[i][j]*Y_over_y[i][j]-Vav[i][j]*Z_over_y[i][j])
                                                    *diff_t/diff_z/Jacobian[i][j];
			Vav_over_y[i][j]=(Z_over_z[i][j]*(Vav[i+1][j]-Vav[i-1][j])/diff_y/2-
					Z_over_y[i][j]*(Vav[i][j+1]-Vav[i][j-1])/diff_z/2)/Jacobian[i][j];
			Wav_over_y[i][j]=(Z_over_z[i][j]*(Wav[i+1][j]-Wav[i-1][j])/diff_y/2-
					Z_over_y[i][j]*(Wav[i][j+1]-Wav[i][j-1])/diff_z/2)/Jacobian[i][j];
			Vav_over_z[i][j]=(Y_over_y[i][j]*(Vav[i][j+1]-Vav[i][j-1])/diff_z/2-
					Y_over_z[i][j]*(Vav[i+1][j]-Vav[i-1][j])/diff_y/2)/Jacobian[i][j];
			Wav_over_z[i][j]=(Y_over_y[i][j]*(Wav[i][j+1]-Wav[i][j-1])/diff_z/2-
					Y_over_z[i][j]*(Wav[i+1][j]-Wav[i-1][j])/diff_y/2)/Jacobian[i][j];

		}
	}
}
void air::save_restart_air(char *restartfile)
{
	int i,j;
	ofstream outfile11(restartfile,ios::app|ios::binary);
	outfile11.setf(ios::scientific,ios::floatfield);
	outfile11.setf(ios::showpos);
	outfile11.precision(6);
	outfile11.width(14);
	for (i=0;i<IMAX;i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			outfile11<<whole_Qv[i][j]<<" ";
			outfile11<<whole_Qw[i][j]<<" ";
			outfile11<<whole_Qp[i][j]<<" ";
		}
	}
	outfile11<<endl;
	outfile11.close();
}

void air::input_restart_air(char *restart_infile,int *point_pois)
{
	int i,j,int_pois;
	int_pois=*point_pois;
	ifstream infile(restart_infile,ios::in|ios::binary);
	infile.seekg(int_pois);
	for (i=0;i<IMAX;i++){
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			infile>>whole_Qv[i][j];
			infile>>whole_Qw[i][j];
			infile>>whole_Qp[i][j];
		}
	}
	infile.ignore(100,'\n');
	int_pois=infile.tellg();
	*point_pois=int_pois;
	infile.close();
}
//--------------------------calculate velocity of V and W for new time (nn)----------------//
void air::cal_fvw(double **temp_v,double **temp_w,
			double **temp_p)
{
	int i,in,j,jn;
	double absorb_coefz1, absorb_coefz2, temp_z, absorb_coefz3, absorb_coefz4, absorb_coefy, absorb_coefy1;
	//  PML boundary for y direction
	double temp_y;
	for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			in=mpi_i1[mpi_iindex]+i;
			jn=mpi_j1[mpi_jindex]+j;
			v_nn[i][j]=-Aav[j]*diff_t*(temp_p[i][j]-temp_p[i-1][j])/(Y[i][j]-Y[i-1][j]);
			w_nn[i][j]=-Aav[j]*diff_t*(temp_p[i][j]-temp_p[i][j-1])/(Z[i][j]-Z[i][j-1]);

			if (jn<JMAX1){		
				absorb_coefz1=PML_AbsorbZmax*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				absorb_coefz3=1.0+25.0*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				temp_z =Aav[j]*(Qp[i][j]-Qp[i-1][j])/(Y[i][j]-Y[i-1][j]);
				v_nn[i][j] -= absorb_coefz1*(temp_z+temp_v[i][j])*diff_t;
				w_nn[i][j] -= absorb_coefz1*temp_w[i][j]*diff_t+(1.0/absorb_coefz3-1.0)*
						Aav[j]*diff_t*(temp_p[i][j]-temp_p[i][j-1])/(Z[i][j]-Z[i][j-1]);

			}
			
			if(jn>=JMAX-JMAX2){
				absorb_coefz2=PML_AbsorbZmax*pow((Z[i][j]-PML_z2)/PML_Height2,PML_alpha);
				absorb_coefz4=1.0+25.0*pow((Z[i][j]-PML_z2)/PML_Height2,PML_alpha);

				temp_z =Aav[j]*(Qp[i][j]-Qp[i-1][j])/(Y[i][j]-Y[i-1][j]);
				v_nn[i][j] -= absorb_coefz2*(temp_z+temp_v[i][j])*diff_t;
							
				w_nn[i][j] -= absorb_coefz2*temp_w[i][j]*diff_t+(1.0/absorb_coefz4-1.0)*
						Aav[j]*diff_t*(temp_p[i][j]-temp_p[i][j-1])/(Z[i][j]-Z[i][j-1]);
			}
			//if (in < IMAX2 && IMAX1 == 26000) { v_nn[i][j] = 0.0; w_nn[i][j] = 0.0 }
			if (in<IMAX2 && IMAX1 >= 26000){

				absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
				absorb_coefz1 = PML_AbsorbZmax*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);

				temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
				v_nn[i][j] -= absorb_coefy * temp_y * diff_t + (1.0 / absorb_coefy1 - 1.0)*	Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (temp_p[i][j] - temp_p[i - 1][j]) -
					Z_over_y[i][j] / diff_z * (temp_p[i][j] - temp_p[i][j - 1])) * diff_t;

				temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
					Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) + temp_w[i][j];

				if (jn >= JMAX - JMAX2 && jn<JMAX){
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz2 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz2 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz2*Qw[i][j];
				}
				else if (jn<JMAX1){
					absorb_coefz4 = 1.0 + 25.0*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz1 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz1 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz1*Qw[i][j];
				}
				w_nn[i][j] -= absorb_coefy*temp_y*diff_t;
			}
			if (in >= IMAX - IMAX2){// && jn >= JMAX1

				absorb_coefy = PML_AbsorbYmax*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefz1 = PML_AbsorbZmax*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);

				temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
				v_nn[i][j] -= absorb_coefy*temp_y*diff_t
					+ (1.0 / absorb_coefy1 - 1.0)*
					Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] * diff_t / diff_y*(temp_p[i][j] - temp_p[i - 1][j]) -
					Z_over_y[i][j] * diff_t / diff_z*(temp_p[i][j] - temp_p[i][j - 1]));

				temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
					Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) + temp_w[i][j];

				if (jn >= JMAX - JMAX2 && jn<JMAX){
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz2 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz2 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz2*Qw[i][j];
				}
				else if (jn<JMAX1){
					absorb_coefz4 = 1.0 + 25.0*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y2 - Y[i][j]) / PML_Width2, PML_alphay);
					temp_y = Qw[i][j] * Vav_over_z[i][j] + temp_v[i][j];
					v_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * absorb_coefz1 * Aav[j] / Jacobian[i][j] * (Z_over_z[i][j] / diff_y * (Qp[i][j] - Qp[i - 1][j]) -
						Z_over_y[i][j] / diff_z * (Qp[i][j] - Qp[i][j - 1])) * diff_t + absorb_coefy * absorb_coefz1 * Qv[i][j] * diff_t;//
					temp_y = Aav[j] / Jacobian[i][j] * (Y_over_y[i][j] / diff_z*(Qp[i][j] - Qp[i][j - 1]) -
						Y_over_z[i][j] / diff_y*(Qp[i][j] - Qp[i - 1][j])) / absorb_coefz4 + temp_w[i][j] + absorb_coefz1*Qw[i][j];
				}
				w_nn[i][j] -= absorb_coefy*temp_y*diff_t;
			}
		}
	}
	if(cr_judge != 0){
		for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
			for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){	
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;		
				if(cr_uc[in][jn] == 1){
					v_nn[i][j]=-Aav1*diff_t*(temp_p[i][j]-temp_p[i-1][j])/(Y[i][j]-Y[i-1][j]);
					w_nn[i][j]=-Aav1*diff_t*(temp_p[i][j]-temp_p[i][j-1])/(Z[i][j]-Z[i][j-1]);
				}	
			}
		}
	}
}
//------------------------------calculate pressure (p) for new time(nn)------------//
void air::cal_fp(double **temp_v,double **temp_w,
			double **temp_p)
{
	int i,in,j,jn;
	double absorb_coefz1, absorb_coefz2, temp_z, absorb_coefz3, absorb_coefz4, absorb_coefy, absorb_coefy1;
	//  PML boundary for y direction
	double temp_y;
// specify the initial condition for circle

	IMAX1=IMAX1+1;
	
	for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			in=mpi_i1[mpi_iindex]+i;
			jn=mpi_j1[mpi_jindex]+j;
			p_nn[i][j]=-speed_air[j]*speed_air[j]/Aav[j]*diff_t*
					((temp_v[i+1][j]-temp_v[i][j])/(Y[i+1][j]-Y[i][j])+
					(temp_w[i][j+1]-temp_w[i][j])/(Z[i][j+1]-Z[i][j]))+
					speed_air[j]*speed_air[j]*diff_t*coef_tao*swp_nn[in][jn]+
					speed_air[j]*speed_air[j]/Aav[j]*diff_t*coef_eta*swvw_nn[in][jn];
			//cr_length=sqrt(pow((Y[i][j]-0.03),2)+pow((Z[i][j]-0.03),2));
			//if (cr_length>=0.005 && cr_length<=0.015)


			if (jn<JMAX1){
				absorb_coefz1=PML_AbsorbZmax*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				absorb_coefz3=1.0+25.0*pow((PML_z1-Z[i][j])/PML_Height1,PML_alpha);
				temp_z=speed_air[j]*speed_air[j]/Aav[j]*(Qv[i+1][j]-Qv[i][j])/(Y[i+1][j]-Y[i][j]);
				p_nn[i][j] -= absorb_coefz1*(temp_z+temp_p[i][j])*diff_t+(1.0/absorb_coefz3-1.0)*
						speed_air[j]*speed_air[j]/Aav[j]*diff_t*(temp_w[i][j+1]-temp_w[i][j])/(Z[i][j+1]-Z[i][j]);
			}
			if (jn>=JMAX-JMAX2){
				absorb_coefz2=PML_AbsorbZmax*pow((Z[i][j]-PML_z2)/PML_Height2,PML_alpha);
				absorb_coefz4=1.0+25.0*pow((Z[i][j]-PML_z2)/PML_Height2,PML_alpha);
				temp_z=speed_air[j]*speed_air[j]/Aav[j]*(Qv[i+1][j]-Qv[i][j])/(Y[i+1][j]-Y[i][j]);
				p_nn[i][j] -= absorb_coefz2*(temp_z+temp_p[i][j])*diff_t+(1.0/absorb_coefz4-1.0)*
						speed_air[j]*speed_air[j]/Aav[j]*diff_t*(temp_w[i][j+1]-temp_w[i][j])/(Z[i][j+1]-Z[i][j]);
			}

			//if (in < IMAX2 && IMAX1 == 26000) { p_nn[i][j] = 0.0 }
			if (in<IMAX2 && IMAX1 >= 26000){

				absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width2, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((PML_y1 - Y[i][j]) / PML_Width2, PML_alphay);

				temp_y = 1.0 / Jacobian[i][j] * speed_air[j] * speed_air[j] / Aav[j] * (Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
					Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) + temp_p[i][j];

				if (jn >= JMAX - JMAX2 && jn<JMAX)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width2, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * speed_air[j] * speed_air[j] / Aav[j] * (Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];

				}
				else if (jn<JMAX1)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z1) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y1 - Y[i][j]) / PML_Width1, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * speed_air[j] * speed_air[j] / Aav[j] * (Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];
				}
				p_nn[i][j] -= absorb_coefy*temp_y*diff_t
					+ (1.0 / absorb_coefy1 - 1.0)*(Z_over_z[i][j] * diff_t / diff_y*(temp_v[i + 1][j] - temp_v[i][j]) - Z_over_y[i][j] * diff_t / diff_z*(temp_v[i][j + 1] - temp_v[i][j]))*speed_air[j] * speed_air[j] / Aav[j] / Jacobian[i][j];
			}
			if (in >= IMAX - IMAX2){// && jn>=JMAX1

				absorb_coefy = PML_AbsorbYmax*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefy1 = 1.0 + 25.0*pow((Y[i][j] - PML_y2) / PML_Width2, PML_alphay);
				absorb_coefz1 = PML_AbsorbZmax*pow((PML_z1 - Z[i][j]) / PML_Height1, PML_alpha);
				absorb_coefz2 = PML_AbsorbZmax*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);

				temp_y = 1.0 / Jacobian[i][j] * speed_air[j] * speed_air[j] / Aav[j] * (Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
					Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) + temp_p[i][j];

				if (jn >= JMAX - JMAX2 && jn<JMAX)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z2) / PML_Height2, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y2 - Y[i][j]) / PML_Width2, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * speed_air[j] * speed_air[j] / Aav[j] * (Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];
				}
				else if (jn<JMAX1)
				{
					absorb_coefz4 = 1.0 + 25.0*pow((Z[i][j] - PML_z1) / PML_Height1, PML_alpha);
					absorb_coefy = PML_AbsorbYmax*pow((PML_y2 - Y[i][j]) / PML_Width2, PML_alphay);
					temp_y = 1.0 / Jacobian[i][j] * speed_air[j] * speed_air[j] / Aav[j] * (Y_over_y[i][j] / diff_z*(Qw[i][j + 1] - Qw[i][j]) -
						Y_over_z[i][j] / diff_y*(Qw[i + 1][j] - Qw[i][j])) / absorb_coefz4 + temp_p[i][j];
				}
				p_nn[i][j] -= absorb_coefy*temp_y*diff_t
					+ (1.0 / absorb_coefy1 - 1.0)*(Z_over_z[i][j] * diff_t / diff_y*(temp_v[i + 1][j] - temp_v[i][j]) - Z_over_y[i][j] * diff_t / diff_z*(temp_v[i][j + 1] - temp_v[i][j]))*speed_air[j] * speed_air[j] / Aav[j] / Jacobian[i][j];
				if (jn >= JMAX - JMAX2 && jn < JMAX)
				{
					p_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * speed_air[j] * speed_air[j] / Aav[j] * absorb_coefz2 / Jacobian[i][j] * (Z_over_z[i][j] / diff_y*(Qv[i + 1][j] - Qv[i][j])
						- Z_over_y[i][j] / diff_z*(Qv[i][j + 1] - Qv[i][j])) * diff_t + absorb_coefy*absorb_coefz2*Qp[i][j] * diff_t;
				}
				else if (jn < JMAX1)
				{
					p_nn[i][j] -= (1.0 / absorb_coefy1 - 1.0) * speed_air[j] * speed_air[j] / Aav[j] * absorb_coefz1 / Jacobian[i][j] * (Z_over_z[i][j] / diff_y*(Qv[i + 1][j] - Qv[i][j])
						- Z_over_y[i][j] / diff_z*(Qv[i][j + 1] - Qv[i][j])) * diff_t + absorb_coefy*absorb_coefz1*Qp[i][j] * diff_t;
				}/**/
			}
		}
	}
	if(cr_judge != 0){
		for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
			for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){	
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if(cr_uc[in][jn] == 1){
					p_nn[i][j]=-speed_air1*speed_air1/Aav1*diff_t*
					((temp_v[i+1][j]-temp_v[i][j])/(Y[i+1][j]-Y[i][j])+
					(temp_w[i][j+1]-temp_w[i][j])/(Z[i][j+1]-Z[i][j]))+
					speed_air1*speed_air1*diff_t*coef_tao1*swp_nn[in][jn]+
					speed_air1*speed_air1/Aav1*diff_t*coef_eta1*swvw_nn[in][jn];
				}
			}
		}
	}
}

void air::cal_delvw()
{
	int i,j;
	for (i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			delvw_nn[i][j]=(v_nn[i+1][j]-v_nn[i][j])/(Y[i+1][j]-Y[i][j])+
							(w_nn[i][j+1]-w_nn[i][j])/(Z[i][j+1]-Z[i][j]);
		}
	}
}
void air::fractional_CD()
{
	int i, j, p, q, Xl, Xr, Yt, Yb, P, Q;
	double temp_x1, temp_z1, temp_x2, temp_z2;

	for (i = 1; i < IMAX - 1; i++) {
		for (j = 1; j < JMAX - 1; j++) {
			if (cr_uc[i][j] != 0) {
				// X(Y) direction:
				temp_x1 = 0.0;	temp_x2 = 0.0;
				if (zone_y[j] == 1 || (zone_y[j] == 2 && Y[i][j] >= Sy[j] && Y[i][j] <= Sy[j] + Ny[j])) { // first piece
					Xl = Sy[j]; Xr = Xl + Ny[j];
					for (p = (i - Xr); p <= (i - Xl); p++) {
						int jn = j - mpi_j1[mpi_jindex];
						P = abs(p);
						if (jn >= 1 && jn < mpi_JMAX[mpi_jindex] - 1) {
							temp_x1 += gp1[P] * p_n[i - p][jn];
							temp_x2 += gp2[P] * delvw_nn[i - p][jn];
						}
						else {
							temp_x1 = 0;
							temp_x2 = 0;
						}
					}
				}
				else if (zone_y[j] == 2 && Y[i][j] >= Sy1[j] && Y[i][j] <= Sy1[j] + Ny1[j]) { // second piece
					Xl = Sy1[j]; Xr = Xl + Ny1[j];
					for (p = (i - Xr); p <= (i - Xl); p++) {
						int jn = j - mpi_j1[mpi_jindex];
						P = abs(p);
						if (jn >= 1 && jn < mpi_JMAX[mpi_jindex] - 1) {
							temp_x1 += gp1[P] * p_n[i - p][jn];
							temp_x2 += gp2[P] * delvw_nn[i - p][jn];
						}
						else {
							temp_x1 = 0;
							temp_x2 = 0;
						}
					}
				}
				// Z direction
				temp_z1 = 0.0;  temp_z2 = 0.0;
				if (zone_z[i] == 1 || (zone_z[i] == 2 && Z[i][j] >= Sz[i] && Z[i][j] <= Sz[i] + Nz[i])) {
					Yb = Sz[i]; Yt = Yb + Nz[i];
					for (q = (j - Yt); q <= (j - Yb); q++) {
						int jn = j - q - mpi_j1[mpi_jindex];
						Q = abs(q);
						if (jn >= 1 && jn < mpi_JMAX[mpi_jindex] - 1) {
							temp_z1 += gq1[Q] * p_n[i][jn];
							temp_z2 += gq2[Q] * delvw_nn[i][jn];
						}
						else {
							temp_z1 = 0;
							temp_z2 = 0;
						}
					}
				}
				else if (zone_z[i] == 2 && Z[i][j] >= Sz1[i] && Z[i][j] <= Sz1[i] + Nz1[i]) {
					Yb = Sz1[i]; Yt = Yb + Nz1[i];
					for (q = (j - Yt); q <= (j - Yb); q++) {
						int jn = j - q - mpi_j1[mpi_jindex];
						Q = abs(q);
						if (jn >= 1 && jn < mpi_JMAX[mpi_jindex] - 1) {
							temp_z1 += gq1[Q] * p_n[i][jn];
							temp_z2 += gq2[Q] * delvw_nn[i][jn];
						}
						else {
							temp_z1 = 0;
							temp_z2 = 0;
						}
					}
				}
				swp_nn[i][j] = (pow_y1 * temp_x1 + pow_z1 * temp_z1);
				swvw_nn[i][j] = (pow_y2 * temp_x2 + pow_z2 * temp_z2);
			}
			else {
				swp_nn[i][j] = 0.0;
				swvw_nn[i][j] = 0.0;
			}
		}
	}
}
double air::Cal_g(int p, double beta)
{
	double w;
	w = pow(-1, p) * tgamma(beta + 1) / (tgamma(0.5 * beta - p + 1) * tgamma(0.5 * beta + p + 1));
	return w;
}
void air::shift_grunwald()
{
	int i, j, in, jn, lz, ly, jy, jz;
	double temp_y, temp_z, temp_y1, temp_z1;

	for (i = 1; i < mpi_IMAX[mpi_iindex] - 1; i++) {
		for (j = 1; j < mpi_JMAX[mpi_jindex] - 1; j++) {
			in = mpi_i1[mpi_iindex] + i;
			jn = mpi_j1[mpi_jindex] + j;
			if (cr_uc[in][jn] != 0) {
				temp_y = 0.0;	temp_z = 0.0;
				temp_y1 = 0.0; temp_z1 = 0.0;
				ly = i - Sy[j];
				lz = j - Sz[i];

				for (jy = 0; jy <= ly; jy++) { temp_y1 += gy[i][j] * delvw_nn[i - jy][j]; }
				for (jy = 0; jy < Ny[j] - ly; jy++) { temp_y1 += gy[i][j] * delvw_nn[i + jy][j]; }

				for (jz = 0; jz <= lz; jz++) { temp_z1 += gz[i][j] * delvw_nn[i][j - jz]; }
				for (jz = 0; jz < Nz[i] - lz; jz++) { temp_z1 += gz[i][j] * delvw_nn[i][j + jz]; }

				for (jy = 0; jy <= ly + 1; jy++) { temp_y += wy[i][j] * p_n[i + 1 - jy][j]; }
				for (jy = 0; jy < Ny[j] - ly + 1; jy++) { temp_y += wy[i][j] * p_n[i - 1 + jy][j]; }

				for (jz = 0; jz <= lz + 1; jz++) { temp_z += wz[i][j] * p_n[i][j + 1 - jz]; }
				for (jz = 0; jz < Nz[i] - lz + 1; jz++) { temp_z += wz[i][j] * p_n[i][j - 1 + jz]; }

				swp_nn[i][j] = (temp_y*pow_y1 + temp_z*pow_z1)*var_miu1;
				swvw_nn[i][j] = (temp_y1*pow_y2 + temp_z1*pow_z2)*var_miu2;
			}
			else {
				swp_nn[i][j] = 0.0;
				swvw_nn[i][j] = 0.0;
			}
		}
	}
}

	/*for (i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			temp_y=0.0;
			temp_z=0.0;
			temp_y1=0.0;
			temp_z1=0.0;
			for (n=0;n<i+2;n++){
				temp_y=temp_y+coef_weighty[0][n]*p_n[i-n+1][j];
			if (n<i+1)temp_y1=temp_y1+coef_weighty[1][n]*delvw_nn[i-n][j];
			}
			for (n=0;n<mpi_IMAX[mpi_iindex]-i+1;n++){
				temp_y=temp_y+coef_weighty[0][n]*p_n[i+n-1][j];
			if(n<mpi_IMAX[mpi_iindex]-i)temp_y1=temp_y1+coef_weighty[1][n]*delvw_nn[i+n][j];
			}
			for (n=0;n<j+2;n++){
				temp_z=temp_z+coef_weightz[0][n]*p_n[i][j-n+1];
			if (n<mpi_JMAX[mpi_jindex]+1)temp_z1=temp_z1+coef_weightz[1][n]*delvw_nn[i][j-n];
			}
			for (n=0;n<mpi_JMAX[mpi_jindex]-j+1;n++){
				temp_z=temp_z+coef_weightz[0][n]*p_n[i][j+n-1];
			if (n<mpi_JMAX[mpi_jindex]-j)temp_z1=temp_z1+coef_weightz[1][n]*delvw_nn[i][j+n];
			}
			swp_nn[i][j]=(temp_y*pow_y1+temp_z*pow_z1)*var_miu1;
			swvw_nn[i][j]=(temp_y1*pow_y2+temp_z1*pow_z2)*var_miu2;
		}
	}*/
void air::cal_weight(double y, int K, double *g) {
	double w;
	w = tgamma(y + 1) / (tgamma(0.5 * y + 1) * tgamma(0.5 * y + 1)); // w0
	for (int k = 0; k < K; k++) {
		g[k] = w;
		w = (1 - (y + 1) / (0.5 * y + k + 1)) * w;
	}
}
double air::cal_weight(double miu,int n)
{
	double temp_weight;
	int i;
	for (i=1; i<=n; i++){
		if (i==1)temp_weight=miu;
		else temp_weight=temp_weight*(miu-i+1)/double(i);
	}
	return pow(-1.0,n)*temp_weight;
}

double air::cal_weight(double y, int l, int s)
{
	double temp_weight;
	int i, j;
	j = l - s;
	temp_weight = 1;
	if (j == 0) return 1;
	else{
		for (i = 1; i <= j; i++){
			temp_weight = temp_weight*(y-i+1)/double(i);			
		}
		return pow(-1.0,j)*temp_weight;		
	}
}

//-----------------------------update boundary condtions of pressure--------------------//
void air::UpdateBC_pressure(boundary_location BC,int time_judge,int time_current)
{
	int i,j,in,jn,rank_send,rank_recv,tag;
	double time_length;
	// specify sender matrix in each block's boundary 
	for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		if (mpi_jindex!=0) pss[i]=p_nn[i][1];
		if (mpi_jindex!=mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-2;
			pns[i]=p_nn[i][jn];
		}
	}
	for (j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
		if (mpi_iindex!=0)pws[j]=p_nn[1][j];
		if (mpi_iindex!=mpi_yarea-1){
			in=mpi_IMAX[mpi_iindex]-2;
			pes[j]=p_nn[in][j];
		}	
	}
	switch(BC){
	case WestBC:
	{
		if (mpi_iindex==0){ 
			for(j=0;j<mpi_JMAX[mpi_jindex];j++) {
				p_nn[0][j]=p_nn[1][j];////left
			}
			// time series boundary condition 
			if (jjmax2!=0 && time_judge==1){
				time_length=time_current*diff_t;
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					jn=mpi_j1[mpi_jindex]+j;
					if (jn>=jjmax1 && jn<jjmax2){
						p_nn[0][j]=cos(8.16814076e6*time_length)*(1.0-cos(1.256637e6*time_length));
							}
						}
					}

		}else{
			tag=1;
			if (mpi_porous!=0)rank_recv=(mpi_iindex-1)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex-1)*mpi_zarea+mpi_jindex;
			MPI_Send(pws, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=2;
			rank_send=rank_recv;
			MPI_Recv(pwr, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		}

	}
	break;
	case EastBC:
	{
		if (mpi_iindex==mpi_yarea-1){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){ 
				in=mpi_IMAX[mpi_iindex]-1;
				p_nn[in][j]=p_nn[in-1][j];//right
			}
		}else{
			tag=1;
			if (mpi_porous!=0)rank_send=(mpi_iindex+1)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+1)*mpi_zarea+mpi_jindex;
			MPI_Recv(per, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=2;
			rank_recv=rank_send;
			MPI_Send(pes, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
	break;
	case SouthBC:
	{
		if (mpi_jindex==0){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++)p_nn[i][0]=p_nn[i][1]; //Bottom
		}else{
			tag=3;
			if (mpi_porous!=0)rank_recv=mpi_iindex*mpi_porous+mpi_jindex-1;
			else rank_recv=mpi_iindex*mpi_zarea+mpi_jindex-1;
			MPI_Send(pss, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=4;
			rank_send=rank_recv;
			MPI_Recv(psr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		}
	}
	break;
	case NorthBC://characteristic boundary condition
	{
		if (mpi_jindex==mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++) p_nn[i][jn]=p_nn[i][jn-1];//upper
		}else{
			tag=3;
			if (mpi_porous!=0)rank_send=mpi_iindex*mpi_porous+mpi_jindex+1;
			else rank_send=mpi_iindex*mpi_zarea+mpi_jindex+1;
			MPI_Recv(pnr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=4;
			rank_recv=rank_send;
			MPI_Send(pns, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
	break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
	
	
	for (i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		if(mpi_jindex!=0) p_nn[i][0]=psr[i];
		if(mpi_jindex!=mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			p_nn[i][jn]=pnr[i];
		}
	}
	for (j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
		if (mpi_iindex!=0)p_nn[0][j]=pwr[j];
		if (mpi_iindex!=mpi_yarea-1){
			in=mpi_IMAX[mpi_iindex]-1;
			p_nn[in][j]=per[j];
		}	
	}
}

//--------------------------update boundary conditions of velocity(V and W)------------------//
void air::UpdateBC_velocity(boundary_location BC)
{
	//v_nn[1][j] and w_nn[i][1] for Salomons scheme is real boundary value, 
	//each time it should be reset.
		// specify sender matrix in each block's boundary 
	int i,j,in,jn,rank_send,rank_recv,tag;
	for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		if (mpi_jindex!=0) {
			vss[i]=v_nn[i][1];
			wss[i]=w_nn[i][1];
		}
		if (mpi_jindex!=mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-2;
			vns[i]=v_nn[i][jn];
			wns[i]=w_nn[i][jn];
		}
	}
	for (j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
		if (mpi_iindex!=0){
			vws[j]=v_nn[1][j];
			wws[j]=w_nn[1][j];
		}
		if (mpi_iindex!=mpi_yarea-1){
			in=mpi_IMAX[mpi_iindex]-2;
			ves[j]=v_nn[in][j];
			wes[j]=w_nn[in][j];
		}
	
	}
	switch(BC){
	case WestBC:
	{
		if (mpi_iindex==0){ 
			for(j=0;j<mpi_JMAX[mpi_jindex];j++) {
				v_nn[0][j]=0.0;w_nn[0][j]=w_nn[1][j];////left
			}
		}else{
			tag=5;
			if (mpi_porous!=0)rank_recv=(mpi_iindex-1)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex-1)*mpi_zarea+mpi_jindex;
			MPI_Send(vws, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=6;
			MPI_Send(wws, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=7;
			rank_send=rank_recv;
			MPI_Recv(vwr, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=8;
			MPI_Recv(wwr, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		}
	}
		break;
	case EastBC:
	{
		if (mpi_iindex==mpi_yarea-1){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){ 
				in=mpi_IMAX[mpi_iindex]-1;
				v_nn[in][j]=0.0;//right
				w_nn[in][j]=w_nn[in-1][j];
			}
		}else{
			tag=5;
			if (mpi_porous!=0)rank_send=(mpi_iindex+1)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+1)*mpi_zarea+mpi_jindex;
			MPI_Recv(ver, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=6;
			MPI_Recv(wer, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=7;
			rank_recv=rank_send;
			MPI_Send(ves, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=8;
			MPI_Send(wes, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
		break;
	case SouthBC:
	{
		if (mpi_jindex==0){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				w_nn[i][0]=0.0;v_nn[i][0]=v_nn[i][1]; //Bottom
			}
		}else{
			tag=9;
			if (mpi_porous!=0)rank_recv=mpi_iindex*mpi_porous+mpi_jindex-1;
			else rank_recv=mpi_iindex*mpi_zarea+mpi_jindex-1;
			MPI_Send(vss, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=10;
			MPI_Send(wss, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=11;
			rank_send=rank_recv;
			MPI_Recv(vsr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=12;
			MPI_Recv(wsr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		}
	}
		break;
	case NorthBC://characteristic boundary condition
	{
		if (mpi_jindex==mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				w_nn[i][jn]=0.0;v_nn[i][jn]=v_nn[i][jn-1];//upper
			}
		}else{
			tag=9;
			if (mpi_porous!=0)rank_send=mpi_iindex*mpi_porous+mpi_jindex+1;
			else rank_send=mpi_iindex*mpi_zarea+mpi_jindex+1;
			MPI_Recv(vnr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=10;
			MPI_Recv(wnr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=11;
			rank_recv=rank_send;
			MPI_Send(vns, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=12;
			MPI_Send(wns, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
	break;
	default:
		cout<<"it is wrong to entrance boundary";
	}

	for (i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		if(mpi_jindex!=0) {
			v_nn[i][0]=vsr[i];
			w_nn[i][0]=wsr[i];
		}
		if(mpi_jindex!=mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			v_nn[i][jn]=vnr[i];
			w_nn[i][jn]=wnr[i];
		}
	}
	for (j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
		if (mpi_iindex!=0){
			v_nn[0][j]=vwr[j];
			w_nn[0][j]=wwr[j];
		}
		if (mpi_iindex!=mpi_yarea-1){
			in=mpi_IMAX[mpi_iindex]-1;
			v_nn[in][j]=ver[j];
			w_nn[in][j]=wer[j];
		}	
	}
}
void air::Update_PML_Qvw()
{
	int i,j,in,jn,jn1,i1,i2;
	int mpi_index,rank_num,rank_send,rank_recv,tag;

	for (i = 0; i<mpi_IMAX[mpi_iindex]; i++){
		for (j = 0; j < mpi_JMAX[mpi_jindex]; j++){
			in = mpi_i1[mpi_iindex] + i;
			jn = mpi_j1[mpi_jindex] + j;
			if (jn < JMAX1 || jn >= JMAX - JMAX2){//
				Qv[i][j] += diff_t*v_nn[i][j];
				Qw[i][j] += diff_t*w_nn[i][j];
			}
			if (in < IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000) {//
				Qv[i][j] += diff_t*v_nn[i][j];
				Qw[i][j] += diff_t*w_nn[i][j];
			}
			if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1){
				Qv[i][j] += diff_t*v_nn[i][j];
				Qw[i][j] += diff_t*w_nn[i][j];
			}
		}
	}
	if (mpi_iindex!=0){
		rank_num=mpi_IMAX[mpi_iindex]*mpi_JMAX[mpi_jindex];
		rank_recv=mpi_jindex;
		if (mpi_iindex!=0)i1=1;
		else i1=0;
		if (mpi_iindex!=mpi_yarea-1)i2=1;
		else i2=0;
		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					mpi_fs[mpi_index]=Qv[i][j];
					mpi_index=mpi_index+1;
				}
				if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					mpi_fs[mpi_index] = Qv[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					mpi_fs[mpi_index] = Qv[i][j];
					mpi_index = mpi_index + 1;
				}
			}
		}
		tag=24;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		
		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					mpi_fs[mpi_index]=Qw[i][j];
					mpi_index=mpi_index+1;
				}
				if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					mpi_fs[mpi_index] = Qw[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					mpi_fs[mpi_index] = Qw[i][j];
					mpi_index = mpi_index + 1;
				}
			}
		}
		tag=25;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		
		rank_send=rank_recv;
		rank_num=IMAX*mpi_JMAX[mpi_jindex];
		tag=26;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					whole_Qv[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
				if (i <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					whole_Qv[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					whole_Qv[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
			}
		}

		tag=27;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					whole_Qw[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
				if (i <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					whole_Qw[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					whole_Qw[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
			}
		}
	}else{
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					whole_Qv[in][j]=Qv[i][j];
					whole_Qw[in][j]=Qw[i][j];
				}
				if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					whole_Qv[in][j] = Qv[i][j];
					whole_Qw[in][j] = Qw[i][j];
				}
				if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					whole_Qv[in][j] = Qv[i][j];
					whole_Qw[in][j] = Qw[i][j];
				}
			}
		}
		for (in=1;in<mpi_yarea;in++){
			if (mpi_porous!=0)rank_send=(mpi_iindex+in)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+in)*mpi_zarea+mpi_jindex;
			jn=mpi_jindex;
			rank_num=mpi_IMAX[in]*mpi_JMAX[jn];
			if (in!=0)i1=1;
			else i1=0;
			if (in!=mpi_yarea-1)i2=1;
			else i2=0;
			mpi_index=0;
			tag=24;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					jn1=mpi_j1[mpi_jindex]+j;
					if (jn1<JMAX1||jn1>=JMAX-JMAX2){
						whole_Qv[i][j]=mpi_fr[mpi_index];
						mpi_index=mpi_index+1;
					}
					if (i <= IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 26000){
						whole_Qv[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 0){
						whole_Qv[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
				}
			}
			mpi_index=0;
			tag=25;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					jn1=mpi_j1[mpi_jindex]+j;
					if (jn1<JMAX1||jn1>=JMAX-JMAX2){
						whole_Qw[i][j]=mpi_fr[mpi_index];
						mpi_index=mpi_index+1;
					}
					if (i <= IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 26000){
						whole_Qw[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 0){
						whole_Qw[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
				}
			}
		}
		for (in=1;in<mpi_yarea;in++){
			jn=mpi_jindex;
			rank_num=IMAX*mpi_JMAX[jn];
			if (mpi_porous!=0) rank_recv=(mpi_iindex+in)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex+in)*mpi_zarea+mpi_jindex;
		
			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					jn1=mpi_j1[jn]+j;
					if (jn1<JMAX1||jn1>=JMAX-JMAX2){
						mpi_fs[mpi_index]=whole_Qv[i][j];
						mpi_index=mpi_index+1;
					}
					if (i <= IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 26000){
						mpi_fs[mpi_index] = whole_Qv[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 0){
						mpi_fs[mpi_index] = whole_Qv[i][j];
						mpi_index = mpi_index + 1;
					}
				}
			}
			tag=26;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			
			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					jn1=mpi_j1[mpi_jindex]+j;
					if (jn1<JMAX1||jn1>=JMAX-JMAX2){
						mpi_fs[mpi_index]=whole_Qw[i][j];
						mpi_index=mpi_index+1;
					}
					if (i <= IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 26000){
						mpi_fs[mpi_index] = whole_Qw[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 0){
						mpi_fs[mpi_index] = whole_Qw[i][j];
						mpi_index = mpi_index + 1;
					}
				}
			}
			tag=27;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
}

void air::Update_PML_Qp()
{
	int i,j,in,jn,jn1,i1,i2;
	int mpi_index,rank_num,rank_send,rank_recv,tag;
	for (i = 0; i < mpi_IMAX[mpi_iindex]; i++){
		for (j = 0; j < mpi_JMAX[mpi_jindex]; j++){
			in = mpi_i1[mpi_iindex] + i;
			jn = mpi_j1[mpi_jindex] + j;
			if (jn < JMAX1)Qp[i][j] += diff_t*p_nn[i][j];//
			if (jn >= JMAX - JMAX2) Qp[i][j] += diff_t*p_nn[i][j];
			if (in < IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
				Qp[i][j] += diff_t*p_nn[i][j];
			}
			if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1){ 
				Qp[i][j] += diff_t*p_nn[i][j];
			}
		}
	}
	if (mpi_iindex!=0){
		rank_num=mpi_IMAX[mpi_iindex]*mpi_JMAX[mpi_jindex];
		rank_recv=mpi_jindex;
		if (mpi_iindex!=0)i1=1;
		else i1=0;
		if (mpi_iindex!=mpi_yarea-1)i2=1;
		else i2=0;
		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					mpi_fs[mpi_index]=Qp[i][j];
					mpi_index=mpi_index+1;
				}
				if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					mpi_fs[mpi_index] = Qp[i][j];
					mpi_index = mpi_index + 1;
				}
				if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					mpi_fs[mpi_index] = Qp[i][j];
					mpi_index = mpi_index + 1;
				}
			}
		}
		tag=28;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

		rank_send=rank_recv;
		rank_num=IMAX*mpi_JMAX[mpi_jindex];
		tag=29;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					whole_Qp[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
				if (i <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					whole_Qp[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
				if (i >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					whole_Qp[i][j] = mpi_fr[mpi_index];
					mpi_index = mpi_index + 1;
				}
			}
		}
	}else{
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					whole_Qp[in][j]=Qp[i][j];
				}
				if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){
					whole_Qp[in][j] = Qp[i][j];
				}
				if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					whole_Qp[in][j] = Qp[i][j];
				}
			}
		}
		for (in=1;in<mpi_yarea;in++){
			if (mpi_porous!=0) rank_send=(mpi_iindex+in)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+in)*mpi_zarea+mpi_jindex;
			
			jn=mpi_jindex;
			rank_num=mpi_IMAX[in]*mpi_JMAX[jn];
			if (in!=0)i1=1;
			else i1=0;
			if (in!=mpi_yarea-1)i2=1;
			else i2=0;
			mpi_index=0;
			tag=28;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					jn1=mpi_j1[mpi_jindex]+j;
					if (jn1<JMAX1||jn1>=JMAX-JMAX2){
						whole_Qp[i][j]=mpi_fr[mpi_index];
						mpi_index=mpi_index+1;
					}
					if (i <= IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 26000){
						whole_Qp[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 0){
						whole_Qp[i][j] = mpi_fr[mpi_index];
						mpi_index = mpi_index + 1;
					}
				}
			}
		}
		for (in=1;in<mpi_yarea;in++){
			jn=mpi_jindex;
			rank_num=IMAX*mpi_JMAX[jn];
			if (mpi_porous!=0) rank_recv=(mpi_iindex+in)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex+in)*mpi_zarea+mpi_jindex;
			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					jn1=mpi_j1[mpi_jindex]+j;
					if (jn1<JMAX1||jn1>=JMAX-JMAX2){
						mpi_fs[mpi_index]=whole_Qp[i][j];
						mpi_index=mpi_index+1;
					}
					if (i <= IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 26000){//
						mpi_fs[mpi_index] = whole_Qp[i][j];
						mpi_index = mpi_index + 1;
					}
					if (i >= IMAX - IMAX2 && jn1<JMAX - JMAX2 && jn1 >= JMAX1 && IMAX1 >= 0){
						mpi_fs[mpi_index] = whole_Qp[i][j];
						mpi_index = mpi_index + 1;
					}
				}
			}
			tag=29;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
}

void air::SetQMove(int N)
{
	int i,j,in,jn;
	int temp_I;
	if (N==0){
		for(i=0;i<mpi_IMAX[mpi_iindex];i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if (jn<JMAX1||jn>=JMAX-JMAX2){
					Qv[i][j]=0.0;
					Qw[i][j]=0.0;
					Qp[i][j]=0.0;
				}
				if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){//
					Qv[i][j] = 0.0;
					Qw[i][j] = 0.0;
					Qp[i][j] = 0.0;
				}
				if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
					Qv[i][j] = 0.0;
					Qw[i][j] = 0.0;
					Qp[i][j] = 0.0;
				}
				temp_I=0;
			}
		}
	}else{
		temp_I=move_frame.IMAX-move_frame.lead_DI;
	}
	
	if(temp_I!=0){
		for(i=0;i<mpi_IMAX[mpi_iindex];i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				if (in<temp_I){
					if (jn<JMAX1||jn>=JMAX-JMAX2){
						Qv[i][j]=whole_Qv[in+move_frame.trail_DI-1][j];
						Qw[i][j]=whole_Qw[in+move_frame.trail_DI-1][j];
						Qp[i][j]=whole_Qp[in+move_frame.trail_DI-1][j];
					}
					if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
					if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){//
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
				}else{
					if (jn<JMAX1||jn>=JMAX-JMAX2){
						Qv[i][j]=0.0;
						Qw[i][j]=0.0;
						Qp[i][j]=0.0;
					}
					if (in >= IMAX - IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 0){
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
					if (in <= IMAX2 && jn<JMAX - JMAX2 && jn >= JMAX1 && IMAX1 >= 26000){// 
						Qv[i][j] = 0.0;
						Qw[i][j] = 0.0;
						Qp[i][j] = 0.0;
					}
				}
			}
		}
	}
}
