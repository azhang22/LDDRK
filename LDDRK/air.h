//++++++++++++++++++++++++++++++filename: air.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/***************************************************************************************/ 
/* Vav,Wav is average speed. Vav_over_y means patial (Vav)/patial (y)
/* Pav is average atmospheric pressure
/* Aav is a reversal of average density
/* adiabatic_coef is adiabatic coefficient
/* Profile_Vav,Profile_Wav point to the average wind speed profile
/**************************************************************************************/ 
#ifndef _AIR_H
#define _AIR_H
#include "mygeneral.h"
#include "calculation_2D.h"

//**************************wind speed profile*****************************//
//z is with respect to local ground level
class SpeedProfile{
private:
	int velocity_method;
	double y0,z0,L,alpha;//m
	double b,c,circulation,vorticity;//b is Mach number for vortex back flow
public:
	SpeedProfile(){};
	SpeedProfile(double velocity_coef,int VM);
	~SpeedProfile();
	double Vav(double x,double y,double z);
	double Wav(double x,double y,double z);
};

class air:public calculation_2D
{
	protected:
		double **Vav_over_y,**Vav_over_z,**Wav_over_y,**Wav_over_z,**Vav,**Wav;
		double **coef1,**coef2;
		double **coef_weighty,**coef_weightz, **gy, **gz, **wy, **wz;
		double *gp1, *gp2, *gq1, *gq2;
		double pow_y1,pow_z1,pow_y2,pow_z2;
		double var_miu1,var_miu2;
		double Pav,adiabatic_coef,coef_tao,coef_eta,coef_y;//sound_speed, Aav;
		int velocity_method;
		SpeedProfile *SP;
		
		double **Qv,**Qw,**Qp;
		double PML_AbsorbZmax,PML_alpha;
		double PML_Height1,PML_z1,PML_Height2,PML_z2;

		// PML boundary prameters for the y direction
		double PML_AbsorbYmax,PML_alphay;
		double PML_Width1,PML_y1,PML_Width2,PML_y2;

		double **whole_Qv,**whole_Qw,**whole_Qp;

		// this is added for porous media hill
		int cr_judge,cr_points;
		double ystar,ystop,zstar,zstop;
		double cr_a,cr_b,cr_c,cr_d,cr_e,cr_f,cr_g,cr_h,cr_i,cr_j,cr_k,cr_l;
		int *Ny, *Nz, *Sy, *Sz, *Ny1, *Nz1, *Sy1, *Sz1, *zone_y, *zone_z;
		// for the hill
		int istar,istop,jstar,jstop;
		double *cr_y,*cr_z;
		int **cr_uc;
		double speed_air1,Aav1,coef_tao1,coef_eta1;
		
		double crpo_eff_density,crpo_porosity,crpo_resistivity;
		double crpo_Kp,crpo_Beta,crpo_Gama;
	public:
		air(){};
		air(scheme_list,DifferenceStep,char *coordi,MovingFrame MF,const double,
			double velocity_coef,int velocity_method1,AirStruct AirPara,double AbsorbZmax,
			PoreStruct PorePara,hillpore hill1,int mpi_rank1,int mpi_size1,int mpi_yarea1,
			int mpi_zarea1,int mpi_porous1);
		~air();
		virtual void cal_fvw(double **temp_v,double **temp_w,
			double **temp_p);
		virtual void cal_fp(double **temp_v,double **temp_w,
			double **temp_p);
		virtual void cal_delvw();
		virtual void shift_grunwald();
		virtual void fractional_CD();
		double Cal_g(int p, double beta);
		double cal_weight(double miu,int n);
		double cal_weight(double y, int l, int s);
		void cal_weight(double y, int K, double *g);
		virtual void save_restart_air(char *restartfile);
		virtual void input_restart_air(char *restartfile,int *point_pois);
		virtual void SetWindProfile(int N);
		virtual void UpdateBC_pressure(boundary_location BC,int time_judge,int time_current);
		virtual void UpdateBC_velocity(boundary_location BC);
		
		virtual void Update_PML_Qvw();
		virtual void Update_PML_Qp();
		virtual void SetQMove(int N);
		virtual void Get_geo();
		virtual void Get_g(int Rank);
};

#endif
