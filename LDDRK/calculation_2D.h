//++++++++++++++++++++++++++++++filename: calculation_2D.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/***************************************************************************************/ 
/* v_n,w_n,p_n will store velocity of y direction, velocity of z direction and pressure 
/* at old time v_nn,w_nn,p_nn will store velocity of y direction, velocity of z direction 
/* and pressure at new time diff_t_air,diff_y_air,diff_z_air is time step, space step 
/* of y direction, space step of z direction Vav_over_y,Vav_over_z,Wav_over_y,Wav_over_z 
/* is differential of average wind speed to axial y and z
/* IMAX,JMAX the max number of points for y and z direction within moving frame
/* IMAX=frame_JMAX reqiure the number is even number so as to be divided as two part when computing
/* whole_IMAX and whole_JMAX are max number of points for whole domain
/* Y_over_y is partial differential of Y to y, where Y is coordinate within physical plane
/* y is coordinate within computational plane.
/* Y_over_z, Z_over_y and Z_over_z have the same meaning as Y_over_y
/*************************************************************************************/ 
#ifndef _CALCULATION_2D_H
#define _CALCULATION_2D_H
#include "mygeneral.h"
#include "string.h"
#include <iostream>
#include <fstream>
using namespace std;
#include "mpi.h"

class calculation_2D 
{
	protected:
		int IMAX,JMAX,JMAX1,JMAX2,IMAX1,IMAX2;
		double diff_t,diff_y,diff_z;//s;m
		scheme_list scheme;
		//enum Col_Method {FB1,FB2,LeapTrap,ABM,RK4,RK2} MyMethod;
		double **v_nn,**w_nn,**p_nn;
		double **v_n,**w_n,**p_n;
//		bool ini_run,alter_run;
		double **Y_over_y,**Y_over_z,**Z_over_y,**Z_over_z,**Jacobian;
		double **Y,**Z;
		double **delvw_nn;
		double gaussian_coef;
		//this is added for moving frame
		double **whole_Y,**whole_Z;
		double **whole_v,**whole_w,**whole_p;
		int whole_IMAX,whole_JMAX,whole_IJKMAX,whole_JMAX1,whole_JMAX2;
		struct MovingFrame move_frame;

		// Different temperature for different speed.
		double *Aav,*speed_air,*temperature_air;
		
		// this is added for reading timeseries at left boundary
		double *time_step,*init_pressure;
		int judge_left,num_data,jjmax1,jjmax2;
		
		//mpi block index: mpi_iindex,mpi_jindex; start and end grids:mpi_i1,mpi_i2,mpi_j1,mpi_j2;
		//Max number of grids:mpi_IMAX, mpi_JMAX; The receiver and send matrix for p and v: pss, 
		//send the pressure data from south boundary; psr receive the pressure data to the sourth boundary. 
		// Other variables are defined defined the same as pss and psr
		int mpi_iindex,mpi_jindex;
		int *mpi_i1,*mpi_i2,*mpi_j1,*mpi_j2,*mpi_IMAX,*mpi_JMAX;
		double *pss,*psr,*pns,*pnr,*pes,*per,*pws,*pwr;
		double *vss,*vsr,*vns,*vnr,*ves,*ver,*vws,*vwr;
		double *wss,*wsr,*wns,*wnr,*wes,*wer,*wws,*wwr;
		double *delvwss,*delvwsr,*delvwns,*delvwnr,*delvwes,*delvwer,*delvwws,*delvwwr;
		double **swp_nn, **swvw_nn;
		// whole domain transfer function for two dimentional variables v_nn,w_nn,p_nn or others;
		double *mpi_fs,*mpi_fr,*mpi_fr1,*mpi_fr2;
		 MPI_Status status;	
		int mpi_rank,mpi_size,mpi_yarea,mpi_zarea,mpi_porous;
//      Left boundary location Y,Z
		double *leftbound_Y,*leftbound_Z,*bottombound_Y;
		double Y_start;

	public:
		calculation_2D();
		calculation_2D(scheme_list,DifferenceStep,char *coordi,MovingFrame move_frame,
			const double,AirStruct,int mpi_rank1,int mpi_size1,int mpi_yarea1,int mpi_zarea1,int mpi_porous1);
		~calculation_2D();
		double cal_InitialPressure(double y,double z);
		void Set_InitialCond(Position source);
		double get_pressure(int ii,int jj,int N);
		int get_whole_IMAX(){return whole_IJKMAX;}
		int get_whole_JMAX(){return whole_JMAX;}

		void get_position(int &ii,int &jj,Position receiver);
		virtual void UpdateBC_pressure(boundary_location BC,int time_judge,int time_current);
		virtual void UpdateBC_velocity(boundary_location BC);
		void UpdateBC_delvw(boundary_location BC);
		virtual void UpdateBC_velocity(boundary_location BC,calculation_2D& porous);
		virtual void UpdateBC_pressure(boundary_location BC,calculation_2D& porous);
		virtual void Update_PML_Qvw(){}
		virtual void Update_PML_Qp(){}
		void UpdateInitialCond(int N);
		void save_restart_cal(char *restartfile);
		void input_restart_cal(char *restart_file,int*point_pois,int N);
		void SetMovingFrame(int N);
		void SetLeftRightMoving(int N);
		void mpi_send_data(void);
	
		virtual void SetWindProfile(int N){}
		virtual void SetQMove(int N){}		
		virtual void save_restart_air(char *restartfile){};
		virtual void save_restart_pml(char *restartfile){};
		virtual void input_restart_air(char *restart_infile,int *point_pois){};
		virtual void input_restart_pml(char *restart_infile,int *point_pois){};
		
		void Cal_velocity(void);
		void Cal_pressure(void);
		void Update_FCD(void);
		void Chg_2D_to_1D(int j1, int j2, double *send_buf, double **u);
		void Chg_1D_to_2D(int j1, int j2, double *rec_buf, double **u);
		virtual void cal_delvw(){};
		virtual void shift_grunwald(){};
		virtual void fractional_CD() {};
		virtual void cal_fvw(double **temp_v,double **temp_w,
			double **temp_p){}
		virtual void cal_fp(double **temp_v,double **temp_w,
			double **temp_p){}
		friend double df_over_dx(double f(double,double,double),double x,
			double y,double z,int Position);
		friend ostream& operator<<(ostream& output_stream, calculation_2D& my_object);
		friend void getOstream(ostream& output_stream, calculation_2D& my_object);
		friend ostream& coordi(ostream& output_stream,calculation_2D& my_object);
};

#endif
