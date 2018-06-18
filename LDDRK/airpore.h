//++++++++++++++++++++++++++++++filename: airpore.h +++++++++++++++++++++++++++++++++++++++//

//----------------------------description of variable----------------------------------//
/***************************************************************************************/ 
/* the number of points FFT_N be a power of 2 (2^FFT_m)
/* p_t_filename and p_f_filename is filename of pressure output within time-domain and 
/* frequency domain respectively.
/* AirStruct is struct data containing air structure parameter such as density.
/* it is the same for PoreStruct for porous media
/* pr,pi record real and imaginary pressure at receiver
/*************************************************************************************/ 
#ifndef _AIRPORE_H
#define _AIRPORE_H
#include "mygeneral.h"

class calculation_2D;

class airpore
{
private:
	int time_domain,FFT_m,FFT_N,out_difft;
	int CaseNo;
	bool MovingFrame;
	double gauss_width,frequency;
	scheme_list scheme;
	struct MovingFrame move_frame;
	double velocity_coef;
	int velocity_method;
	int restart,restart_out;
	char restart_infile[20];

	calculation_2D *AirMedia;
	calculation_2D *west_pore,*south_pore,*north_pore,*WestSouth_pore,*WestNorth_pore,*east_pore, *EastSouth_pore;

	boundary air_boundary;
	AirStruct AirPara;
	PoreStruct PorePara,AbsorbPara,SideAbsorbPara;
	Position source,receiver;
	BoundaryFrame AirFrame,NorthFrame,SouthFrame,WestFrame,WestSouthFrame,
		WestNorthFrame,EastFrame,EastSouthFrame,EastNorthFrame;
	DifferenceStep AirStep,WestStep,NorthStep,SouthStep,WestSouthStep,
		WestNorthStep,EastStep, EastSouthStep,EastNorthStep;
	
	// add the speed of moving frame
	double speed_moveframe;

	// judge time boundary at left side
    int time_judge;
	// Mpi physics domain decomposition 
	int mpi_iindex,mpi_jindex,mpi_yarea,mpi_zarea,mpi_rank,mpi_size;
	hillpore hill1; // add the porous media for hill;

public:
	airpore();
	airpore(char *input_file);
	~airpore();
	void get_output(int mpi_rank,int mpi_size);
	void save_restartfile(char *restartfile);
	void input_restartfile(int *point_pois,int MF_count);
	void output_coordi();
	void SetInitialCond();
	void SetMovingFrame(int MF_count);
	void Cal_pressure(int time_judge,int temp_time);
	void Cal_velocity();
	void UpdateInitialCond(int MF_count);
	void get_data_contour(int n);
	void mpisend_data_contour();
	void get_FFT_y1(int out_type);
	void get_FFT_y11(int out_type);
	void get_FFT_y12(int out_type);
	void get_FFT_y13(int out_type);
	void get_FFT_y14(int out_type);
	void get_FFT_y15(int out_type);
	void get_FFT_y16(int out_type);
	void get_FFT_y(int out_type);
	void get_FFT_y(int out_type,int No);
	void get_FFT_y();
	void get_FFT_z(void);
	void transfer_output(void);
};

#endif

