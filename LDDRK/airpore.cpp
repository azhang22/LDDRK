//++++++++++++++++++++++++++++filename: airpore.cpp ++++++++++++++++++++++++++//

//-------------------------------airpore class---------------------------//
/*******************************************************************************/ 
/* the class construct computational environment to calculate
/* the coupling of air and porous 
/* Note: AirStep.diff_y and AirStep.diff_z which is read through input.txt is the space step
/* of the domain above transient domain. Space step of whole air domain is calculated again
/* inside the fuction void airpore::get_coordi().
/* move_frame.IMAX and move_frame.diff_I is the grid points of moving frame and moving 
/* grid points each time
/*******************************************************************************/ 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "airpore.h"
#include "air.h"
#include "porous.h"
#include "time.h"
#include <direct.h>  

//****************************airpore code*******************************//
airpore::airpore(char *input_file)
{
	ifstream infile(input_file);
	infile>>gauss_width;infile.ignore(100,'\n');
	infile>>move_frame.IMAX>>move_frame.lead_DI>>move_frame.trail_DI>>move_frame.Judge;infile.ignore(100,'\n');
	int scheme_type;
	infile>>scheme_type;infile.ignore(100,'\n');
	switch(scheme_type)
	{
	case 1:
		scheme=FB_v_p;break;
	case 2:
		scheme=FB_p_v;break;
	case 3:
		scheme=FB_vp;break;
	case 4:
		scheme=LeapTrap;break;
	case 5:
		scheme=RK4;break;
	case 6:
		scheme=RK2;break;
	case 7:
		scheme=ABM;break;
	default:
		scheme=FB_p_v;
	}

	int boundary_choice[4];
	int i;
	for(i=0;i<4;i++) infile>>boundary_choice[i];
	infile.ignore(100,'\n');
	infile.ignore(100,'\n');
	for(i=0;i<4;i++){
		switch(boundary_choice[i]){
		case 1://rigid
			{
				switch(i){
		case 0:
			air_boundary.air_north=rigid;break;
		case 1:
			air_boundary.air_south=rigid;break;
		case 2:
			air_boundary.air_west=rigid;break;
		case 3:
			air_boundary.air_east=rigid;break;
		default:
			cout<<"please entrance your choice again!";
				}
			}
			break;
		case 2://absorbing layer
			{
				switch(i){
		case 0:
			air_boundary.air_north=absorbing;break;
		case 1:
			air_boundary.air_south=absorbing;break;
		case 2:
			air_boundary.air_west=absorbing;break;
		case 3:
			air_boundary.air_east=absorbing;break;
		default:
			cout<<"please entrance your choice again!";
				}
			}
			break;
		case 3://porous media
			{
				switch(i){
		case 0:
			air_boundary.air_north=porous_media;break;
		case 1:
			air_boundary.air_south=porous_media;break;
		case 2:
			air_boundary.air_west=porous_media;break;
		case 3:
			air_boundary.air_east=porous_media;break;
		default:
			cout<<"please entrance your choice again!";
				}
			}
			break;
		default:
			cout<<"please entrance your choice again!";
		}
	}
	if(air_boundary.air_west==rigid)
	{
		if(air_boundary.air_east==rigid)
		{
			if(air_boundary.air_north==rigid)
			{
				CaseNo=1;
			}
			else
			{//absorbing
				CaseNo=2;
			}
		}
		else
		{
			if(air_boundary.air_south==absorbing)
			{
				CaseNo=6;
			}
			else
			{//absorbing
				CaseNo=5;
			}
		}
	}
	else
	{//absorbing
		if(air_boundary.air_south==absorbing)
		{
			CaseNo=4;
		}
		else
		{
			CaseNo=3;
		}
	}
	infile>>time_domain>>out_difft>>FFT_m>>frequency;infile.ignore(100,'\n');
	FFT_N=1;
	for (i=0;i<FFT_m;i++) FFT_N *= 2;
	//air properties
	infile>>velocity_coef;infile>>velocity_method;infile.ignore(100,'\n');
	//infile>>AirPara.adiabatic_coef>>AirPara.Pav>>AirPara.sound_speed;infile.ignore(100,'\n'); 
	infile>>AirPara.adiabatic_coef>>AirPara.Pav>>AirPara.temperature1>>AirPara.temperature2
												>>AirPara.coef_y>>AirPara.coef_tao>>AirPara.coef_eta;
	infile.ignore(100,'\n');infile.ignore(100,'\n');

	speed_moveframe=331.3*sqrt(1+AirPara.temperature1/273.15);	
	//AirPara.Aav=pow(AirPara.sound_speed,2)/AirPara.adiabatic_coef/AirPara.Pav;
	//absorbing layer structure
	infile>>AbsorbPara.Cs>>AbsorbPara.porosity>>AbsorbPara.resistivity;infile.ignore(100,'\n');
//	AbsorbPara.eff_density=AbsorbPara.Cs/AbsorbPara.porosity;  ///AirPara.Aav;
	AbsorbPara.Kp=1.0/AirPara.adiabatic_coef/AirPara.Pav;
	SideAbsorbPara=AbsorbPara;SideAbsorbPara.resistivity=AbsorbPara.resistivity/4.0;
	//porous media structure
	infile>>PorePara.Cs>>PorePara.porosity>>PorePara.resistivity;infile.ignore(100,'\n');
	//PorePara.eff_density=PorePara.Cs/PorePara.porosity/AirPara.Aav;
	PorePara.eff_density=PorePara.Cs/PorePara.porosity;
	PorePara.Kp=1.0/AirPara.adiabatic_coef/AirPara.Pav;

	//source and receiver 
	infile>>source.y>>source.z>>receiver.y>>receiver.z;infile.ignore(100,'\n');

	//air space dimension
	infile>>AirStep.diff_t;infile.ignore(100,'\n');
	infile>>AirStep.diff_y>>AirStep.diff_z;infile.ignore(100,'\n');
	infile>>AirFrame.left>>AirFrame.right>>AirFrame.upper>>AirFrame.lower;infile.ignore(100,'\n');

	//absorbing boundary setting
	infile>>NorthStep.diff_y>>NorthStep.diff_z;infile.ignore(100,'\n');
	infile>>NorthFrame.left>>NorthFrame.right>>NorthFrame.upper>>NorthFrame.lower;
	infile.ignore(100,'\n');
	infile>>SouthFrame.left>>SouthFrame.right>>SouthFrame.upper>>SouthFrame.lower;
	infile.ignore(100,'\n');
	
	//west frame absorbing boundary
	infile>>WestStep.diff_y>>WestStep.diff_z;infile.ignore(100,'\n');
	infile>>WestFrame.left>>WestFrame.right>>WestFrame.upper>>WestFrame.lower;
	infile.ignore(100,'\n');
	infile>>WestSouthFrame.left>>WestSouthFrame.right>>WestSouthFrame.upper>>WestSouthFrame.lower;
	infile.ignore(100,'\n');
	infile>>WestNorthFrame.left>>WestNorthFrame.right>>WestNorthFrame.upper>>WestNorthFrame.lower;
	infile.ignore(100,'\n');
	
	//east frame absorbing boundary
	infile>>EastStep.diff_y>>EastStep.diff_z;infile.ignore(100,'\n');
	infile>>EastFrame.left>>EastFrame.right>>EastFrame.upper>>EastFrame.lower;
	infile.ignore(100,'\n');
	infile>>EastSouthFrame.left>>EastSouthFrame.right>>EastSouthFrame.upper>>EastSouthFrame.lower;
	infile.ignore(100,'\n');
	infile>>EastNorthFrame.left>>EastNorthFrame.right>>EastNorthFrame.upper>>EastNorthFrame.lower;
	infile.ignore(100,'\n');

	SouthStep.diff_t=WestStep.diff_t=NorthStep.diff_t=
	WestSouthStep.diff_t=WestNorthStep.diff_t=EastStep.diff_t=
	EastSouthStep.diff_t=EastNorthStep.diff_t=AirStep.diff_t;
	SouthStep=NorthStep;
	WestNorthStep.diff_y=WestStep.diff_y;WestNorthStep.diff_z=NorthStep.diff_y;
	EastNorthStep.diff_y=EastStep.diff_y;EastNorthStep.diff_z=NorthStep.diff_y;
	EastSouthStep.diff_y=EastStep.diff_y;EastSouthStep.diff_z=SouthStep.diff_z;
	
	// MPI: Domain blocks yarea*zarea
	infile>>mpi_yarea>>mpi_zarea;infile.ignore(100,'\n');
	
	// output.data
    infile.ignore(100,'\n');infile.ignore(100,'\n');
    infile>>restart>>restart_out>>restart_infile;
	
	// curve setting,adds the porous mediea hill, 
	/* curve_judge: judge the curve or not; curve_coefa,curve_coefb, curve_coefc
	and curve_coefd: The coefficient of equations:y=a*x^3+b*x^2+c*x+d;
	*/
	infile.ignore(100,'\n');infile.ignore(100,'\n');
	infile>>hill1.curve_judge>>hill1.curve_coefa>>hill1.curve_coefb>>
			  hill1.curve_coefc>>hill1.curve_coefd>>hill1.curve_coefe>>
			  hill1.curve_coeff>>hill1.curve_coefg>>hill1.curve_coefh>>
			  hill1.curve_coefi>>hill1.curve_coefj>>hill1.curve_coefk>>
			  hill1.curve_coefl>>hill1.curve_points;
	infile.ignore(100,'\n');
	// The domain we need to calculate for immersed boundary
	infile>>hill1.curve_ystar>>hill1.curve_ystop>>hill1.curve_zstar>>hill1.curve_zstop;

	//porous media setting
	infile.ignore(100,'\n');infile.ignore(100,'\n');
	if(air_boundary.air_south==porous_media)
	{
		infile>>SouthFrame.left>>SouthFrame.right>>SouthFrame.upper>>SouthFrame.lower;
		infile.ignore(100,'\n');
		infile>>SouthStep.diff_y>>SouthStep.diff_z;
		infile.ignore(100,'\n');
		EastSouthStep.diff_z=SouthStep.diff_z;
	}
	else
	{
		infile.ignore(100,'\n');infile.ignore(100,'\n');
		infile.ignore(100,'\n');infile.ignore(100,'\n');
	}
	infile.close();
	//test time step for moving frame according to cp*dx/dt=m
	//AirStep.diff_t=dx*m in input.txt
	/*
	if((AirStep.diff_t/AirStep.diff_z)>0.05){
		AirStep.diff_t=AirStep.diff_t/AirPara.sound_speed;
		SouthStep.diff_t=WestStep.diff_t=NorthStep.diff_t=
			WestSouthStep.diff_t=WestNorthStep.diff_t=EastStep.diff_t=AirStep.diff_t;
	}
	*/
}


airpore::~airpore()
{
}


void airpore::Cal_velocity()
{
	AirMedia->Cal_velocity();	
	switch(CaseNo){
	case 1:
		{	
			AirMedia->UpdateBC_velocity(NorthBC);		
			AirMedia->UpdateBC_velocity(EastBC);
			AirMedia->UpdateBC_velocity(WestBC);
			if(air_boundary.air_south==rigid){
				AirMedia->UpdateBC_velocity(SouthBC);
			}else{
				if (mpi_jindex!=0)AirMedia->UpdateBC_velocity(SouthBC);
				else{
					south_pore->Cal_velocity();
					south_pore->UpdateBC_velocity(SouthBC);
					south_pore->UpdateBC_velocity(EastBC);
					south_pore->UpdateBC_velocity(WestBC);
					AirMedia->UpdateBC_velocity(SouthBC,*south_pore);
				}
			}
			AirMedia->Update_PML_Qvw();
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::Cal_pressure(int time_judge,int temp_time)
{
	AirMedia->Cal_pressure();
	switch(CaseNo){
	case 1:
		{
			AirMedia->UpdateBC_pressure(NorthBC,time_judge,temp_time);
			AirMedia->UpdateBC_pressure(EastBC,time_judge,temp_time);
			AirMedia->UpdateBC_pressure(WestBC,time_judge,temp_time);
			if(air_boundary.air_south==rigid){
				AirMedia->UpdateBC_pressure(SouthBC,time_judge,temp_time);
			}else{
				if (mpi_jindex!=0)AirMedia->UpdateBC_pressure(SouthBC,time_judge,temp_time);
				else{
					south_pore->Cal_pressure();
					south_pore->UpdateBC_pressure(SouthBC,time_judge,temp_time);
					south_pore->UpdateBC_pressure(EastBC,time_judge,temp_time);
					south_pore->UpdateBC_pressure(WestBC,time_judge,temp_time);
					AirMedia->UpdateBC_pressure(SouthBC,*south_pore);
				}
			}
			AirMedia -> Update_PML_Qp();

		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::get_output(int mpi_rank1,int mpi_size1)
{
	//initialize all computational domains
	mpi_rank=mpi_rank1;
	mpi_size=mpi_size1;
	mpi_iindex=mpi_rank/mpi_zarea;
	mpi_jindex=mpi_rank%mpi_zarea;
	SetInitialCond();
	//output coordinate (combine all domains into one coordinate system)
	//output_coordi();//for plot3d in Tecplot

	//main calculation
	int *temp_i,*temp_j,IMAX,JMAX,NO_Rec_y;
	int i,j,mm,index_pressure;
	double **pressure_temp;
	mm=7;
	pressure_temp=new double *[mm];
	for(i=0;i<mm;i++) pressure_temp[i]=new double[time_domain];
	for (i=0;i<mm;i++){
		for (j=0;j<time_domain;j++){
			pressure_temp[i][j]=0.0;
		}
	}
	IMAX=AirMedia->get_whole_IMAX();
	JMAX=AirMedia->get_whole_JMAX();
	if(receiver.y==0){
		NO_Rec_y=1;
	}else{
		NO_Rec_y=mm;
	}
	temp_i=new int [NO_Rec_y];temp_j=new int [NO_Rec_y];
	for(i=0;i<NO_Rec_y;i++){
		int ii,jj;
		struct Position rec;
		if (i == 0){ rec.y = 0.001; rec.z = 0.010; }
		if (i == 1){ rec.y = 0.003; rec.z = 0.010; }
		if (i == 2){ rec.y = 0.005; rec.z = 0.010; }
		if (i == 3){ rec.y = 0.010; rec.z = 0.010; }
		if (i == 4){ rec.y = 0.015; rec.z = 0.010; }
		if (i == 5){ rec.y = 0.020; rec.z = 0.010; }
		if (i == 6){ rec.y = 0.025; rec.z = 0.010; }
//		rec.y=receiver.y/NO_Rec_y*(i+1);rec.z=receiver.z/NO_Rec_y*(i+1);
		AirMedia->get_position(ii,jj,rec);
		if(ii==0&&jj==0){ii=1000000;jj=1000000;}
		temp_i[i]=ii;temp_j[i]=jj;
	}
	//output of p_t(y) -- z is fixed z=temp_j	
	//if pty results files exists, remove them
   	if (restart==0 && mpi_rank==0){
		ofstream ofout("pty1.dat",ios::out|ios::binary);
		if (ofout.is_open()) remove("pty1.dat");

		ofstream ofout1("pty2.dat",ios::out|ios::binary);
		if (ofout1.is_open())remove("pty2.dat");

		ofstream ofout2("pty3.dat",ios::out|ios::binary);
		if (ofout2.is_open()) remove("pty3.dat");

		ofstream ofout3("pty4.dat",ios::out|ios::binary);
		if (ofout3.is_open()) remove("pty4.dat");

		ofstream ofout4("pty5.dat",ios::out|ios::binary);
		if (ofout4.is_open()) remove("pty5.dat");

		ofstream ofout5("pty6.dat",ios::out|ios::binary);
		if (ofout5.is_open()) remove("pty6.dat");

		ofstream ofout6("pty7.dat",ios::out|ios::binary);
		if (ofout6.is_open()) remove("pty7.dat");
		
		ofout.close();ofout1.close();ofout2.close();
		ofout3.close();ofout4.close();ofout5.close();
		ofout6.close();
		_mkdir("Ultrasound_IB_Grunwald_recoding");
		//_mkdir("Ultrasound_IB_Grunwald_recoding\\Strip_v002");
	}
	if (restart==0){
		index_pressure=0;	
	}
	
	int MF_count=0,temp_time=0;//MF_count to count how many frames when moving
	int time_series1,time_series2,time_series3;
	int time_series4,time_series5,time_series6;
	time_series1=450000;time_series2=450000;time_series3=450000;
	time_series4=450000;time_series5=450000;time_series6=450000;
	
	int int_pois,*point_pois;
	if(restart==1){
		char restart_temp1[20];
		sprintf_s(restart_temp1,"%d",mpi_rank);	
		ifstream infile;
                strcat_s(restart_infile,restart_temp1);

		strcat_s(restart_infile,".dat");
		infile.open(restart_infile,ios::in|ios::binary);
		infile>>MF_count>>temp_time>>index_pressure;
		
		infile.ignore(100,'\n');
		int_pois=infile.tellg();
		point_pois=&int_pois;
		infile.close();

		ifstream infile1("pty1.dat",ios::in|ios::binary);		
		ifstream infile2("pty2.dat",ios::in|ios::binary);
		ifstream infile3("pty3.dat",ios::in|ios::binary);
		ifstream infile4("pty4.dat",ios::in|ios::binary);
		ifstream infile5("pty5.dat",ios::in|ios::binary);
		ifstream infile6("pty6.dat",ios::in|ios::binary);
		ifstream infile7("pty7.dat",ios::in|ios::binary);
		
		if (mpi_rank==0){
			for(j=0;j<index_pressure;j++){
				if(j<time_series1){
					infile1>>pressure_temp[0][j];
					infile1.ignore(100,'\n');
				}
				if (j<time_series2){
					infile2>>pressure_temp[1][j];
					infile2.ignore(100,'\n');
				}
				if (j<time_series3){
					infile3>>pressure_temp[2][j];
					infile3.ignore(100,'\n');
				}
			
				if(j<time_series4){
					infile4>>pressure_temp[3][j];
					infile4.ignore(100,'\n');
				}
			
				if (j<time_series5){
	           			infile5>>pressure_temp[4][j];
			        	infile5.ignore(100,'\n');
            			}
				if (j<time_series6){
					infile6>>pressure_temp[5][j];
					infile6.ignore(100,'\n');
				}
	        		infile7>>pressure_temp[6][j];
		    		infile7.ignore(100,'\n');
			}
		}
		infile1.close();infile2.close();infile3.close();
                infile4.close();infile5.close();infile6.close();
		infile7.close();
	}

	char pty_infile[20]=" ";
	strcpy_s(pty_infile,"pty1.dat");
	if (mpi_rank!=0)sprintf_s(pty_infile,"%d",mpi_rank);

	ofstream outfile(pty_infile,ios::app | ios::binary);
	outfile.setf(ios::scientific,ios::floatfield);
	outfile.precision(6);
	if (mpi_rank!=0){
		outfile.close();
		remove(pty_infile);
	}
	if (time_series1<=index_pressure && mpi_rank==0) outfile.close();
	
	char pty_infile1[20]=" ";
	strcpy_s(pty_infile1,"pty2.dat");
	if (mpi_rank!=0)sprintf_s(pty_infile1,"%d",mpi_rank);

	ofstream outfile1(pty_infile1,ios::app | ios::binary);
	outfile1.setf(ios::scientific,ios::floatfield);
	outfile1.precision(6);
	if (mpi_rank!=0){
		outfile1.close();
		remove(pty_infile1);
	}
	if (time_series2<=index_pressure && mpi_rank==0) outfile1.close();

	char pty_infile2[20]=" ";
	strcpy_s(pty_infile2,"pty3.dat");
	if (mpi_rank!=0)sprintf_s(pty_infile2,"%d",mpi_rank);

   	ofstream outfile2(pty_infile2,ios::app | ios::binary);
   	outfile2.setf(ios::scientific,ios::floatfield);
	outfile2.precision(6);
	if (mpi_rank!=0){
		outfile2.close();
		remove(pty_infile2);
	}
	if (time_series3<=index_pressure && mpi_rank==0) outfile2.close();
	char pty_infile3[20]=" ";
	strcpy_s(pty_infile3,"pty4.dat");
	if (mpi_rank!=0)sprintf_s(pty_infile3,"%d",mpi_rank);

   	 ofstream outfile3(pty_infile3,ios::app | ios::binary);
   	 outfile3.setf(ios::scientific,ios::floatfield);
   	 outfile3.precision(6);
	if (mpi_rank!=0){
		outfile3.close();
		remove(pty_infile3);
	}
	if (time_series4<=index_pressure && mpi_rank==0) outfile3.close();

	char pty_infile4[20]=" ";
	strcpy_s(pty_infile4,"pty5.dat");
	if (mpi_rank!=0)sprintf_s(pty_infile4,"%d",mpi_rank);

	ofstream outfile4(pty_infile4,ios::app | ios::binary);
    	outfile4.setf(ios::scientific,ios::floatfield);
    	outfile4.precision(6);
	if (mpi_rank!=0){
		outfile4.close();
		remove(pty_infile4);
	}
	if (time_series5<=index_pressure && mpi_rank==0)outfile4.close();
	
	char pty_infile5[20]=" ";
	strcpy_s(pty_infile5,"pty6.dat");
	if (mpi_rank!=0)sprintf_s(pty_infile5,"%d",mpi_rank);

  	ofstream outfile5(pty_infile5,ios::app | ios::binary);
   	outfile5.setf(ios::scientific,ios::floatfield);
    	outfile5.precision(6);
	if (mpi_rank!=0){
		outfile5.close();
		remove(pty_infile5);
	}
	if (time_series6<=index_pressure && mpi_rank==0)outfile5.close();
	
	char pty_infile6[20]=" ";
	strcpy_s(pty_infile6,"pty7.dat");
	if (mpi_rank!=0)sprintf_s(pty_infile6,"%d",mpi_rank);

	ofstream outfile6(pty_infile6,ios::app | ios::binary);
   	outfile6.setf(ios::scientific,ios::floatfield);
    	outfile6.precision(6);
	if (mpi_rank!=0){
		outfile6.close();
		remove(pty_infile6);
	}

	if (temp_time<=int((5e-6+1.0e-20)/AirStep.diff_t)){  //judge to use the time boundary
		time_judge=1;
	}else{
		time_judge=0;
	}
	do {
		//when count=0, computation will cover whole frame; next, computation will cover half of frame
		int MF_limit;
		double MF_lefttime1;
		double MF_lefttime2;  
		if(MF_count==(int)((IMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.e-6)||(IMAX==move_frame.IMAX)){
			MF_limit=time_domain-temp_time;
		}else if(MF_count==0){
			MF_limit=(int)((move_frame.IMAX-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t*7.0/10.0);
		}else{
			MF_limit=(int)((move_frame.trail_DI-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t);
			MF_lefttime1=0.0;
			MF_lefttime2=0.0;
			if (MF_count%10==0){		
				MF_lefttime1=(move_frame.trail_DI-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t-MF_limit;
				MF_lefttime1=MF_lefttime1*10.0;
			}
			if(MF_count%100==0){
				MF_lefttime2=(move_frame.trail_DI-1)*AirStep.diff_y/speed_moveframe/AirStep.diff_t-MF_limit;
				MF_lefttime2=MF_lefttime2*100-int(MF_lefttime2*10)*10;
			}
			MF_limit=MF_limit+int(MF_lefttime1)+int(MF_lefttime2);
			if ((temp_time+MF_limit)>time_domain)MF_limit=time_domain-temp_time;
		}

		if (restart==1){
			input_restartfile(point_pois,MF_count);
			restart=0;
		}
		//update Y,Z and v,w,p for moving frame(MF)
		SetMovingFrame(MF_count);
		int n;
		for(n=0;n<MF_limit;n++){
			//print out data
			double temp;
			index_pressure=index_pressure+1;
			if(n+temp_time+1<=time_series1){
				temp=AirMedia->get_pressure(temp_i[0],temp_j[0],MF_count);
				if (mpi_rank==0){
					if (fabs(temp)<1.0e-100) temp=0.0;
					outfile<<temp<<endl;
					cout << "Time step = " << index_pressure << ", #1 reading is " << temp;
					pressure_temp[0][index_pressure-1]=temp;
					if (n+temp_time+1==time_series1)outfile.close();
				}
			}
			if(n+temp_time+1<=time_series2){
				temp=AirMedia->get_pressure(temp_i[1],temp_j[1],MF_count);
				if (mpi_rank==0){
					if (fabs(temp)<1.0e-100) temp=0.0;
					outfile1<<temp<<endl;
					pressure_temp[1][index_pressure-1]=temp;
					if(n+temp_time+1==time_series2)outfile1.close();
				}
			}
			if(n+temp_time+1<=time_series3){
				temp=AirMedia->get_pressure(temp_i[2],temp_j[2],MF_count);
				if(mpi_rank==0){
					if (fabs(temp)<1.0e-100) temp=0.0;
					outfile2 << temp << endl;
					pressure_temp[2][index_pressure-1]=temp;
					if (n+temp_time+1==time_series3)outfile2.close();
				}
			}
			if(n+temp_time+1<=time_series4){
				temp=AirMedia->get_pressure(temp_i[3],temp_j[3],MF_count);
				if (mpi_rank==0){
					if (fabs(temp)<1.0e-100) temp=0.0;
					outfile3<<temp<<endl;
					pressure_temp[3][index_pressure-1]=temp;
					if (n+temp_time+1==time_series4)outfile3.close();
				}
			}
			if(n+temp_time+1<=time_series5){
				temp=AirMedia->get_pressure(temp_i[4],temp_j[4],MF_count);
				if (mpi_rank==0){
					if (fabs(temp)<1.0e-100) temp=0.0;
					outfile4<<temp<<endl;
					pressure_temp[4][index_pressure-1]=temp;
					if (n+temp_time+1==time_series5)outfile4.close();
				}
			}
			if(n+temp_time+1<=time_series6){
				temp=AirMedia->get_pressure(temp_i[5],temp_j[5],MF_count);
				if (mpi_rank==0){
					if (fabs(temp)<1.0e-100) temp=0.0;
					outfile5<<temp<<endl;
					pressure_temp[5][index_pressure-1]=temp;
					if (n+temp_time+1==time_series6)outfile5.close();
				}
			}
			temp=AirMedia->get_pressure(temp_i[6],temp_j[6],MF_count);
			if (mpi_rank==0){
				if (fabs(temp)<1.0e-100) temp=0.0;
				outfile6<<temp<<endl;
				pressure_temp[6][index_pressure-1]=temp;
			}
			clock_t t1,t2;
			t1 = clock();
			//calculate pressure and velocity to advance 
			if((scheme==FB_p_v)||(scheme==LeapTrap)){
				Cal_pressure(time_judge,temp_time+n);
				Cal_velocity();
			}else if(scheme==FB_v_p){
				Cal_velocity();
				Cal_pressure(time_judge,temp_time+n);
			}else{
				Cal_pressure(time_judge,temp_time+n); 
				Cal_velocity();
				UpdateInitialCond(MF_count);
				Cal_pressure(time_judge,temp_time+n);
				Cal_velocity();
			}
			t2 = clock();
			if(mpi_rank == 0) cout << ", one time step using time: " << (double(t2 - t1)) / CLOCKS_PER_SEC << " s." << endl;
			//update data for next time
			UpdateInitialCond(MF_count);
			if((n+temp_time+1)%out_difft==0){//output of contour
				if (mpi_iindex==0 && mpi_jindex!=0)mpisend_data_contour();
				if (mpi_rank==0) get_data_contour(n+temp_time);
			}

			if (temp_time+n+1>=int((5e-6+1.0e-20)/AirStep.diff_t) && time_judge==1) time_judge=0;
		}
		MF_count=MF_count+1;
		temp_time=MF_limit+temp_time;
		if(MF_count%restart_out==0){
				char restartfile[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/ct",restart_temp[20];
				sprintf_s(restart_temp,"%d",mpi_rank);
				strcat_s(restartfile,restart_temp);
				int restart_nn;
				restart_nn=MF_count/restart_out;
				if(restart_nn<10)strcat_s(restartfile,"00");
				if((restart_nn>=10)&&(restart_nn<100)) strcat_s(restartfile,"0");
				sprintf_s(restart_temp,"%d",restart_nn);
				strcat_s(restartfile,restart_temp);
				strcat_s(restartfile,".dat");
				ofstream outfile11(restartfile,ios::app|ios::binary);
				outfile11<<MF_count<<" "<<temp_time<<" "<<index_pressure<<endl;
				outfile11.close();
				save_restartfile(restartfile);
				if (mpi_rank==0){
					char restart_pty1[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/pty1",restart_pty2[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/pty2";
					char restart_pty3[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/pty3",restart_pty4[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/pty4";
					char restart_pty5[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/pty5",restart_pty6[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/pty6";
					char restart_pty7[200]="\\20171202_2D_Ultrasound_IB_Grunwald_recoding/Strip_v002/pty7";

					sprintf_s(restart_temp,"%d",restart_nn);
					strcat_s(restart_pty1,restart_temp);
					strcat_s(restart_pty1,".dat");
					strcat_s(restart_pty2,restart_temp);
					strcat_s(restart_pty2,".dat");
					strcat_s(restart_pty3,restart_temp);
					strcat_s(restart_pty3,".dat");
					strcat_s(restart_pty4,restart_temp);
					strcat_s(restart_pty4,".dat");
					strcat_s(restart_pty5,restart_temp);
					strcat_s(restart_pty5,".dat");
					strcat_s(restart_pty6,restart_temp);
					strcat_s(restart_pty6,".dat");
					strcat_s(restart_pty7,restart_temp);
					strcat_s(restart_pty7,".dat");

					ofstream outfile12(restart_pty1,ios::out|ios::binary);
					outfile12.setf(ios::scientific,ios::floatfield);
					outfile12.precision(6);
					ofstream outfile13(restart_pty2,ios::out|ios::binary);
					outfile13.setf(ios::scientific,ios::floatfield);
					outfile13.precision(6);
					ofstream outfile14(restart_pty3,ios::out|ios::binary);
					outfile14.setf(ios::scientific,ios::floatfield);
					outfile14.precision(6);
					ofstream outfile15(restart_pty4,ios::out|ios::binary);
					outfile15.setf(ios::scientific,ios::floatfield);
					outfile15.precision(6);
					ofstream outfile16(restart_pty5,ios::out|ios::binary);
					outfile16.setf(ios::scientific,ios::floatfield);
					outfile16.precision(6);
					ofstream outfile17(restart_pty6,ios::out|ios::binary);
					outfile17.setf(ios::scientific,ios::floatfield);
					outfile17.precision(6);
					ofstream outfile18(restart_pty7,ios::out|ios::binary);
					outfile18.setf(ios::scientific,ios::floatfield);
					outfile18.precision(6);
					for(j=0;j<index_pressure;j++){
						if (j<time_series1)outfile12<<pressure_temp[0][j]<<endl;
						if (j<time_series2)outfile13<<pressure_temp[1][j]<<endl;
						if (j<time_series3)outfile14<<pressure_temp[2][j]<<endl;
						if (j<time_series4)outfile15<<pressure_temp[3][j]<<endl;
						if (j<time_series5)outfile16<<pressure_temp[4][j]<<endl;
						if (j<time_series6)outfile17<<pressure_temp[5][j]<<endl;
						outfile18<<pressure_temp[6][j]<<endl;

					}
					outfile12.close();outfile13.close();outfile14.close();
					outfile15.close();outfile16.close();outfile17.close();
					outfile18.close();
				}
		}

	}while((MF_count<=(int)((IMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.e-6))&&
		(temp_time<time_domain));
	outfile6.close();//outfile2.close();
	//delete new variabe to vacate memory
	delete AirMedia;
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if(air_boundary.air_south==porous_media) delete south_pore;
		}
		break;
	default:
		cout<<"it is wrong to delete object";
	}
	delete[] temp_i; delete[] temp_j;
	for (i=0;i<mm;i++){delete[] pressure_temp[i];}
	delete[] pressure_temp;
}


void airpore::UpdateInitialCond(int MF_count)
{
	AirMedia->UpdateInitialCond(MF_count);
	switch(CaseNo){
			case 1:
			if (mpi_jindex==0){
				if(air_boundary.air_south==porous_media) south_pore->UpdateInitialCond(MF_count);
			}
				break;
			default:
				cout<<"it is wrong to entrance boundary";
	}
}
void airpore::save_restartfile(char *restartfile)
{
	AirMedia->save_restart_cal(restartfile);
	AirMedia->save_restart_air(restartfile);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if (air_boundary.air_south==porous_media) south_pore->save_restart_cal(restartfile);
		}
		break;
	case 2:
		{
			north_pore->save_restart_cal(restartfile);
			if(AbsorbPara.Cs<0){
				north_pore->save_restart_air(restartfile);
				north_pore->save_restart_pml(restartfile);
			}
			if(air_boundary.air_south!=rigid) south_pore->save_restart_cal(restartfile);
			if((AbsorbPara.Cs<0)&&(air_boundary.air_south==absorbing)){
				south_pore->save_restart_air(restartfile);
				south_pore->save_restart_pml(restartfile);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}
void airpore::input_restartfile(int *point_pois,int MF_count)
{
	AirMedia->input_restart_cal(restart_infile,point_pois,MF_count);
	AirMedia->input_restart_air(restart_infile,point_pois);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if (air_boundary.air_south==porous_media) south_pore->input_restart_cal(restart_infile,point_pois,MF_count);
		}
		break;
	case 2:
		{
			north_pore->input_restart_cal(restart_infile,point_pois,MF_count);
			if(AbsorbPara.Cs<0){
				north_pore->input_restart_air(restart_infile,point_pois);
				north_pore->input_restart_pml(restart_infile,point_pois);
			}
			if(air_boundary.air_south!=rigid) south_pore->input_restart_cal(restart_infile,point_pois,MF_count);
			if((AbsorbPara.Cs<0)&&(air_boundary.air_south==absorbing)){
				south_pore->input_restart_air(restart_infile,point_pois);
				south_pore->input_restart_pml(restart_infile,point_pois);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}

}


void airpore::SetMovingFrame(int MF_count)
{
	AirMedia->SetMovingFrame(MF_count);
	AirMedia->SetWindProfile(MF_count);
	AirMedia->SetQMove(MF_count);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if(air_boundary.air_south==porous_media) south_pore->SetMovingFrame(MF_count);
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}


void airpore::SetInitialCond()
{
	char filename[100];
	strcpy_s(filename,"coordi_air.dat");
	AirMedia= new air(scheme,AirStep,filename,move_frame,gauss_width,
		velocity_coef,velocity_method,AirPara,AbsorbPara.resistivity,PorePara,hill1,
		mpi_rank,mpi_size,mpi_yarea,mpi_zarea,0);
	AirMedia->Set_InitialCond(source);
	switch(CaseNo){
	case 1:
		if (mpi_jindex==0){
			if(air_boundary.air_south==porous_media){
				strcpy_s(filename,"coordi_pore.dat");
				south_pore=new porous(scheme,SouthStep,filename,move_frame,gauss_width,PorePara,AirPara,
																mpi_rank,mpi_size,mpi_yarea,1,mpi_zarea);			
				south_pore->Set_InitialCond(source);
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void airpore::get_FFT_y1(int out_type)
{ 
	int i;
	double dt;
	dt=AirStep.diff_t;
	ifstream infile("pty7.dat",ios::in|ios::binary);
	double *ppr,*pr,*pi;
   	ppr=new double [time_domain];
	pi=new double [FFT_N];
	pr=new double [FFT_N];
	for (i=0;i<time_domain;i++)ppr[i]=0.0;
	for (i=0;i<FFT_N;i++){
		pr[i]=0.0;
		pi[i]=0.0;
	}
	for (i=0;i<time_domain;i++){
		infile>>ppr[i];
		infile.ignore(100,'\n');
   	 }
    	char temp[50],p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,")7.dat");
	strcat_s(p_f_file,temp);	

	if(receiver.y!=0){
		for(i=0;i<time_domain;i++) {pr[i]=ppr[i];pi[i]=0;}
		for(i=time_domain;i<FFT_N;i++) {pr[i]=0;pi[i]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
	delete[] pi;delete[] pr;delete[] ppr;
}
void airpore::get_FFT_y11(int out_type)
{ 
	int i;
	double dt;
	dt=AirStep.diff_t;
	ifstream infile("pty6.dat",ios::in|ios::binary);
	double *ppr,*pr,*pi;
   	ppr=new double [time_domain];
	pi=new double [FFT_N];
	pr=new double [FFT_N];
	for (i=0;i<time_domain;i++)ppr[i]=0.0;
	for (i=0;i<FFT_N;i++){
		pr[i]=0.0;
		pi[i]=0.0;
	}
	for (i=0;i<time_domain;i++){
		infile>>ppr[i];
		infile.ignore(100,'\n');
   	 }
    	char temp[50],p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,")6.dat");
	strcat_s(p_f_file,temp);	

	if(receiver.y!=0){
		for(i=0;i<time_domain;i++) {pr[i]=ppr[i];pi[i]=0;}
		for(i=time_domain;i<FFT_N;i++) {pr[i]=0;pi[i]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
	delete[] pi;delete[] pr;delete[] ppr;
}

void airpore::get_FFT_y12(int out_type)
{ 
	int i;
	double dt;
	dt=AirStep.diff_t;
	ifstream infile("pty5.dat",ios::in|ios::binary);
	double *ppr,*pr,*pi;
   	ppr=new double [time_domain];
	pi=new double [FFT_N];
	pr=new double [FFT_N];
	for (i=0;i<time_domain;i++)ppr[i]=0.0;
	for (i=0;i<FFT_N;i++){
		pr[i]=0.0;
		pi[i]=0.0;
	}
	for (i=0;i<time_domain;i++){
		infile>>ppr[i];
		infile.ignore(100,'\n');
   	 }
    	char temp[50],p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,")5.dat");
	strcat_s(p_f_file,temp);	

	if(receiver.y!=0){
		for(i=0;i<time_domain;i++) {pr[i]=ppr[i];pi[i]=0;}
		for(i=time_domain;i<FFT_N;i++) {pr[i]=0;pi[i]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
	delete[] pi;delete[] pr;delete[] ppr;
}

void airpore::get_FFT_y13(int out_type)
{ 
	int i;
	double dt;
	dt=AirStep.diff_t;
	ifstream infile("pty4.dat",ios::in|ios::binary);
	double *ppr,*pr,*pi;
   	ppr=new double [time_domain];
	pi=new double [FFT_N];
	pr=new double [FFT_N];
	for (i=0;i<time_domain;i++)ppr[i]=0.0;
	for (i=0;i<FFT_N;i++){
		pr[i]=0.0;
		pi[i]=0.0;
	}
	for (i=0;i<time_domain;i++){
		infile>>ppr[i];
		infile.ignore(100,'\n');
   	 }
    	char temp[50],p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,")4.dat");
	strcat_s(p_f_file,temp);	

	if(receiver.y!=0){
		for(i=0;i<time_domain;i++) {pr[i]=ppr[i];pi[i]=0;}
		for(i=time_domain;i<FFT_N;i++) {pr[i]=0;pi[i]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
	delete[] pi;delete[] pr;delete[] ppr;
}
void airpore::get_FFT_y14(int out_type)
{ 
	int i;
	double dt;
	dt=AirStep.diff_t;
	ifstream infile("pty3.dat",ios::in|ios::binary);
	double *ppr,*pr,*pi;
   	ppr=new double [time_domain];
	pi=new double [FFT_N];
	pr=new double [FFT_N];
	for (i=0;i<time_domain;i++)ppr[i]=0.0;
	for (i=0;i<FFT_N;i++){
		pr[i]=0.0;
		pi[i]=0.0;
	}
	for (i=0;i<time_domain;i++){
		infile>>ppr[i];
		infile.ignore(100,'\n');
   	 }
    	char temp[50],p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,")3.dat");
	strcat_s(p_f_file,temp);	

	if(receiver.y!=0){
		for(i=0;i<time_domain;i++) {pr[i]=ppr[i];pi[i]=0;}
		for(i=time_domain;i<FFT_N;i++) {pr[i]=0;pi[i]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
	delete[] pi;delete[] pr;delete[] ppr;
}

void airpore::get_FFT_y15(int out_type)
{ 
	int i;
	double dt;
	dt=AirStep.diff_t;
	ifstream infile("pty2.dat",ios::in|ios::binary);
	double *ppr,*pr,*pi;
   	ppr=new double [time_domain];
	pi=new double [FFT_N];
	pr=new double [FFT_N];
	for (i=0;i<time_domain;i++)ppr[i]=0.0;
	for (i=0;i<FFT_N;i++){
		pr[i]=0.0;
		pi[i]=0.0;
	}
	for (i=0;i<time_domain;i++){
		infile>>ppr[i];
		infile.ignore(100,'\n');
   	 }
    	char temp[50],p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,")2.dat");
	strcat_s(p_f_file,temp);	

	if(receiver.y!=0){
		for(i=0;i<time_domain;i++) {pr[i]=ppr[i];pi[i]=0;}
		for(i=time_domain;i<FFT_N;i++) {pr[i]=0;pi[i]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
	delete[] pi;delete[] pr;delete[] ppr;
}

void airpore::get_FFT_y16(int out_type)
{ 
	int i;
	double dt;
	dt=AirStep.diff_t;
	ifstream infile("pty1.dat",ios::in|ios::binary);
	double *ppr,*pr,*pi;
   	ppr=new double [time_domain];
	pi=new double [FFT_N];
	pr=new double [FFT_N];
	for (i=0;i<time_domain;i++)ppr[i]=0.0;
	for (i=0;i<FFT_N;i++){
		pr[i]=0.0;
		pi[i]=0.0;
	}
	for (i=0;i<time_domain;i++){
		infile>>ppr[i];
		infile.ignore(100,'\n');
   	 }
    	char temp[50],p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,")1.dat");
	strcat_s(p_f_file,temp);	

	if(receiver.y!=0){
		for(i=0;i<time_domain;i++) {pr[i]=ppr[i];pi[i]=0;}
		for(i=time_domain;i<FFT_N;i++) {pr[i]=0;pi[i]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
	delete[] pi;delete[] pr;delete[] ppr;
}
void airpore::get_FFT_y(int out_type)
//require big memory 
//out_type=1, only output FFT for single point;if not, output FFT for p(r,) and for single point
{
	int i,j;
	//read data from files of coordinate
	ifstream myfile("coordi_air.dat",ios::in|ios::binary);
	if(!myfile){
		cout<<"cannot open file:"<<" for read!!!"<<endl;
		return;
	}
	int IMAX,JMAX,KMAX;
	double **Y,**Z;
	myfile>>IMAX>>JMAX>>KMAX;
	Y=new double* [IMAX];Z=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		Y[i]=new double [JMAX];Z[i]=new double [JMAX];
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Y[i][j];
		}
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Z[i][j];
		}
	}
	myfile.close();
//	receiver.y=receiver.y/2.0;
//	receiver.z=receiver.z/2.0;

	//get temp_i,temp_j
	int temp_i,temp_j;
	for(i=0;i<IMAX;i++){
		for(j=0;j<JMAX;j++){
			if(((receiver.y-Y[i][j])<(AirStep.diff_y/8))&&
				((receiver.z-Z[i][j])<(AirStep.diff_z/8)))
			{
				temp_i=i;temp_j=j;
				i=IMAX;j=JMAX;
			}
		}
	}

	//get p_t file to FFT
	int times,whole_IMAX;
	double dt;
	ifstream infile("pty7.dat",ios::in|ios::binary);
	infile>>times>>whole_IMAX>>dt;infile.ignore(100,'\n');infile.ignore(100,'\n');
	//ppr to record pressure at position (r,temp_j) for all times;ppr-real part,
	//pi-imaginary part,pr-record pressure at (temp_i,temp_j) for all times
    
	double **ppr,*pi,*pr;
	ppr=new double* [time_domain];pi=new double [FFT_N];pr=new double [FFT_N];
	for(i=0;i<time_domain;i++){
		ppr[i]=new double [whole_IMAX];
	}
	for(i=0;i<time_domain;i++){
		for(j=0;j<whole_IMAX;j++){
			infile>>ppr[i][j];
		}
	}
	infile.close();

	//compute FFT
	int n,nn;
	//FFT for different position like p(r,0.5)
	char temp[50],p_d_file[100]="p_y(";
	if(out_type!=1)
	{
		sprintf_s(temp,"%d",(int)frequency);
		strcat_s(temp,"Hz).dat");
		strcat_s(p_d_file,temp);
		ofstream outfile1(p_d_file,ios::out | ios::trunc | ios::binary);
		nn=(int)(frequency*FFT_N*dt);
		for(i=0;i<whole_IMAX;i++){
			for(n=0;n<time_domain;n++){ pr[n]=ppr[n][i];pi[n]=0;}
			for(n=time_domain;n<FFT_N;n++) {pr[n]=0;pi[n]=0;}
			FFT(1,FFT_m,pr,pi);
			outfile1<<Y[i][temp_j]<<" "<<Z[i][temp_j]<<" ";
			Complex comp;
			comp=Complex(pr[nn],pi[nn]);
			outfile1<<abs(comp)<<" "<<arg(comp)<<endl;
		}
		outfile1.close();
	}

	//filename for p_t and p_f
	char p_t_file[100]="p_t(",p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,").dat");
	strcat_s(p_t_file,temp);
	strcat_s(p_f_file,temp);
	//output p_t for some position
	ofstream outfile2(p_t_file,ios::out | ios::trunc | ios::binary);
	if(receiver.y!=0){
		for(n=0;n<time_domain;n++) {
			outfile2<<dt*n<<" "<<ppr[n][temp_i]<<endl;
		}
	}
	outfile2.close();
	//FFT for single point
	if(receiver.y!=0){
		for(n=0;n<time_domain;n++) {pr[n]=ppr[n][temp_i];pi[n]=0;}
		for(n=time_domain;n<FFT_N;n++) {pr[n]=0;pi[n]=0;}
		FFT(1,FFT_m,pr,pi);
		FFT_output(FFT_m,pr,pi,p_f_file,dt);
	}
}


void airpore::get_FFT_y(int out_type,int No)
//require big memory 
//out_type=1, only output FFT for single point;if not, output FFT for p(r,) and for single point
//No is the number of recording data along z. see NO_Rec_y in get_output().
{
	int i,j;
	//read data from files of coordinate
	ifstream myfile("coordi_air.dat",ios::in|ios::binary);
	if(!myfile){
		cout<<"cannot open file:"<<" for read!!!"<<endl;
		return;
	}
	int IMAX,JMAX,KMAX;
	double **Y,**Z;
	myfile>>IMAX>>JMAX>>KMAX;
	Y=new double* [IMAX];Z=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		Y[i]=new double [JMAX];Z[i]=new double [JMAX];
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Y[i][j];
		}
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Z[i][j];
		}
	}
	myfile.close();
	//get temp_i,temp_j
	int temp_i,temp_j;
	for(i=0;i<IMAX;i++){
		for(j=0;j<JMAX;j++){
			if(((receiver.y-Y[i][j])<(AirStep.diff_y/8))&&
				((receiver.z-Z[i][j])<(AirStep.diff_z/8))){
					temp_i=i;temp_j=j;
					i=IMAX;j=JMAX;
				}
		}
	}
	//get p_t file to FFT
	int times,whole_IMAX;
	double dt;
	ifstream infile("pty.dat",ios::in|ios::binary);
	infile>>times>>whole_IMAX>>dt;infile.ignore(100,'\n');infile.ignore(100,'\n');
	//ppr to record pressure at position (r,temp_j) for all times;ppr-real part,
	//pi-imaginary part,pr-record pressure at (temp_i,temp_j) for all times
	double **ppr,*pi,*pr;
	ppr=new double* [time_domain];pi=new double [FFT_N];pr=new double [FFT_N];
	for(i=0;i<time_domain;i++){
		ppr[i]=new double [whole_IMAX/No+1];
	}
	int n,nn;
	char temp[50],p_d_file[100]="p_y";
	sprintf_s(temp,"%d",(int)frequency);
	strcat_s(temp,"Hz.dat");
	strcat_s(p_d_file,temp);
	ofstream outfile1(p_d_file,ios::out | ios::trunc | ios::binary);
	nn=(int)(frequency*FFT_N*dt);
	//FFT for different position like p(r,0.5)
	for(int jjj=0;jjj<No;jjj++){
		if(jjj!=(No-1)){
			IMAX=whole_IMAX/No+1;
		}else{
			IMAX=whole_IMAX-IMAX*(No-1);
		}
		for(i=0;i<time_domain;i++){
			for(j=0;j<IMAX;j++){
				infile>>ppr[i][j];infile.ignore(whole_IMAX,'\n');
			}
		}
		for(i=0;i<IMAX;i++){
			for(n=0;n<time_domain;n++){ pr[n]=ppr[n][i];pi[n]=0;}
			for(n=time_domain;n<FFT_N;n++) {pr[n]=0;pi[n]=0;}
			FFT(1,FFT_m,pr,pi);
			outfile1<<Y[i][temp_j]<<" "<<Z[i][temp_j]<<" ";
			Complex comp;
			comp=Complex(pr[nn],pi[nn]);
			outfile1<<abs(comp)<<" "<<arg(comp)<<endl;
		}
	}
	infile.close();outfile1.close();
}

void airpore::get_FFT_y()//no big memory requirement
{
	int i,j;
	//read data from files of coordinate
	ifstream myfile("coordi_air.dat",ios::in|ios::binary);
	if(!myfile){
		cout<<"cannot open file:"<<" for read!!!"<<endl;
		return;
	}
       
   	int IMAX,JMAX,KMAX;
	double **Y,**Z;
	myfile>>IMAX>>JMAX>>KMAX;
	Y=new double* [IMAX];Z=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		Y[i]=new double [JMAX];Z[i]=new double [JMAX];
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Y[i][j];
		}
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Z[i][j];
		}
	}
	myfile.close();
	//get temp_i,temp_j
	int temp_i,temp_j;

	for(i=0;i<IMAX;i++){
		for(j=0;j<JMAX;j++){
			if(((receiver.y-Y[i][j])<(AirStep.diff_y/8))&&
				((receiver.z-Z[i][j])<(AirStep.diff_z/8))){
					temp_i=i;temp_j=j;
					i=IMAX;j=JMAX;
				}
		}
	}
	//get p_t file to FFT
	int times,whole_IMAX;
	double dt;
	ifstream infile;
	infile.open("pty1.dat",ios::in|ios::binary);
	infile>>times>>whole_IMAX>>dt;infile.ignore(100,'\n');infile.ignore(100,'\n');
	//get data position property
	int infile_pos,inter_pos,temp_pos1,temp_pos2;
	double temp_data;
	infile_pos=infile.tellg();infile>>temp_data;temp_pos1=infile.tellg();
	infile>>temp_data;temp_pos2=infile.tellg();inter_pos=temp_pos2-temp_pos1;
       // printf("%d,%d,%d,%d,%f\n",infile_pos,inter_pos,temp_pos1,temp_pos2,temp_data);//add term
	//ppr to record pressure at position (r,temp_j) for all times;ppr-real part,
	//pi-imaginary part,pr-record pressure at (temp_i,temp_j) for all times
        double *pi,*pr;
	pi=new double [FFT_N];pr=new double [FFT_N];
	//filename for signal p along y-direction 
	char temp[50],p_d_file[100]="p_y";
	sprintf_s(temp,"%d",(int)frequency);
	strcat_s(temp,"Hz.dat");
	strcat_s(p_d_file,temp);
	//filename for single point p_f
	char p_f_file[100]="p_f(";
	sprintf_s(temp,"%d,%d",(int)receiver.y,(int)receiver.z);
	strcat_s(temp,").dat");
	strcat_s(p_f_file,temp);
	int n,nn;
	for(n=time_domain;n<FFT_N;n++) {pr[n]=0;pi[n]=0;}
	//FFT for different position like p(r,0.5)
	ofstream outfile(p_d_file,ios::out | ios::trunc | ios::binary);
	nn=(int)(frequency*FFT_N*dt);

	for(i=0;i<whole_IMAX;i++){
		int ini_pos;
		ini_pos=infile_pos;
		if(i!=0){
			infile.close();
			infile.open("p_t(y).dat",ios::in|ios::binary);
		}
		printf("i=%d/n",i);// add term
		for(n=0;n<time_domain;n++){
			infile.seekg(i*inter_pos+ini_pos);
		//------add term
			if (n==0||n==1||n==2){
                           printf("i=%d\n",i);
			   printf("n=%d,%d,%d,%d\n",n,inter_pos,ini_pos,i*inter_pos+ini_pos);
                          };			
                       
		//------add term
         		infile>>pr[n];pi[n]=0;
		//------add term
                 	if (n==0||n==1||n==2||n==3||n==4||n==5||n==6){
			   printf("n=%d\n", n);
                        }
		//------add term	
           
			infile.ignore(whole_IMAX*inter_pos,'\n');
			ini_pos=infile.tellg();
		}
		FFT(1,FFT_m,pr,pi);
		if((receiver.y!=0)&&(i==temp_i)){//FFT for single point
			FFT_output(FFT_m,pr,pi,p_f_file,dt);
		}
		outfile<<Y[i][temp_j]<<" "<<Z[i][temp_j]<<" ";
		Complex comp;
		comp=Complex(pr[nn],pi[nn]);
		outfile<<abs(comp)<<" "<<arg(comp)<<endl;
	}
	infile.close();outfile.close();
}

void airpore::transfer_output(void)
{
	//get p_t file to FFT
	int times,whole_IMAX;
	double dt;
	ifstream infile;
	infile.open("pty1.dat",ios::in|ios::binary);
	infile>>times>>whole_IMAX>>dt;infile.ignore(100,'\n');infile.ignore(100,'\n');

	ofstream outfile("pty_modi.dat",ios::out | ios::trunc | ios::binary);
	outfile<<times<<" "<<whole_IMAX<<" "<<dt<<endl;
	outfile<<100<<" "<<100<<endl;
	outfile.setf(ios::scientific,ios::floatfield);
	outfile.setf(ios::showpos);
	outfile.precision(6);
	outfile.width(14);
	for(int n=0;n<time_domain;n++){
		for(int i=0;i<whole_IMAX;i++){
			double a;
			infile>>a;
			outfile<<a<<" ";
		}
		outfile<<endl;
	}
	infile.close();outfile.close();
}

void airpore::get_FFT_z(void)
//require big memory 
//output FFT for p(,r) and for single point
{
	int i,j,k;
	//read data from files of coordinate
	ifstream myfile("coordi_air.dat",ios::in|ios::binary);
	if(!myfile){
		cout<<"cannot open file:"<<" for read!!!"<<endl;
		return;
	}
	int IMAX,JMAX,KMAX;
	double **Y,**Z;
	myfile>>IMAX>>JMAX>>KMAX;
	Y=new double* [IMAX];Z=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		Y[i]=new double [JMAX];Z[i]=new double [JMAX];
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Y[i][j];
		}
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<IMAX;i++){
			myfile>>Z[i][j];
		}
	}
	myfile.close();

	//get temp_i,temp_j
	int *temp_i,temp_j,NO_Rec_y;
	if(receiver.y==0){
		NO_Rec_y=1;
	}else{
		NO_Rec_y=1;
	}
	temp_i=new int [NO_Rec_y];
	for(k=0;k<NO_Rec_y;k++){
		struct Position rec;
		rec.y=receiver.y/NO_Rec_y*(k+1);rec.z=receiver.z;
		for(i=0;i<IMAX;i++){
			for(j=0;j<JMAX;j++){
				if(((receiver.y-Y[i][j])<(AirStep.diff_y/8))&&
					((receiver.z-Z[i][j])<(AirStep.diff_z/8))){
						temp_i[k]=i;temp_j=j;
						i=IMAX;j=JMAX;
					}
			}
		}
	}

	//get p_t file to FFT
	int times,whole_JMAX;
	double dt;
	ifstream infile("ptz.dat",ios::in|ios::binary);
	infile>>times>>whole_JMAX>>dt;infile.ignore(100,'\n');infile.ignore(100,'\n');
	//ppr to record pressure at position (r,temp_j) for all times;ppr-real part,
	//pi-imaginary part,pr-record pressure at (temp_i,temp_j) for all times
	double ***ppr,*pi,*pr;
	pi=new double [FFT_N];pr=new double [FFT_N];
	ppr=new double** [NO_Rec_y];
	for(k=0;k<NO_Rec_y;k++){
		ppr[k]=new double* [time_domain];
		for(i=0;i<time_domain;i++){
			ppr[k][i]=new double [whole_JMAX];
		}
	}
	for(i=0;i<time_domain;i++){
		for(k=0;k<NO_Rec_y;k++){
			for(j=0;j<whole_JMAX;j++){
				infile>>ppr[k][i][j];
			}
		}
	}
	infile.close();

	//compute FFT
	int n,nn;
	//FFT for different position like p(0.5,z)
	char temp[50],p_d_file[100]="p_z";
	sprintf_s(temp,"%d",(int)frequency);
	strcat_s(temp,"Hz.dat");
	strcat_s(p_d_file,temp);
	ofstream outfile1(p_d_file,ios::out | ios::trunc | ios::binary);
	nn=(int)(frequency*FFT_N*dt);

	for(k=0;k<NO_Rec_y;k++){
		outfile1<<"Y Z receiver("<<receiver.y/NO_Rec_y*(k+1)<<","<<receiver.z<<")-amplitude -phase"<<" ";
	}
	outfile1<<endl;
	for(i=0;i<whole_JMAX;i++){
		for(k=0;k<NO_Rec_y;k++){
			for(n=0;n<time_domain;n++){ pr[n]=ppr[k][n][i];pi[n]=0;}
			for(n=time_domain;n<FFT_N;n++) {pr[n]=0;pi[n]=0;}
			FFT(1,FFT_m,pr,pi);
			outfile1<<Y[temp_i[k]][i]<<" "<<Z[temp_i[k]][i]<<" ";
			Complex comp;
			comp=Complex(pr[nn],pi[nn]);
			outfile1<<abs(comp)<<" "<<arg(comp)<<endl;
		}
		outfile1<<endl;
	}
	outfile1.close();
}

void airpore::mpisend_data_contour(){		
	AirMedia->mpi_send_data();
}

//the following is the output format for the command "load data files" in Tecplot
void airpore::get_data_contour(int n)
{
	char p_contour[200]="Ultrasound_IB_Grunwald_recoding\\pt",temp[20];
	int nn;
	nn=(n+1)/out_difft;
	if(nn<10) strcat_s(p_contour,"00");
	if((nn>=10)&&(nn<100)) strcat_s(p_contour,"0");
	int ii;
	ii=sprintf_s(temp,"%d",nn);
	strcat_s(p_contour,temp);strcat_s(p_contour,".dat");

	//output
	ofstream fp_p_t(p_contour,ios::out | ios::trunc | ios::binary);
	fp_p_t.setf(ios::scientific,ios::floatfield);
	fp_p_t.precision(6);

	switch(CaseNo){
	case 1:
		{
			if(air_boundary.air_south==rigid){
				fp_p_t<<"ZONE T = \"Air Zone\",";
				fp_p_t<<*AirMedia;
			}else{

				fp_p_t<<"ZONE T = \"Air Zone\",";
				fp_p_t<<*AirMedia;
				fp_p_t<<"ZONE T = \"south pore Zone\",";
				fp_p_t<<*south_pore;
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
	fp_p_t.close();
}
