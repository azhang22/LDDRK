//+++++++++++++++++++++++filename: calculation_2D.cpp +++++++++++++++++++++++++++++//

//-----------------------fathe class for air and porous class------------------//
/*******************************************************************************/ 
/* air and porous class will inherit from the class calculation_2D. 
/*******************************************************************************/ 
#include <stdio.h>
#include <math.h>
#include "calculation_2D.h"

calculation_2D::calculation_2D(scheme_list scheme1,DifferenceStep Step,char *coordi,
					MovingFrame MF,const double gauss_width,AirStruct AirPara,int mpi_rank1,
					int mpi_size1,int mpi_yarea1,int mpi_zarea1,int mpi_porous1)
{
	int mpi_ny,mpi_nz;
	mpi_rank=mpi_rank1;
	mpi_size=mpi_size1;
	mpi_yarea=mpi_yarea1;
	mpi_zarea=mpi_zarea1;
	mpi_porous=mpi_porous1;
	scheme=scheme1;
	move_frame=MF;
	diff_t=Step.diff_t;
	diff_y=Step.diff_y;
	diff_z=Step.diff_z;
	gaussian_coef=4.0*log(2.0)/pow(gauss_width,2);
	//gaussian_coef=4.0*log(2.0)/pow(7.5e-3,2);
	//gaussian_coef=4.0*log(2.0)/pow(0.6,2);
	//gaussian_coef=4.0*log(2.0)/pow(5.0e-3,2);
	//gaussian_coef=40*7*7/pow(550*diff_y,2);
	//gaussian_coef=40*7*7/pow(400*diff_y,2);

	int i,j;//,k
	//read data from files
	ifstream myfile(coordi,ios::in|ios::binary);
	if(!myfile){
		cout<<"cannot open file:"<<coordi<<" for read!!!"<<endl;
		return;
	}
	int whole_KMAX;
	myfile>>whole_IMAX>>whole_JMAX>>whole_JMAX1>>whole_JMAX2>>whole_KMAX>>whole_IJKMAX;
	//the following setting is for one moving frame
	int N_limit;
	if(MF.Judge==1){
		if ((whole_IJKMAX-move_frame.IMAX)%(move_frame.lead_DI-1)!=0){
			N_limit=int((whole_IJKMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.0e-10);
			whole_IJKMAX=move_frame.IMAX+(N_limit+1)*(move_frame.lead_DI-1);
		}
		IMAX=move_frame.IMAX;
		JMAX=whole_JMAX;
		JMAX1=whole_JMAX1;
		JMAX2=whole_JMAX2;
		IMAX2=whole_KMAX;
		IMAX1=0;
	}else{
		IMAX=whole_IMAX;JMAX=whole_JMAX;
		JMAX1=whole_JMAX1;JMAX2=whole_JMAX2;
		IMAX2=whole_KMAX;
		IMAX1=0;
		whole_IJKMAX=IMAX;
	}
	// mpi:specify the domain decomposition
	mpi_ny=int(IMAX/mpi_yarea);
	mpi_nz=int(JMAX/mpi_zarea);
	mpi_iindex=mpi_rank/mpi_zarea;
	mpi_jindex=mpi_rank%mpi_zarea;
	if (mpi_porous!=0){
		mpi_iindex=mpi_rank/mpi_porous;
		mpi_jindex=mpi_rank%mpi_porous;
	}
	mpi_i1=new int[mpi_yarea];
	mpi_i2=new int[mpi_yarea];
	mpi_j1=new int[mpi_zarea];
	mpi_j2=new int[mpi_zarea];
	mpi_IMAX=new int[mpi_yarea];
	mpi_JMAX=new int[mpi_zarea];
	for (i=0;i<mpi_yarea;i++){
		mpi_i1[i]=i*mpi_ny;
		if (i!=mpi_yarea-1) mpi_i2[i]=(i+1)*mpi_ny+1;
		else mpi_i2[i]=IMAX-1;
		mpi_IMAX[i]=mpi_i2[i]-mpi_i1[i]+1; 
	}
	for (j=0;j<mpi_zarea;j++){
		mpi_j1[j]=j*mpi_nz;
		if (j!=mpi_zarea-1) mpi_j2[j]=(j+1)*mpi_nz+1;
		else mpi_j2[j]=JMAX-1;
		mpi_JMAX[j]=mpi_j2[j]-mpi_j1[j]+1;
	}
	whole_Y=new double* [IMAX];whole_Z=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		whole_Y[i]=new double [mpi_JMAX[mpi_jindex]];
		whole_Z[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	leftbound_Y=new double[JMAX];leftbound_Z=new double[JMAX];
	bottombound_Y = new double[whole_IMAX];
	int jn;
	double temp_Y,temp_Z;
	for(j=0;j<JMAX;j++){
		for(i=0;i<whole_IMAX;i++){
			myfile>>temp_Y;
			if (j>=mpi_j1[mpi_jindex]&&j<=mpi_j2[mpi_jindex]){
				jn=j-mpi_j1[mpi_jindex];
				whole_Y[i][jn]=temp_Y;
			}
			if (i == 0) leftbound_Y[j] = temp_Y;
			if (j == 0) bottombound_Y[i] = temp_Y;
		}
	}
	for(j=0;j<JMAX;j++){
		for(i=0;i<whole_IMAX;i++){
			myfile>>temp_Z;
			if (j>=mpi_j1[mpi_jindex]&&j<=mpi_j2[mpi_jindex]){
				jn=j-mpi_j1[mpi_jindex];
				whole_Z[i][jn]=temp_Z;
			}
			if (i==0)leftbound_Z[j]=temp_Z;
		}
	}
	myfile.close();
	whole_v=new double* [IMAX];
	whole_w=new double* [IMAX];
	whole_p=new double* [IMAX];
	for(i=0;i<IMAX;i++){
		whole_v[i]=new double [mpi_JMAX[mpi_jindex]];
		whole_w[i]=new double [mpi_JMAX[mpi_jindex]];
		whole_p[i]=new double [mpi_JMAX[mpi_jindex]];
	}

	// mpi: The specify sender and receiver matrix for pressure;
	pss=new double[mpi_IMAX[mpi_iindex]];psr=new double[mpi_IMAX[mpi_iindex]];
	pns=new double[mpi_IMAX[mpi_iindex]];pnr=new double[mpi_IMAX[mpi_iindex]];
	pes=new double[mpi_JMAX[mpi_jindex]];per=new double[mpi_JMAX[mpi_jindex]];
	pws=new double[mpi_JMAX[mpi_jindex]];pwr=new double[mpi_JMAX[mpi_jindex]];	
	// v,w velocity;
	vss=new double[mpi_IMAX[mpi_iindex]];vsr=new double[mpi_IMAX[mpi_iindex]];
	vns=new double[mpi_IMAX[mpi_iindex]];vnr=new double[mpi_IMAX[mpi_iindex]];
	ves=new double[mpi_JMAX[mpi_jindex]];ver=new double[mpi_JMAX[mpi_jindex]];
	vws=new double[mpi_JMAX[mpi_jindex]];vwr=new double[mpi_JMAX[mpi_jindex]];
	delvwss=new double[mpi_IMAX[mpi_iindex]];delvwsr=new double[mpi_IMAX[mpi_iindex]];
	delvwns=new double[mpi_IMAX[mpi_iindex]];delvwnr=new double[mpi_IMAX[mpi_iindex]];
	delvwes=new double[mpi_JMAX[mpi_jindex]];delvwer=new double[mpi_JMAX[mpi_jindex]];
	delvwws=new double[mpi_JMAX[mpi_jindex]];delvwwr=new double[mpi_JMAX[mpi_jindex]];	
	wss=new double[mpi_IMAX[mpi_iindex]];wsr=new double[mpi_IMAX[mpi_iindex]];
	wns=new double[mpi_IMAX[mpi_iindex]];wnr=new double[mpi_IMAX[mpi_iindex]];
	wes=new double[mpi_JMAX[mpi_jindex]];wer=new double[mpi_JMAX[mpi_jindex]];
	wws=new double[mpi_JMAX[mpi_jindex]];wwr=new double[mpi_JMAX[mpi_jindex]];	
	//mpi_fps;mpi_fpr;mpi_fvs;mpi_fvr;mpi_fws;mpi_fwr
	mpi_fs=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+1)];
	mpi_fr=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+1)];
	mpi_fr1=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+1)];
	mpi_fr2=new double[IMAX*(mpi_JMAX[mpi_zarea-1]+1)];

	//add the temperature distribution in the air
	speed_air=new double[mpi_JMAX[mpi_jindex]];
	temperature_air=new double[mpi_JMAX[mpi_jindex]];
	Aav=new double[JMAX];
	double slope_temperature;

	if (leftbound_Z[JMAX-1]<1.e-6){  // temperature at the ground
		for (j=0;j<mpi_JMAX[mpi_jindex];j++){
			temperature_air[j]=AirPara.temperature1;
			speed_air[j]=331.3*sqrt(1+temperature_air[j]/273.15);
			Aav[j]=speed_air[j]*speed_air[j]/AirPara.adiabatic_coef/AirPara.Pav;
		}
	}else{
		
		slope_temperature=(AirPara.temperature2-AirPara.temperature1)/fabs(leftbound_Z[JMAX-1]-leftbound_Z[0]);
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			jn=mpi_j1[mpi_jindex]+j;
			temperature_air[j]=slope_temperature*(leftbound_Z[jn]-leftbound_Z[0])+AirPara.temperature1;	
			speed_air[j]=331.3*sqrt(1+temperature_air[j]/273.15);
			Aav[j]=speed_air[j]*speed_air[j]/AirPara.adiabatic_coef/AirPara.Pav;
		}
		/*speed_air1=331.3*sqrt(1+AirPara.temperature1/273.15);
		speed_air2=331.3*sqrt(1+AirPara.temperature2/273.15);
		slope_sound=(speed_air2-speed_air1)/fabs(leftbound_Z[JMAX-1]-leftbound_Z[0]);
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			jn=mpi_j1[mpi_jindex]+j;
			speed_air=slope_sound*(leftbound_Z[jn]-leftbound_Z[0])+speed_air1;	
			Aav[j]=speed_air*speed_air/AirPara.adiabatic_coef/AirPara.Pav;
		}*/
	}	
	Y=new double* [mpi_IMAX[mpi_iindex]];Z=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		Y[i]=new double [mpi_JMAX[mpi_jindex]];Z[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	Y_over_y=new double* [mpi_IMAX[mpi_iindex]];Y_over_z=new double* [mpi_IMAX[mpi_iindex]];
	Z_over_y=new double* [mpi_IMAX[mpi_iindex]];Z_over_z=new double* [mpi_IMAX[mpi_iindex]];
	Jacobian=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		Y_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Y_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
		Z_over_y[i]=new double [mpi_JMAX[mpi_jindex]];Z_over_z[i]=new double [mpi_JMAX[mpi_jindex]];
		Jacobian[i]=new double [mpi_JMAX[mpi_jindex]];
	}
	//allocate memory for velocity and pressure
	v_n=new double* [mpi_IMAX[mpi_iindex]];w_n=new double* [mpi_IMAX[mpi_iindex]];
	p_n=new double* [mpi_IMAX[mpi_iindex]];
	v_nn=new double* [mpi_IMAX[mpi_iindex]];w_nn=new double* [mpi_IMAX[mpi_iindex]];
	p_nn=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		v_n[i]=new double [mpi_JMAX[mpi_jindex]];w_n[i]=new double [mpi_JMAX[mpi_jindex]];
		p_n[i]=new double [mpi_JMAX[mpi_jindex]];
		v_nn[i]=new double [mpi_JMAX[mpi_jindex]];w_nn[i]=new double [mpi_JMAX[mpi_jindex]];
		p_nn[i]=new double [mpi_JMAX[mpi_jindex]];
	}

	// add the timeseries data
	jjmax1=501;
	jjmax2=701;
	judge_left=1;
	if (judge_left==1){
		if (whole_JMAX==1201){
			jjmax1=0;
			jjmax2=1201;
		}
		if(whole_JMAX==401){
			jjmax1=0;
			jjmax2=400;
		}
		
		}
	delvw_nn=new double* [mpi_IMAX[mpi_iindex]];
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		delvw_nn[i]=new double[mpi_JMAX[mpi_jindex]];
	}

	swp_nn = new double*[IMAX]; swvw_nn = new double*[IMAX];
	for (i = 0; i < IMAX; i++) {
		swp_nn[i] = new double[JMAX]; swvw_nn[i] = new double[JMAX];
	}
	for (i = 0; i < IMAX; i++) {
		for (j = 0; j < JMAX; j++) {
			swp_nn[i][j] = 0.0; swvw_nn[i][j] = 0.0;
		}
	}
	cout << "rank = " << mpi_rank << ", mpi_JMAX = " << mpi_JMAX[mpi_jindex] << ", mpi_j1 = " << mpi_j1[mpi_jindex] << ", mpi_j2 = " << mpi_j2[mpi_jindex] << endl;
}

calculation_2D::~calculation_2D()
{
	int i;
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		delete[] p_n[i];delete[] w_n[i];delete[] v_n[i];
		delete[] p_nn[i];delete[] w_nn[i];delete[] v_nn[i];
		delete[] Y_over_y[i];delete[] Y_over_z[i];
		delete[] Z_over_y[i];delete[] Z_over_z[i];
		delete[] Y[i];delete[] Z[i];delete[] Jacobian[i];delete[] delvw_nn[i];
	}
	delete[] v_n;delete[] w_n;delete[] p_n;
	delete[] v_nn;delete[] w_nn;delete[] p_nn;
	delete[] Y_over_y;delete[] Y_over_z;
	delete[] Z_over_y;delete[] Z_over_z;
	delete[] Y;delete[] Z;delete[] Jacobian;delete[] delvw_nn;
	for(i=0;i<IMAX;i++){
		delete[] whole_v[i];delete[] whole_w[i];delete[] whole_p[i];
		delete[] whole_Y[i];delete[] whole_Z[i];
		delete[] swp_nn[i]; delete[] swvw_nn[i];
	}
	delete[] whole_v;delete[] whole_w;delete[] whole_p;
	delete[] whole_Y;delete[] whole_Z;
	delete[] swp_nn; delete[] swvw_nn;
	delete[] speed_air;delete[] temperature_air; delete[] Aav;
	delete[] time_step; delete[] init_pressure;
	delete[] pss; delete[] psr; delete[] pns; delete[] pnr;
	delete[] pes; delete[] per; delete[] pws; delete[] pwr;
	delete[] vss; delete[] vsr; delete[] vns; delete[] vnr;
	delete[] ves; delete[] ver; delete[] vws; delete[] vwr;
	delete[] wss; delete[] wsr; delete[] wns; delete[] wnr;
	delete[] wes; delete[] wer; delete[] wws; delete[] wwr;
	delete[] delvwss; delete[] delvwsr; delete[] delvwns; delete[] delvwnr;
	delete[] delvwes; delete[] delvwer; delete[] delvwws; delete[] delvwwr;
	delete[] mpi_fs;delete[] mpi_fr;delete[] mpi_fr1;delete[] mpi_fr2;
	delete[] mpi_i1;delete[] mpi_i2;delete[] mpi_j1;
	delete[] mpi_j2;delete[] mpi_IMAX;delete[] mpi_JMAX;
	delete[] leftbound_Y;delete[] leftbound_Z;
	delete[] bottombound_Y;
}
double calculation_2D::cal_InitialPressure(double y,double z)
{
	double r;
	r=sqrt(pow(y,2)+pow(z,2));
	return (1.0*exp(-gaussian_coef*r*r));//this is from Salomons paper(2002) p=exp(-40*r*r)
}

void calculation_2D::Set_InitialCond(Position source)
{
		
	int i,j;
	for(j=0;j<mpi_JMAX[mpi_jindex];j++){
		for(i=0;i<whole_IMAX;i++){
			whole_v[i][j]=0.0;
			whole_w[i][j]=0.0;
			if (judge_left==1) whole_p[i][j]=0.0;
			else {
				whole_p[i][j]=cal_InitialPressure(whole_Y[i][j]-source.y,
				whole_Z[i][j]-source.z);
			}
		}
	}
}


void calculation_2D::get_position(int &ii,int &jj,Position receiver)
{	
	int i, j;

	for (i=0;i<whole_IJKMAX;i++){
		for (j=0;j<JMAX;j++){
			if (((receiver.y-(leftbound_Y[j]+i*diff_y))<(diff_y/8.0))&&
				((receiver.z-leftbound_Z[j])<(diff_z/8.0))){
				ii=i;jj=j;
				return;
			}
		}
	}

// for fixed domain we can use the fllowing method
/*		ii=0;
	jj=0;
	for (j=0;j<JMAX;j++){
		if((receiver.z-leftbound_Z[j])<(diff_z/8.0)) {
			jj=j;
			break;
		}
	}
	for (i=0;i<IMAX;i++){
		if ((receiver.y-bottombound_Y[i])<(diff_y/8.0)){
			ii=i;
			return;
		}
	}*/
}

void calculation_2D::UpdateInitialCond(int N)
{
	//transfer data under time n+1 to data under time n
	int i,j,in,jn,tag,rank_send,rank_recv,rank_num;
	int i1,i2,mpi_index;//j2,j1,
	for(i=0;i<mpi_IMAX[mpi_iindex];i++){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			v_n[i][j]=v_nn[i][j];
			w_n[i][j]=w_nn[i][j];
			p_n[i][j]=p_nn[i][j];
		}
	}
	//update v,w,p for whole domain
	if (mpi_iindex!=0){
		rank_num=mpi_IMAX[mpi_iindex]*mpi_JMAX[mpi_jindex];
		rank_recv=mpi_jindex;
		if (mpi_iindex!=0)i1=1;
		else i1=0;
		if (mpi_iindex!=mpi_yarea-1)i2=1;
		else i2=0;
		if (N!=0){
			mpi_index=0;
			for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
				for(j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=Y[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=16;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

			mpi_index=0;
			for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
				for(j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=Z[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=17;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				mpi_fs[mpi_index]=v_nn[i][j];
				mpi_index=mpi_index+1;
			}
		}
		tag=18;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		
		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				mpi_fs[mpi_index]=w_nn[i][j];
				mpi_index=mpi_index+1;
			}
		}
		tag=19;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		
		mpi_index=0;
		for(i=i1;i<mpi_IMAX[mpi_iindex]-i2;i++){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){
				mpi_fs[mpi_index]=p_nn[i][j];
				mpi_index=mpi_index+1;
			}
		}
		tag=20;
		MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		
		rank_send=rank_recv;
		rank_num=IMAX*mpi_JMAX[mpi_jindex];
		tag=21;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				whole_v[i][j]=mpi_fr[mpi_index];
				mpi_index=mpi_index+1;
			}
		}

		tag=22;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				whole_w[i][j]=mpi_fr[mpi_index];
				mpi_index=mpi_index+1;
			}
		}

		tag=23;
		MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		mpi_index=0;
		for (i=0;i<IMAX;i++){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++){
				whole_p[i][j]=mpi_fr[mpi_index];
				mpi_index=mpi_index+1;
			}
		}
	}else{
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				whole_v[in][j]=v_nn[i][j];
				whole_w[in][j]=w_nn[i][j];
				whole_p[in][j]=p_nn[i][j];
				if (N!=0){
					whole_Y[in][j]=Y[i][j];
					whole_Z[in][j]=Z[i][j];
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
			if (N!=0){
				mpi_index=0;
				tag=16;
				MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
				for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
					for(j=0;j<mpi_JMAX[jn];j++){
						whole_Y[i][j]=mpi_fr[mpi_index];
						mpi_index=mpi_index+1;
					}
				}
				mpi_index=0;
				tag=17;
				MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
				for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
					for(j=0;j<mpi_JMAX[jn];j++){
						whole_Z[i][j]=mpi_fr[mpi_index];
						mpi_index=mpi_index+1;
					}
				}
			}
			mpi_index=0;
			tag=18;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					whole_v[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
			}
			mpi_index=0;
			tag=19;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					whole_w[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
			}
			mpi_index=0;
			tag=20;
			MPI_Recv(mpi_fr, rank_num, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			for(i=mpi_i1[in]+i1;i<mpi_i2[in]+1-i2;i++){
				for(j=0;j<mpi_JMAX[jn];j++){
					whole_p[i][j]=mpi_fr[mpi_index];
					mpi_index=mpi_index+1;
				}
			}
		}
		for (in=1;in<mpi_yarea;in++){
			jn=mpi_jindex;
			rank_num=IMAX*mpi_JMAX[jn];
			if (mpi_porous!=0)rank_recv=(mpi_iindex+in)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex+in)*mpi_zarea+mpi_jindex;
			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=whole_v[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=21;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			
			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=whole_w[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=22;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

			mpi_index=0;
			for (i=0;i<IMAX;i++){
				for (j=0;j<mpi_JMAX[mpi_jindex];j++){
					mpi_fs[mpi_index]=whole_p[i][j];
					mpi_index=mpi_index+1;
				}
			}
			tag=23;
			MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
}


double calculation_2D::get_pressure(int ii,int jj,int N)
{
	int i,il,in,jn,rank_send,rank_recv,tag;
	double p_rec;
	rank_send = 0;
	rank_recv = 0;
	if ((ii>N*(move_frame.trail_DI-1))&&(ii<IMAX-1+N*(move_frame.trail_DI-1))){
		il=ii-N*(move_frame.trail_DI-1);
		if(il>mpi_i1[mpi_iindex]&& il<mpi_i2[mpi_iindex]&& 
			jj>mpi_j1[mpi_jindex]&&jj<mpi_j2[mpi_jindex]){		
			in=il-mpi_i1[mpi_iindex];
			jn=jj-mpi_j1[mpi_jindex];
			if (mpi_rank!=0){
				tag=0;
				rank_recv=0;
				MPI_Send(&p_n[in][jn], 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
				return 0.0;
			}else return p_n[in][jn];
		}else if (mpi_rank==0){
			for (i=1;i<mpi_size;i++){
				in=i/mpi_zarea;
				jn=i%mpi_zarea;
				if(il>mpi_i1[in]&& il<mpi_i2[in]&&
                			jj>mpi_j1[jn]&&jj<mpi_j2[jn]){
					rank_send=i;
					i=mpi_size;
				}
			}
			tag=0;
                        MPI_Recv(&p_rec,1,MPI_DOUBLE,rank_send,tag,  MPI_COMM_WORLD, &status);
			return p_rec;
		}else return 0.0;
	}else return 0.0;
}

void calculation_2D::SetMovingFrame(int N)
//N is total frame number starting from 0
{
	int mpi_ny,i;
	if(N==0){
		IMAX=whole_IMAX;
	}else{
		IMAX=move_frame.IMAX+N*(move_frame.lead_DI-move_frame.trail_DI);
	}
// Redefine the number of grids in y drection 
	mpi_ny=int(IMAX/mpi_yarea);
	for (i=0;i<mpi_yarea;i++){
		mpi_i1[i]=i*mpi_ny;
		if (i!=mpi_yarea-1)mpi_i2[i]=(i+1)*mpi_ny+1;
		else mpi_i2[i]=IMAX-1;
		mpi_IMAX[i]=mpi_i2[i]-mpi_i1[i]+1;
	};
	SetLeftRightMoving(N);

}

void calculation_2D::SetLeftRightMoving(int N)
{
	int i,j,in,jn;
	//transfer coordinate
	if (N==0){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				Y[i][j]=whole_Y[in][j];
				Z[i][j]=whole_Z[in][j];
			}
		}
		Y_start = bottombound_Y[0];
	}else{
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				jn=mpi_j1[mpi_jindex]+j;
				Y[i][j] = bottombound_Y[in] + N*(move_frame.trail_DI - 1)*diff_y;
				Z[i][j]=leftbound_Z[jn];
			}
		}
		Y_start = bottombound_Y[0] + N*(move_frame.trail_DI - 1)*diff_y;
	}

	//calculate partial differential of coordinate transformation using central difference
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
				Y_over_y[i][j]=(Y[i+1][j]-Y[i-1][j])/diff_y/2;
				Z_over_y[i][j]=(Z[i+1][j]-Z[i-1][j])/diff_y/2;
				Y_over_z[i][j]=(Y[i][j+1]-Y[i][j-1])/diff_z/2;
				Z_over_z[i][j]=(Z[i][j+1]-Z[i][j-1])/diff_z/2;
				Jacobian[i][j]=Y_over_y[i][j]*Z_over_z[i][j]-Y_over_z[i][j]*Z_over_y[i][j];
		}
	}
	//update v,w,p for moving frame
	int temp_I;
	if (N==0){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				v_n[i][j]=whole_v[in][j];
				w_n[i][j]=whole_w[in][j];
				p_n[i][j]=whole_p[in][j];
			}
		}
		temp_I=0;
   	 }else temp_I=move_frame.IMAX-move_frame.lead_DI;
	// temp_I is the overlapping part after frame is moved in y direction
	if (temp_I!=0){
		for(j=0;j<mpi_JMAX[mpi_jindex];j++){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				in=mpi_i1[mpi_iindex]+i;
				if (in<temp_I){
					v_n[i][j]=whole_v[in+move_frame.trail_DI-1][j];
					w_n[i][j]=whole_w[in+move_frame.trail_DI-1][j];
					p_n[i][j]=whole_p[in+move_frame.trail_DI-1][j];
				}else{
					v_n[i][j]=0.0;
					w_n[i][j]=0.0;
					p_n[i][j]=0.0;								
				}
			}
		}
		// specfy the left boundary condition which is used to avoid numerical error in the case wind shear. 
		if (mpi_iindex==0){
			for (j=0;j<mpi_JMAX[mpi_jindex];j++)p_n[0][j]=p_n[1][j];
		}
	}
}
void calculation_2D::save_restart_cal(char *restartfile)
{
	int i,j;
	ofstream outfile11(restartfile,ios::app|ios::binary);
	outfile11.setf(ios::scientific,ios::floatfield);
	outfile11.setf(ios::showpos);
	outfile11.precision(6);
	outfile11.width(14);
	for(j=0;j<mpi_JMAX[mpi_jindex];j++){
		for(i=0;i<IMAX;i++){
			outfile11<<whole_Y[i][j]<<" ";
			outfile11<<whole_Z[i][j]<<" ";
			outfile11<<whole_v[i][j]<<" ";
			outfile11<<whole_w[i][j]<<" ";
			outfile11<<whole_p[i][j]<<" ";
		}
	}
	outfile11<<endl;
	outfile11.close();
}

void calculation_2D::input_restart_cal(char *restart_infile,int *point_pois,int N)
{
	int i,j,int_pois;
	N=N-1;
	if(N==(int)((whole_IJKMAX-move_frame.IMAX)/(move_frame.lead_DI-1)+1.0e-6)){  
		IMAX=whole_IJKMAX-N*(move_frame.trail_DI-1);
	}else{
		IMAX=move_frame.IMAX+N*(move_frame.lead_DI-move_frame.trail_DI);
	}
	int_pois=*point_pois;
	ifstream infile(restart_infile,ios::in|ios::binary);
	infile.seekg(int_pois);
	for (j=0;j<mpi_JMAX[mpi_jindex];j++){
		for (i=0;i<IMAX;i++){
			infile>>whole_Y[i][j];
			infile>>whole_Z[i][j];
			infile>>whole_v[i][j];
			infile>>whole_w[i][j];
			infile>>whole_p[i][j];
		}
	}
	infile.ignore(100,'\n');
	int_pois=infile.tellg();
	*point_pois=int_pois;
	infile.close();
}

//-----------------------------update boundary condtions of pressure--------------------//
void calculation_2D::UpdateBC_pressure(boundary_location BC,int time_judge,int time_current)
{
	int i,j,in,jn,rank_send,rank_recv,tag;
	double time_length;//,k_slope
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

void calculation_2D::UpdateBC_delvw(boundary_location BC)
{
	int i,j,in,jn,rank_send,rank_recv,tag;
	for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		if (mpi_jindex!=0) {
			delvwss[i]=delvw_nn[i][1];
		}
		if (mpi_jindex!=mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-2;
			delvwns[i]=delvw_nn[i][jn];
		}
	}
	for (j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
		if (mpi_iindex!=0){
			delvwws[j]=delvw_nn[1][j];
		}
		if (mpi_iindex!=mpi_yarea-1){
			in=mpi_IMAX[mpi_iindex]-2;
			delvwes[j]=delvw_nn[in][j];
		}
	
	}
	switch(BC){
	case WestBC:
	{
		if (mpi_iindex==0){ 
			for(j=0;j<mpi_JMAX[mpi_jindex];j++) {
				delvw_nn[0][j]=delvw_nn[1][j];////left
			}
			
		}else{
			tag=31;
			if (mpi_porous!=0)rank_recv=(mpi_iindex-1)*mpi_porous+mpi_jindex;
			else rank_recv=(mpi_iindex-1)*mpi_zarea+mpi_jindex;
			MPI_Send(delvwws, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=32;
			rank_send=rank_recv;
			MPI_Recv(delvwwr, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		}

	}
	break;
	case EastBC:
	{
		if (mpi_iindex==mpi_yarea-1){
			for(j=0;j<mpi_JMAX[mpi_jindex];j++){ 
				in=mpi_IMAX[mpi_iindex]-1;
				delvw_nn[in][j]=delvw_nn[in-1][j];//right
			}
		}else{
			tag=31;
			if (mpi_porous!=0)rank_send=(mpi_iindex+1)*mpi_porous+mpi_jindex;
			else rank_send=(mpi_iindex+1)*mpi_zarea+mpi_jindex;
			MPI_Recv(delvwer, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=32;
			rank_recv=rank_send;
			MPI_Send(delvwes, mpi_JMAX[mpi_jindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
	break;
	case SouthBC:
	{
		if (mpi_jindex==0){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++)delvw_nn[i][0]=delvw_nn[i][1]; //Bottom
		}else{
			tag=33;
			if (mpi_porous!=0)rank_recv=mpi_iindex*mpi_porous+mpi_jindex-1;
			else rank_recv=mpi_iindex*mpi_zarea+mpi_jindex-1;
			MPI_Send(delvwss, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
			tag=34;
			rank_send=rank_recv;
			MPI_Recv(delvwsr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
		}
	}
	break;
	case NorthBC://characteristic boundary condition
	{
		if (mpi_jindex==mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			for(i=0;i<mpi_IMAX[mpi_iindex];i++) delvw_nn[i][jn]=delvw_nn[i][jn-1];//upper
		}else{
			tag=33;
			if (mpi_porous!=0)rank_send=mpi_iindex*mpi_porous+mpi_jindex+1;
			else rank_send=mpi_iindex*mpi_zarea+mpi_jindex+1;
			MPI_Recv(delvwnr, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_send, tag, MPI_COMM_WORLD, &status);
			tag=34;
			rank_recv=rank_send;
			MPI_Send(delvwns, mpi_IMAX[mpi_iindex] - 1, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
		}
	}
	break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
	
	
	for (i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
		if(mpi_jindex!=0) delvw_nn[i][0]=delvwsr[i];
		if(mpi_jindex!=mpi_zarea-1){
			jn=mpi_JMAX[mpi_jindex]-1;
			delvw_nn[i][jn]=delvwnr[i];
		}
	}
	for (j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
		if (mpi_iindex!=0)delvw_nn[0][j]=delvwwr[j];
		if (mpi_iindex!=mpi_yarea-1){
			in=mpi_IMAX[mpi_iindex]-1;
			delvw_nn[in][j]=delvwer[j];
		}	
	}
}

//--------------------------update boundary conditions of velocity(V and W)------------------//
void calculation_2D::UpdateBC_velocity(boundary_location BC)
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
				w_nn[i][0]=0.0;v_nn[i][0]=v_nn[i][1];
			} //Bottom
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

void calculation_2D::UpdateBC_velocity(boundary_location BC,calculation_2D& porous)
{
	//update absorbing boundary
	int AirDy_CouplingDy,AirDz_CouplingDz;
	AirDy_CouplingDy=(int)(diff_y/porous.diff_y);
	AirDz_CouplingDz=(int)(diff_z/porous.diff_z);
	int i;
	switch(BC){
	case SouthBC:
		if(AirDz_CouplingDz<1.5){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				v_nn[i][0]=porous.v_nn[i][porous.mpi_JMAX[mpi_jindex]-2];
				porous.v_nn[i][porous.mpi_JMAX[mpi_jindex]-1]=v_nn[i][1];
				w_nn[i][0]=porous.w_nn[i][porous.mpi_JMAX[mpi_jindex]-2];
				porous.w_nn[i][porous.mpi_JMAX[mpi_jindex]-1]=w_nn[i][1];
			}
		}else{//transient interface
			for(i=0;i<(mpi_IMAX[mpi_iindex]-1);i++){
				for(int k=0;k<(AirDy_CouplingDy);k++){
					double temp_v;
					temp_v=(v_nn[i+1][1]-v_nn[i][1])/(AirDy_CouplingDy)*k+v_nn[i][1];
					porous.v_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-1]=
						(temp_v-porous.v_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2])/
						(AirDz_CouplingDz+1)+porous.v_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2];
					double temp_w;
					temp_w=(w_nn[i+1][1]-w_nn[i][1])/(AirDy_CouplingDy)*k+w_nn[i][1];
					porous.w_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-1]=
						(temp_w-porous.w_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2])/
						(AirDz_CouplingDz+1)+porous.w_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2];
				}
			}
			for(i=1;i<(mpi_IMAX[mpi_iindex]-1);i++){
				v_nn[i][0]=porous.v_nn[i*AirDy_CouplingDy][porous.mpi_JMAX[mpi_jindex]-1];
				w_nn[i][0]=porous.w_nn[i*AirDy_CouplingDy][porous.mpi_JMAX[mpi_jindex]-1];
			}
		}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}

void calculation_2D::UpdateBC_pressure(boundary_location BC,calculation_2D& porous)
{
	//update absorbing boundary
	int AirDy_CouplingDy,AirDz_CouplingDz;
	AirDy_CouplingDy=(int)(diff_y/porous.diff_y);
	AirDz_CouplingDz=(int)(diff_z/porous.diff_z);
	int i;
	switch(BC){
	case SouthBC:
	{	
		if(AirDz_CouplingDz<1.5){
			for(i=0;i<mpi_IMAX[mpi_iindex];i++){
				p_nn[i][0]=porous.p_nn[i][porous.mpi_JMAX[mpi_jindex]-2];
				porous.p_nn[i][porous.mpi_JMAX[mpi_jindex]-1]=p_nn[i][1];
			}
		}else{//transient interface
			for(i=0;i<mpi_IMAX[mpi_iindex]-1;i++){
				for(int k=0;k<AirDy_CouplingDy;k++){
					double temp_p;
					temp_p=(p_nn[i+1][1]-p_nn[i][1])/(AirDy_CouplingDy)*k+p_nn[i][1];
					porous.p_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-1]=
						(temp_p-porous.p_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2])/
						(AirDz_CouplingDz+1)+porous.p_nn[i*AirDy_CouplingDy+k][porous.mpi_JMAX[mpi_jindex]-2];
				}
			}
			for(i=1;i<(mpi_IMAX[mpi_iindex]-1);i++){
				p_nn[i][0]=porous.p_nn[i*AirDy_CouplingDy][porous.mpi_JMAX[mpi_jindex]-1];
			}
		}
	}
		break;
	default:
		cout<<"it is wrong to entrance boundary";
	}
}


//cal_fvw() and cal_fp() to get v_nn=fvw*diff_t;p_nn=fp*diff_t//
void calculation_2D::Cal_velocity(void)
{
	int i,j;
	switch(scheme){
	case FB_v_p://forward - backward scheme (1st)
	{	
		cal_fvw(v_n,w_n,p_n);
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
				v_nn[i][j]=v_n[i][j]+v_nn[i][j];
				w_nn[i][j]=w_n[i][j]+w_nn[i][j];
			}
		}
	}
		break;
	case FB_p_v://forward - backward scheme (1st)

		cal_fvw(v_n,w_n,p_nn);
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
				v_nn[i][j]=v_n[i][j]+v_nn[i][j];
				w_nn[i][j]=w_n[i][j]+w_nn[i][j];
			}
		}
		break;
	default:
		printf("wrong with scheme name");
	}//end switch
}

void calculation_2D::Cal_pressure(void)
{
	int i,j;
	//double temp_p;
	switch(scheme){
	case FB_v_p://forward - backward scheme (1st)
	{
		cal_delvw();
		UpdateBC_delvw(NorthBC);
		UpdateBC_delvw(SouthBC);
		UpdateBC_delvw(WestBC);
		UpdateBC_delvw(EastBC);
		fractional_CD();
		Update_FCD();
		//shift_grunwald();
		cal_fp(v_nn,w_nn,p_n);
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
				p_nn[i][j]=p_n[i][j]+p_nn[i][j];
			}
		}
	}
	break;
	case FB_p_v://forward - backward scheme (1st)
	{
		cal_fp(v_n,w_n,p_n);
		UpdateBC_delvw(NorthBC);
		UpdateBC_delvw(SouthBC);
		UpdateBC_delvw(WestBC);
		UpdateBC_delvw(EastBC);
		shift_grunwald();
		for(j=1;j<mpi_JMAX[mpi_jindex]-1;j++){
			for(i=1;i<mpi_IMAX[mpi_iindex]-1;i++){
				p_nn[i][j]=p_n[i][j]+p_nn[i][j];
			}
		}
	}
	break;
	default:
		printf("wrong with scheme name");
	}//end switch
}

void calculation_2D::Update_FCD(void)
{
	for (int rank = 0; rank < mpi_zarea; rank++) {
		int j1 = mpi_j1[rank] + 1;
		int j2 = mpi_j2[rank] - 1;
		int size = (mpi_IMAX[0] - 2) * (mpi_JMAX[rank] - 2);
		double *send_buf = new double[size];
		double *rec_buf = new double[size];
		for (int i = 0; i < size; i++) { send_buf[i] = 0; rec_buf[i] = 0; }
		//cout << endl;
		Chg_2D_to_1D(j1, j2, send_buf, swp_nn);
		MPI_Reduce(send_buf, rec_buf, size, MPI_DOUBLE, MPI_SUM, rank, MPI_COMM_WORLD);
		Chg_1D_to_2D(j1, j2, rec_buf, swp_nn);

		for (int i = 0; i < size; i++) { send_buf[i] = 0; rec_buf[i] = 0; }
		Chg_2D_to_1D(j1, j2, send_buf, swvw_nn);
		MPI_Reduce(send_buf, rec_buf, size, MPI_DOUBLE, MPI_SUM, rank, MPI_COMM_WORLD);
		Chg_1D_to_2D(j1, j2, rec_buf, swvw_nn);
		delete[] send_buf; delete[] rec_buf;
	}
}

void calculation_2D::Chg_2D_to_1D(int j1, int j2, double * send_buf, double ** u)
{
	int rec = 0;
	for (int i = 1; i < mpi_IMAX[mpi_iindex] - 1; i++) {
		for (int j = j1; j <= j2; j++) {
			send_buf[rec] = u[i][j];
			rec++;
		}
	}
}

void calculation_2D::Chg_1D_to_2D(int j1, int j2, double * rec_buf, double ** u)
{
	int rec = 0;
	for (int i = 1; i < mpi_IMAX[mpi_iindex] - 1; i++) {
		for (int j = j1; j <= j2; j++) {
			u[i][j] = rec_buf[rec];
			rec++;
		}
	}
}

void calculation_2D::mpi_send_data(void)
{
	int i,j,j1,j2;
	int rank_num,mpi_index,tag,rank_recv;
	mpi_index=0;
	rank_recv=0;
	rank_num=IMAX*mpi_JMAX[mpi_jindex];
	if (mpi_jindex!=0)j1=1;
	else j1=0;
	if (mpi_jindex!=mpi_zarea-1)j2=1;
	else j2=0;
	for (j=j1;j<mpi_JMAX[mpi_jindex]-j2;j++){
		for (i=0;i<IMAX;i++){
			if (i%3==0||i==IMAX-1){
				mpi_fs[mpi_index]=whole_Y[i][j];
				mpi_index=mpi_index+1;
			}
		}
	}
	tag=13;
	MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

	mpi_index=0;
	for (j=j1;j<mpi_JMAX[mpi_jindex]-j2;j++){
		for (i=0;i<IMAX;i++){
			if (i%3==0||i==IMAX-1){
				mpi_fs[mpi_index]=whole_Z[i][j];
				mpi_index=mpi_index+1;
			}
        }
    }
	tag=14;
	MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);

	mpi_index=0;
	for (j=j1;j<mpi_JMAX[mpi_jindex]-j2;j++){
		for (i=0;i<IMAX;i++){
			if (i%3==0||i==IMAX-1){
				mpi_fs[mpi_index]=whole_p[i][j];
				mpi_index=mpi_index+1;
			}
        }
    }
	tag=15;
	MPI_Send(mpi_fs, rank_num, MPI_DOUBLE, rank_recv, tag, MPI_COMM_WORLD);
}
void getOstream(ostream& output_stream, calculation_2D& my_object)
{
	int IMAX1,IMAX,JMAX,i,j,jn,j1,j2;
	int mpi_index,rank_num,rank_send,tag;
	IMAX=my_object.IMAX;JMAX=my_object.JMAX;
	if ((IMAX-1)%3==0)IMAX1=int((IMAX-1+1.0e-6)/3.0)+1;
	else IMAX1=int((IMAX-1+1.0e-6)/3.0)+2;
	output_stream<<"I="<<IMAX1<<",J="<<JMAX<<",F=POINT"<<endl;
	if (my_object.mpi_zarea-1!=0)j2=1;
	else j2=0;
	for(j=0;j<my_object.mpi_JMAX[0]-j2;j++){
		for(i=0;i<IMAX;i++){
			if (i%3==0||i==IMAX-1){
				output_stream<<(my_object.whole_Y[i][j])<<" "<<(my_object.whole_Z[i][j])<<" "
					<<my_object.whole_p[i][j]<<" "<<endl;
			}
		}
	}
	for (rank_send=1;rank_send<my_object.mpi_zarea;rank_send++){
		jn=my_object.mpi_JMAX[rank_send];
		if (rank_send!=0)j1=1;
		else j1=0;
		if (rank_send!=my_object.mpi_zarea-1)j2=1;
		else j2=0;
		mpi_index=0;
		rank_num=IMAX*jn;
		tag=13;
		MPI_Recv(my_object.mpi_fr,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&my_object.status);
		tag=14;
		MPI_Recv(my_object.mpi_fr1,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&my_object.status);
		tag=15;
		MPI_Recv(my_object.mpi_fr2,rank_num,MPI_DOUBLE,rank_send,tag, MPI_COMM_WORLD,&my_object.status);
		for(j=j1;j<jn-j2;j++){
			for(i=0;i<IMAX;i++){
				if (i%3==0||i==IMAX-1){
					output_stream<<(my_object.mpi_fr[mpi_index])<<" "<<(my_object.mpi_fr1[mpi_index])<<" "
							<<my_object.mpi_fr2[mpi_index]<<" "<<endl;
					mpi_index=mpi_index+1;
				}
			}
		}
	}
}

ostream& operator<<(ostream& output_stream, calculation_2D& my_object)
{
	getOstream(output_stream,my_object);
	return output_stream;
}
