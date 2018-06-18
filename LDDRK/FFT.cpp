#include <math.h>
#include <stdio.h>
#include "mygeneral.h"

/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
void FFT(short int dir,long m,double x[],double y[])
{
   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
	return ;
}

void FFT_output(long m,double x[],double y[],char* filename,double diff_t)
{
	Complex comp;
	int n;
	/* Calculate the number of points */
   n = 1;
   for (int i=0;i<m;i++) 
      n *= 2;
	//output to file
   FILE *fp_p_f;
	if((fopen_s(&fp_p_f,filename,"w"))== NULL){
		printf( "can't creat or write file\n" );
		return;
	}
	//fprintf(fp_p_f,"ZONE I=%d,F=POINT \n",n/2);
	comp=Complex(x[0],y[0]);
	fprintf(fp_p_f,"%10.10e %10.10e \n",0.0,abs(comp));//,arg(comp)
	for(int jj=1;jj<=(n/2-1);jj++){
		comp=Complex(x[jj],y[jj]);
		fprintf(fp_p_f,"%10.10e %10.10e \n",jj/diff_t/n,abs(comp));//,arg(comp)
	}
	comp=Complex(x[n/2],y[n/2]);
	fprintf(fp_p_f,"%10.10e %10.10e \n",1.0/2/diff_t,abs(comp));//,arg(comp)
	//fprintf(fp_p_f,"%10.10f %10.10f \n",-1.0/2/diff_t,sqrt(x[n/2]*x[n/2]+y[n/2]*y[n/2]));
	//for(n=(FFT_N/2+1);n<=(FFT_N-1);n++){
	//	fprintf(fp_p_f,"%10.10f %10.10f \n",(n-FFT_N)/diff_t_air/FFT_N,sqrt(pr[n]*pr[n]+pi[n]*pi[n]));
	//}
	fclose(fp_p_f);
}

//interpolation of Aitken method
double enatk(double *x,double *y,int n,double t)
{ 
	double eps=1.0e-6;
	int i,j,k,m,l;
    double z,xx[10],yy[10];
    z=0.0;
    if (n<1) return(z);
    if (n==1) { z=y[0]; return(z);}
    m=10;
    if (m>n) m=n;
    if (t<=x[0]) k=1;
    else if (t>=x[n-1]) k=n;
    else
      { k=1; j=n;
        while ((k-j!=1)&&(k-j!=-1))
          { l=(k+j)/2;
            if (t<x[l-1]) j=l;
            else k=l;
          }
        if (fabs(t-x[l-1])>fabs(t-x[j-1])) k=j;
      }
    j=1; l=0;
    for (i=1;i<=m;i++)
      { k=k+j*l;
        if ((k<1)||(k>n))
          { l=l+1; j=-j; k=k+j*l;}
        xx[i-1]=x[k-1]; yy[i-1]=y[k-1];
        l=l+1; j=-j;
      }
    i=0;
    do
      { i=i+1; z=yy[i];
        for (j=0;j<=i-1;j++)
          z=yy[j]+(t-xx[j])*(yy[j]-z)/(xx[j]-xx[i]);
        yy[i]=z;
      }
    while ((i!=m-1)&&(fabs(yy[i]-yy[i-1])>eps));
    return z;
}

