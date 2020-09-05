#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "groups.h"

#define HALO_MAX 1.0E+16

float g7_ngal,
  g7_msub;

/*external functions
 */
float qromo(float (*func)(float), float a, float b,
	     float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);
float zbrent(float (*func)(float), float x1, float x2, float tol);

/* Local functions
 */
float func_ngal(float mag);
float func_nhalo(float m);
float subhalo_abundance(float m);
float subhalo_abundance2(float m);
float halo_abundance(float m);
float halo_abundance2(float m);
float func_subhalo1(float mhost);
float func_match_nhalo(float mmin);
float func_match_ngal(float mass);
float func_match_nhost(float mass);

float density2halo(float galaxy_density)
{
  g7_ngal = galaxy_density;
  return exp(zbrent(func_match_ngal,log(1.0E+7),log(HALO_MAX),1.0E-5));
  
}
float density2host_halo(float galaxy_density)
{
  g7_ngal = galaxy_density;
  return exp(zbrent(func_match_nhost,log(1.0E+7),log(HALO_MAX),1.0E-5));  
}

/* keep a running tab of calls to this function 
 * (until reset with z<0) in bins of Delta z=0.1
 * and do the shamming within each bin.
 */
float density2host_halo_zbins(float z)
{
#define NZBIN 200
  int i, iz;
  float rlo, rhi, dz, dzmin;
  static int flag = 1;
  static int zcnt[NZBIN];
  static float volume[NZBIN], zlo[NZBIN], zhi[NZBIN];

  // if first call, get the volume in each dz bin
  if(flag)
    {
      for(i=0;i<NZBIN;++i)
	{
	  zlo[i] = i*1./(NZBIN-1);
	  zhi[i] = zlo[i] + 0.05;
	  if(i==0) rlo = 0;
	  else rlo = distance_redshift(zlo[i]);
	  rhi = distance_redshift(zhi[i]);
	  volume[i] = 4./3.*PI*(rhi*rhi*rhi-rlo*rlo*rlo)*FRAC_AREA;
	  //fprintf(stderr,"volume[%d]= %e %f\n",i,volume[i],r);
	}
      flag = 0;
    }
  // if negative redshift, reset the counters;
  if(z<0)
    {
      //fprintf(stderr,"Resetting sham counts\n");
      for(i=0;i<NZBIN;++i)
	zcnt[i] =0;
    }
  if(z<0)return 0;

  // what bins does this galaxy belong to?
  dzmin = 1;
  for(i=0;i<NZBIN;++i)
    {
      if(z>=zlo[i] && z<=zhi[i]) zcnt[i]++;
      dz = fabs(z-(zhi[i]+zlo[i])/2);
      if(dz<dzmin) {
	dzmin = dz;
	iz = i;
	//fprintf(stderr," %d %f\n",i,dzmin);
      }
    }
  //fprintf(stderr,"%f %d %d %e %f %f %f\n",z,iz,zcnt[iz],zcnt[iz]/volume[iz],zlo[iz],zhi[iz],dzmin);
  return density2host_halo(zcnt[iz]/volume[iz]);

#undef NZBIN
}
/* Using a vmax correction for galaxies that can't make it
 * to the end of the redshift bin.
 */
float density2host_halo_zbins3(float z, float vmax)
{
#define NZBIN 200
  int i, iz;
  float rlo, rhi, dz, dzmin,vv;
  static int flag = 1, negcnt[NZBIN];
  static double zcnt[NZBIN];
  static float volume[NZBIN], zlo[NZBIN], zhi[NZBIN],
    vhi[NZBIN], vlo[NZBIN];

  // if first call, get the volume in each dz bin
  if(flag)
    {
      for(i=0;i<NZBIN;++i)
	{
	  zlo[i] = i*1./NZBIN;
	  zhi[i] = zlo[i] + 0.05;
	  if(i==0) rlo = 0;
	  else rlo = distance_redshift(zlo[i]);
	  rhi = distance_redshift(zhi[i]);
	  volume[i] = 4./3.*PI*(rhi*rhi*rhi-rlo*rlo*rlo)*FRAC_AREA;
	  vhi[i] = 4./3.*PI*rhi*rhi*rhi*FRAC_AREA;
	  vlo[i] = 4./3.*PI*rlo*rlo*rlo*FRAC_AREA;
	  //fprintf(stderr,"volume[%d]= %e %f\n",i,volume[i],r);
	}
      flag = 0;
    }
  // if negative redshift, reset the counters;
  if(z<0)
    {
      //fprintf(stderr,"Resetting sham counts\n");
      for(i=0;i<NZBIN;++i)
	zcnt[i] = negcnt[i] = 0;
      return 0;
    }

  if(z>100) {
    for(i=0;i<NZBIN;++i)
      if(negcnt[i])
	fprintf(stderr,"%d %f %d\n",i,zhi[i]-0.025,negcnt[i]);
    return 0;
  }
  
  // what bins does this galaxy belong to?
  dzmin = 1;
  for(i=0;i<NZBIN;++i)
    {
      if(z>=zlo[i] && z<zhi[i]) {
	if(vmax>vhi[i])vv = volume[i];
	else vv = vmax-vlo[i];
	if(vv<0) { vv = volume[i]; negcnt[i]++; }
	zcnt[i] += 1/vv;
	//fprintf(stdout,"> %d %e %e %e %e\n",i,vv,vmax,vlo[i],vhi[i]);
	//fflush(stdout);
      }
      dz = fabs(z-(zhi[i]+zlo[i])/2);
      if(dz<dzmin) {
	dzmin = dz;
	iz = i;
      }
    }
  //fprintf(stdout,"%f %d %e %e %f %f %f\n",z,iz,zcnt[iz],vmax,zlo[iz],zhi[iz],dzmin);
  //fflush(stdout);
  return density2host_halo(zcnt[iz]);

#undef NZBIN
}

/* Backup of previous version, 
 * non-overlapping redshift bins.
 */
float density2host_halo_zbins2(float z)
{
#define NZBIN 20
  int i, iz;
  float r;
  static int flag = 1;
  static int zcnt[NZBIN];
  static float volume[NZBIN];

  // if first call, get the volume in each dz bin
  if(flag)
    {
      for(i=0;i<NZBIN;++i)
	{
	  r = distance_redshift((i+1.)/NZBIN);
	  volume[i] = 4./3.*PI*r*r*r*FRAC_AREA;
	  if(i>0)volume[i] = volume[i] - volume[i-1];
	  //fprintf(stderr,"volume[%d]= %e %f\n",i,volume[i],r);
	}
      flag = 0;
    }
  // if negative redshift, reset the counters;
  if(z<0)
    {
      //fprintf(stderr,"Resetting sham counts\n");
      for(i=0;i<NZBIN;++i)
	zcnt[i] =0;
    }
  if(z<0)return 0;
  iz = (int)(z*NZBIN);
  zcnt[iz]++;
  //fprintf(stderr,"%f %d %d %e\n",z,iz,zcnt[iz],zcnt[iz]/volume[iz]);
  return density2host_halo(zcnt[iz]/volume[iz]);
#undef NZBIN
}


float func_match_ngal(float mass)
{
  static int flag=1, n=100, prev_cosmo=-1, call_count=1;
  static float *mh, *ms, *mx, *nh, mlo, mhi, dlogm, mmax;
  int i;
  float a, maglo, maghi, dmag, m, n1, n2;

  
  if(flag)
    {
      flag = 0;
      mh = vector(1,n);
      nh = vector(1,n);
      mx = vector(1,n);

      mlo = 1.0E+8;
      mhi = HALO_MAX;
      dlogm = log(mhi/mlo)/n;

      for(i=1;i<=n;++i)
	{
	  mh[i] = exp((i-0.5)*dlogm)*mlo;
	  n1 = qromo(halo_abundance2,log(mh[i]),log(HALO_MAX),midpnt);
	  if(mh[i]<HALO_MAX/1.0)
	    n2 = qromo(subhalo_abundance2,log(mh[i]),log(HALO_MAX/1.0),midpnt);
	  else
	    n2 = 0;
	  if(n2<0)n2=0;
	  nh[i] = log(n1+n2);
	  //nh[i] = log(qromo(func_nhalo,log(mh[i]),log(HALO_MAX),midpnt));
	  //printf("SHAM %e %e %e %e %e\n",(mh[i]),exp(nh[i]),n1,n2,subhalo_abundance2(log(mh[i])));
	  mh[i] = log(mh[i]);
	  fflush(stdout);
	}
      spline(mh,nh,n,1.0E+30,1.0E+30,mx);
    }
  splint(mh,nh,mx,n,mass,&a);
  return exp(a)-g7_ngal;
}

float func_match_nhost(float mass)
{
  static int flag=1, n=100, prev_cosmo=-1, call_count=1;
  static float *mh, *ms, *mx, *nh, mlo, mhi, dlogm, mmax;
  int i;
  float a, maglo, maghi, dmag, m, n1, n2;

  
  if(flag)
    {
      flag = 0;
      mh = vector(1,n);
      nh = vector(1,n);
      mx = vector(1,n);

      mlo = 1.0E+8;
      mhi = HALO_MAX;
      dlogm = log(mhi/mlo)/n;

      for(i=1;i<=n;++i)
	{
	  mh[i] = exp((i-0.5)*dlogm)*mlo;
	  n1 = qromo(halo_abundance2,log(mh[i]),log(HALO_MAX),midpnt);
	  nh[i] = log(n1);
	  mh[i] = log(mh[i]);
	  fflush(stdout);
	}
      spline(mh,nh,n,1.0E+30,1.0E+30,mx);
    }
  splint(mh,nh,mx,n,mass,&a);
  return exp(a)-g7_ngal;
}


// recall that the input mmin is already in log units
float func_match_nhalo(float mmin)
{
  return qromo(func_nhalo,mmin,log(HALO_MAX),midpnt) - g7_ngal;
}

float func_nhalo(float m)
{
  m = exp(m);
  //printf("%e %e\n",m,subhalo_abundance(m));
  return (halo_abundance(m)+subhalo_abundance(m))*m;
}

/* integrate over the parent halo mass function
 * to get the density of subhalos of mass m
 */
float subhalo_abundance(float m)
{  
  g7_msub = m;
  return qromo(func_subhalo1,log(m),log(HALO_MAX),midpnt);
}

float subhalo_abundance2(float m)
{  
  g7_msub = exp(m);
  //printf("%e\n",g7_msub);
  return qromo(func_subhalo1,(m),log(HALO_MAX),midpnt)*g7_msub;
}


/* NB the 1.5 factor is to fit the results from the z=0 200 box
 */
float func_subhalo1(float mhost)
{
  double x;
  mhost = exp(mhost);
  x= pow(g7_msub/mhost,-0.7)*exp(-9.9*pow(g7_msub/mhost,2.5))*0.3*halo_abundance(mhost)*mhost/g7_msub;
  //printf("%e %e\n",mhost,x)
  return x;
  //  return pow(g7_msub/mhost,-0.8)*exp(-g7_msub/mhost*1.25)*0.2*halo_abundance(mhost)*mhost/g7_msub;
}
 
float halo_abundance2(float m)
{
  m=exp(m);
  return halo_abundance(m)*m;
}

float halo_abundance(float m)
{
  int i;
  FILE *fp;
  float a;
  static int n=0;
  static float *x, *y, *z;
  char aa[1000];
  
  if(!n)
    {
      fp = openfile("halo_mass_function.dat");
      //fp = openfile("wmap1.massfunc");
      //fp = openfile("s8_0.7.massfunc");
      n = filesize(fp);
      x = vector(1,n);
      y = vector(1,n);
      z = vector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f",&x[i],&y[i]);
	  x[i] = log(x[i]);
	  y[i] = log(y[i]);
	  fgets(aa,1000,fp);
	}
      spline(x,y,n,1.0E+30,1.0E+30,z);
      fclose(fp);
    }
  splint(x,y,z,n,log(m),&a);
  return exp(a);
}

float func_ngal(float mag)
{
  float phi1=0.0156;
  float phi2=0.0062;
  float mstar=-20.04;
  float a1=-0.17;
  float a2=-1.52;

  // temp! replacing with blanton 03 for a second...
  return 0.4*2.30258*exp(-pow(10.0,-0.4*(mag +20.44)))*(1.49e-2*pow(10.0,-0.4*(mag + 20.44)*(-1.05+1)));

  // blanton's 05 lowL LF
  return 0.4*2.30258*exp(-pow(10.0,-0.4*(mag - mstar)))*
    (phi1*pow(10.0,-0.4*(mag - mstar)*(a1+1)) + phi2*pow(10.0,-0.4*(mag-mstar)*(a2+1)));

}
