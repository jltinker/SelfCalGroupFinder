// Initialization //

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include "nrutil.h"
#include "kdtree.h"
#include "groups.h"

// Definitions
#define OMEGA_M 0.25
#define PI 3.141592741
#define RHO_CRIT 2.775E+11
#define DELTA_HALO 200
#define SPEED_OF_LIGHT 3.0E+5
#define c_on_H0 2997.92
#define BIG_G 4.304E-9 /* BIG G in units of (km/s)^2*Mpc/M_sol */
#define G0 (1.0/sqrt(2.0*3.14159))
#define ROOT2 1.41421
#define Q0 2.0
#define Q1 -1.0
#define QZ0 0.1
#define THIRD (1.0/3.0)
#define ANG (PI/180.0)
#define RT2PI 2.50663

#define NBINS 5 // TNG has 6 bins

#define NRANDOM 1000000

/* Global for the random numbers
 */
float UNIFORM_RANDOM[NRANDOM];
float GAUSSIAN_RANDOM[NRANDOM];
int IRAN_CURRENT[100];

/* Globals for the halos
 */
struct halo *HALO;
int NHALO;
float BOX_SIZE = 250.0;
float BOX_EPSILON = 0.01;

/* Globals for the tabulated HODs
 */
double ncenr[NBINS][200], nsatr[NBINS][200],
  ncenb[NBINS][200], nsatb[NBINS][200], nhalo[NBINS][200];
int IMAG, BLUE_FLAG;

float REDSHIFT = 0.0,
  CVIR_FAC=1.0;

/* local functions
 */
float NFW_position(float mass, float x[], int thisTask);
float NFW_velocity(float mass, float v[], int thisTask);
float NFW_density(float r, float rs, float ps);
double poisson_prob(int n, float nave);
int poisson_deviate(float nave, int thisTask);
float halo_concentration(float mass);
float N_sat(float m, int imag, int blue_flag);
float N_cen(float m, int imag, int blue_flag);
void boxwrap_galaxy(float xh[], float xg[]);
float tabulated_gaussian_random(int thisTask);
float tabulated_uniform_random(int thisTask);


/* Given the group results, calculate the mean
 * expected Lsat50 values as a function of central 
 * luminosity (or mass). The code reads in a pre-tabulated
 * list matching halo mass to Lsat50 (from the C250 box).
 */
void lsat_model()
{
  FILE *fp;
  int i,j,k,n,nt, nhb[150], nhr[150], im, ihm, ix;
  float *mx, *lx, lsatb[150], lsatr[150], *m2x, x;
  int **nhrx, **nhbx;
  float **lsatbx, **lsatrx;
  int **nhrx2, **nhbx2;
  float **lsatbx2, **lsatrx2;
  int NBINS2;
  float dpropx, M0_PROPX=8.75, DM_PROPX=0.5;
  double lsat;
  int nsat;
  nhrx = imatrix(0,5,-20,20);
  nhbx = imatrix(0,5,-20,20);
  lsatbx = matrix(0,5,-20,20);
  lsatrx = matrix(0,5,-20,20);

  nhrx2 = imatrix(0,5,-20,20);
  nhbx2 = imatrix(0,5,-20,20);
  lsatbx2 = matrix(0,5,-20,20);
  lsatrx2 = matrix(0,5,-20,20);

  // this is binning for the mocks (rnage = -2,2)
  NBINS2 = 10;
  dpropx = 0.2;

  /* for the actual SDSS fitting, the
   * number of propx bins is 11, going from -2 to 2, 
   * with binwidth of 0.4
   */
  NBINS2 = 5;
  dpropx = 0.4;

  /* This is the binning for FIT3,
   * where we're doing 1 mass bin in c50/90
   */
  NBINS2 = 12;
  dpropx = 0.2;
  M0_PROPX = 0;
  DM_PROPX = 100; // evertyhing should be in one mass bin
  
  for(i=0;i<=5;++i)
    for(j=-20;j<=20;++j)
      nhbx[i][j] = nhrx[i][j] = lsatbx[i][j] = lsatrx[i][j] = 0;
  
  for(i=0;i<=5;++i)
    for(j=-20;j<=20;++j)
      nhbx2[i][j] = nhrx2[i][j] = lsatbx2[i][j] = lsatrx2[i][j] = 0;
  
  fp = openfile("lsat_lookup.dat");
  nt = filesize(fp);
  mx = vector(1,nt);
  lx = vector(1,nt);
  m2x = vector(1,nt);
  for(i=1;i<=nt;++i)
    fscanf(fp,"%f %f",&mx[i],&lx[i]);
  fclose(fp);

  spline(mx,lx,nt,1.0E+30,1.0E+30,m2x);

  for(i=0;i<150;++i)
    nhr[i] = lsatr[i] = nhb[i] = lsatb[i]  = 0;

  for(i=0;i<NGAL;++i)
    {
      if(GAL[i].psat>0.5)continue;
      // binning for lsat-vs-lum
      im = (int)(log10(GAL[i].mstellar)/0.1+0.5);
      splint(mx,lx,m2x,nt,log10(GAL[i].mass),&x);
      if(GAL[i].color>0.8)
	{
	  nhr[im]++;	  
	  lsatr[im]+=pow(10.0,x);
	}
      else
	{
	  nhb[im]++;
	  lsatb[im]+=pow(10.0,x);
	}
      // let's print this out just for kick
      //printf("LSAT %d %d %f %f %f %f\n",i,im,log10(GAL[i].mass),x,log10(GAL[i].mstellar),GAL[i].color);
      
      if(!SECOND_PARAMETER)continue;

      // if FIT3 binning, only do z<0.15
      if(GAL[i].redshift>0.15)continue;
      
      // binning for lsat-vs-propx at fixed lum
      im = (int)((log10(GAL[i].mstellar)-M0_PROPX)/DM_PROPX);
      //printf("BINX %d %f %f\n",im,log10(GAL[i].mstellar),(log10(GAL[i].mstellar)-M0_PROPX-DM_PROPX/2)/DM_PROPX);
      if(im<0 || im>=5)continue;
      ix = (int)floor((GAL[i].propx+dpropx/2)/dpropx); //for the mocks
      // does this work for negative bins?
      if(ix<-NBINS2 || ix>NBINS2)goto NEXT_PROP;
      if(GAL[i].color>0.8)
	{
	  nhrx[im][ix]++;
	  lsatrx[im][ix]+=pow(10.0,x);
	  //lsatrx[im][ix]+=GAL[i].mass;
	}
      if(GAL[i].color<0.8)
	{
	  nhbx[im][ix]++;
	  lsatbx[im][ix]+=pow(10.0,x);
	  //lsatbx[im][ix]+=GAL[i].mass;
	}
      
    NEXT_PROP:
      if(SECOND_PARAMETER==1)continue;
      ix = (int)((GAL[i].propx2-0.1)*5);
      if(ix<-10 || ix>10)continue;
      if(GAL[i].color>0.8)
	{
	  nhrx2[im][ix]++;
	  lsatrx2[im][ix]+=pow(10.0,x);
	}
      if(GAL[i].color<0.8)
	{
	  nhbx2[im][ix]++;
	  lsatbx2[im][ix]+=pow(10.0,x);
	}
	  
    }

  // output this to a pre-specified file
  // (plus we knoe the limits of the data)
  fp = fopen("lsat_groups.out","w");
  if(STELLAR_MASS)
    {
      for(i=91;i<=113;++i) // UM
	fprintf(fp,"%e %e %e\n",i/10.0,log10(lsatr[i]/nhr[i]),log10(lsatb[i]/nhb[i]));
      fclose(fp);
    }
  else
    {
      for(i=88;i<=107;++i) //C250
      //for(i=88;i<=119;++i) //TNG
	fprintf(fp,"%e %e %e\n",i/10.0,log10(lsatr[i]/nhr[i]),log10(lsatb[i]/nhb[i]));
      fclose(fp);
    }
  if(SECOND_PARAMETER==0)return;
  // check if we're doing the FIT3 single bin
  if(DM_PROPX>10) {
    fp = fopen("lsat_groups_propx_red.out","w");
    // get the mean lsat to do internal normalization
    lsat = 0;
    nsat = 0;
    for(i=-NBINS2;i<=NBINS2;++i)
      {
	lsat += lsatrx[0][i];
	nsat += nhrx[0][i];
      }
    for(i=-NBINS2;i<=NBINS2;++i)
      fprintf(fp,"%.1f %e\n",i*dpropx,lsatrx[0][i]/(nhrx[0][i]+1.0E-10)*nsat/lsat);
    fclose(fp);
    
    fp = fopen("lsat_groups_propx_blue.out","w");
    // get the mean lsat to do internal normalization
    lsat = 0;
    nsat = 0;
    for(i=-NBINS2;i<=NBINS2;++i)
      {
	lsat += lsatbx[0][i];
	nsat += nhbx[0][i];
      }
    for(i=-NBINS2;i<=NBINS2;++i)
      fprintf(fp,"%.1f %e\n",i*dpropx,lsatbx[0][i]/(nhbx[0][i]+1.0E-10)*nsat/lsat);
    fclose(fp);
    return;
  }
  // print out the correlations with propx
  // at fixed stellar mass
  fp = fopen("lsat_groups_propx_red.out","w");
  for(i=-NBINS2;i<=NBINS2;++i)
    {
      fprintf(fp,"%.1f %e %e %e %e %e\n",i*dpropx,lsatrx[0][i]/(nhrx[0][i]+1.0E-10),
	      lsatrx[1][i]/(nhrx[1][i]+1.0E-10),
	      lsatrx[2][i]/(nhrx[2][i]+1.0E-10),
	      lsatrx[3][i]/(nhrx[3][i]+1.0E-10),
	      lsatrx[4][i]/(nhrx[4][i]+1.0E-10));
    }
  fclose(fp);

  fp = fopen("lsat_groups_propx_blue.out","w");
  for(i=-NBINS2;i<=NBINS2;++i)
    {
      fprintf(fp,"%.1f %e %e %e %e %e\n",i*dpropx,lsatbx[0][i]/(nhbx[0][i]+1.0E-10),
	      lsatbx[1][i]/(nhbx[1][i]+1.0E-10),
	      lsatbx[2][i]/(nhbx[2][i]+1.0E-10),
	      lsatbx[3][i]/(nhbx[3][i]+1.0E-10),
	      lsatbx[4][i]/(nhbx[4][i]+1.0E-10));
    }
  fclose(fp);

  if(SECOND_PARAMETER==1)return;

  // print out the correlations with propx
  // at fixed stellar mass
  fp = fopen("lsat_groups_propx2_red.out","w");
  for(i=-10;i<=10;++i)
    {
      fprintf(fp,"%.1f %e %e %e %e %e\n",i/5.0,lsatrx2[0][i]/(nhrx2[0][i]+1.0E-10),
	      lsatrx2[1][i]/(nhrx2[1][i]+1.0E-10),
	      lsatrx2[2][i]/(nhrx2[2][i]+1.0E-10),
	      lsatrx2[3][i]/(nhrx2[3][i]+1.0E-10),
	      lsatrx2[4][i]/(nhrx2[4][i]+1.0E-10));
    }
  fclose(fp);

  fp = fopen("lsat_groups_propx2_blue.out","w");
  for(i=-10;i<=10;++i)
    {
      fprintf(fp,"%.1f %e %e %e %e %e\n",i/5.0,lsatbx2[0][i]/(nhbx2[0][i]+1.0E-10),
	      lsatbx2[1][i]/(nhbx2[1][i]+1.0E-10),
	      lsatbx2[2][i]/(nhbx2[2][i]+1.0E-10),
	      lsatbx2[3][i]/(nhbx2[3][i]+1.0E-10),
	      lsatbx2[4][i]/(nhbx2[4][i]+1.0E-10));
    }
  fclose(fp);

  return;
}


/* tabluate the HODs for magnitude bins (specified below)
 * for red/blue subsamples. This code assumes a flux-limited 
 * sample, so it partitions the data into volume-limited
 * samples appropriate for each mag bin (once again, the redshift
 * limits are specific below.
 *
 * The results are left in globals, which are then called to 
 * populate the simulations for measuring clustering.
 */
void tabulate_hods()
{
  FILE *fp;
  int i, j, im, igrp, ibin;
  float mag;
  // these are for MXXL-BGS
  //float magbins[NBINS] = { -17, -18, -19, -20, -21 };
  //float maxz[NBINS] = { 0.0633186, 0.098004, 0.150207, 0.227501, 0.340158 };
  // these are for MXXL-SDSS
  //float magbins[NBINS] = { -17, -18, -19, -20, -21 };
  //float maxz[NBINS] = { 0.0292367, 0.0458043, 0.0713047, 0.110097, 0.16823 };
  // these are for SHAM-SDSS (r=17.5)
  float magbins[NBINS] = { -17, -18, -19, -20, -21 };
  float maxz[NBINS] = {0.02586, 0.0406, 0.06336, 0.0981, 0.1504 };
  float volume[NBINS];
  // these are for TNG300
  //float magbins[NBINS] = { -17, -18, -19, -20, -21, -22 };
  //float maxz[NBINS] = { 0.0633186, 0.098004, 0.150207, 0.227501, 0.340158, 0.5 };
  float w0 = 1.0, w1 = w1=0.0;

  if(FLUXLIM==0)
    {
      maxz[0] = MAXREDSHIFT;
      maxz[1] = MAXREDSHIFT;
      maxz[2] = MAXREDSHIFT;
      maxz[3] = MAXREDSHIFT;
      maxz[4] = MAXREDSHIFT;
      //maxz[4] = MAXREDSHIFT;
    }
  if(STELLAR_MASS)
    {
      magbins[0] = 9.0;
      magbins[1] = 9.5;
      magbins[2] = 10.0;
      magbins[3] = 10.5;
      magbins[4] = 11.0;
    }
  
  for(i=0;i<5;++i)
    {
      volume[i] = 4./3.*PI*pow(distance_redshift(maxz[i]),3.0)*FRAC_AREA;
      for(j=0;j<200;++j)
	ncenr[i][j] = nsatr[i][j] = nhalo[i][j] = 
	  ncenb[i][j] = nsatb[i][j] = 0;
    }
  fprintf(stderr,"Tabulating HODs...\n");
  
  for(i=0;i<NGAL;++i)
    {
      // what is host halo mass?
      if(GAL[i].psat>0.5)
	{
	  igrp = GAL[i].igrp;
	  im = log10(GAL[igrp].mass)/0.1;
	}
      else
	{
	  im = log10(GAL[i].mass)/0.1;
	  for(j=0;j<NBINS;++j)
	    if(maxz[j]>GAL[i].redshift) {
	      w0 = 1/volume[j];
	      if(GAL[i].vmax<volume[j])
		w0 = 1/GAL[i].vmax;
	      nhalo[j][im]+=w0;
	    }
	}

      // check the magnitude of the galaxy
      mag = -2.5*log10(GAL[i].mstellar)+4.65;      
      ibin = (int)(fabs(mag))-17;
      //printf("BOO %f %f %d\n",log10(GAL[i].mstellar),mag,ibin);
      // check if we are using stellar mass
      if(STELLAR_MASS) {
	mag = log10(GAL[i].mstellar)*2;
	ibin = (int)(mag)-18;
      }
      
      if(ibin<0 || ibin>=NBINS)continue;
      if(GAL[i].redshift>maxz[ibin])continue;
      if(im<0 || im>=200)fprintf(stderr,"err> %d %e\n",im,GAL[i].mass);

      // vmax-weight everything
      w0 = 1/volume[ibin];
      if(GAL[i].vmax<volume[ibin])
	w0 = 1/GAL[i].vmax;

      
      if(GAL[i].color>0.8) //red
	{
	  if(GAL[i].psat>0.5){
	    nsatr[ibin][im]+=w0;
	  } else {
	    ncenr[ibin][im]+=w0;
	  }
	}
      else //blue
	{
	  if(GAL[i].psat>0.5){
	    nsatb[ibin][im]+=w0;
	  } else {
	    ncenb[ibin][im]+=w0;
	  }
	}
	  
    }
  //fprintf(stderr,"printing out\n");
  // print out the tabulated hods
  if(1)
    {
      fp = fopen("hod.out","w");
      for(i=100;i<155;++i)
	{
	  fprintf(fp,"HOD %f ",i/10.0);
	  for(j=0;j<NBINS;++j)
	    fprintf(fp,"%e %e ",ncenr[j][i]*1./(nhalo[j][i]+1.0E-20),
		   nsatr[j][i]*1./(nhalo[j][i]+1.0E-20));
	  for(j=0;j<NBINS;++j)
	    fprintf(fp,"%e %e ",ncenb[j][i]*1./(nhalo[j][i]+1.0E-20),
		   nsatb[j][i]*1./(nhalo[j][i]+1.0E-20));
	  fprintf(fp,"MF %f ",i/10.0);
	  for(j=0;j<NBINS;++j)
	    fprintf(fp,"%e ",nhalo[j][i]);
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

  for(i=90;i<200;++i)
    {
      for(j=0;j<NBINS;++j)
	{
	  ncenr[j][i] = log10(ncenr[j][i]*1./(nhalo[j][i]+1.0E-20)+1.0E-10);
	  nsatr[j][i] = log10(nsatr[j][i]*1./(nhalo[j][i]+1.0E-20)+1.0E-10);
	  ncenb[j][i] = log10(ncenb[j][i]*1./(nhalo[j][i]+1.0E-20)+1.0E-10);
	  nsatb[j][i] = log10(nsatb[j][i]*1./(nhalo[j][i]+1.0E-20)+1.0E-10);
	}
    }
}

/* Do the same as above, but now giving
 * each halo an actual value of Lsat from the simulation.
 * Thus, we test if scatter (and fraction of Lsat=0's)
 * makes any difference.
 * --------------------------------
 * RESULTS: it makes no difference.
 */
void lsat_model_scatter()
{
  FILE *fp;
  int i,j,k,id, n,nt, nhb[150], nhr[150], im, ihm;
  float *mx, *lx, *m2x, x, lsat;
  double lsatb[150], lsatr[150];
  int *indx;
  float *mvir;

  indx = ivector(1,NHALO);
  mvir = vector(1,NHALO);
  for(i=1;i<=NHALO;++i)
    {
      indx[i] = i;
      mvir[i] = HALO[i].mass;
    }
  sort2(NHALO,mvir,indx);
  fprintf(stderr,"lsat> done sort\n");
  

  for(i=0;i<150;++i)
    nhr[i] = lsatr[i] = nhb[i] = lsatb[i]  = 0;

  for(i=0;i<NGAL;++i)
    {
      if(GAL[i].psat>0.5)continue;
      im = (int)(log10(GAL[i].mstellar)/0.1+0.5);
      id = search(NHALO, mvir, GAL[i].mass);
      lsat = HALO[indx[id]].lsat;
      //printf("LSAT %d %e %e %e\n",id,GAL[i].mass, HALO[indx[id]].mass, lsat);
      
      if(GAL[i].color>0.8)
	{
	  nhr[im]++;	  
	  lsatr[im]+=lsat;
	}
      else
	{
	  nhb[im]++;
	  lsatb[im]+=lsat;
	}
    }

  // output this to a pre-specified file
  // (plus we knoe the limits of the data)
  fp = fopen("lsat_groups.out","w");
  for(i=88;i<=106;++i)
    fprintf(fp,"%e %e %e\n",i/10.0,log10(lsatr[i]/nhr[i]),log10(lsatb[i]/nhb[i]));
  fclose(fp);
  return;
}

void populate_simulation_omp(int imag, int blue_flag, int thisTask)
{
  //static int flag=1;
  int i;
  FILE *fp;
  long IDUM3=-555, iseed=555;
  FILE *outf;
  char fname[1000];
  int j,n1,imag_offset, imag_mult, istart, iend;
  float nsat, ncen, mass, xg[3],vg[3],xh[3], logm, bfit;
  struct drand48_data drand_buf;
  double r;
  

  
  if(imag<0)
    {
      srand48(555);
      fprintf(stderr,"popsim> reading halo data...\n");
      //flag = 0;
      //fp = openfile("/export/sirocco1/tinker/SIMULATIONS/BOLSHOI/hosthalo_z0.0_M1e10.dat");

      if(!(fp = fopen("/export/sirocco2/tinker/SIMULATIONS/C250_2560/hosthalo_z0.0_M1e10_Lsat.dat","r"))) {	 
	  fp = fopen("/mount/sirocco2/tinker/SIMULATIONS/C250_2560/hosthalo_z0.0_M1e10_Lsat.dat","r");
	}
      NHALO = filesize(fp);
      HALO = calloc(NHALO, sizeof(struct halo));
      for(i=0;i<NHALO;++i)
	{
	  fscanf(fp,"%f %f %f %f %f %f %f %f", &HALO[i].mass,
		 &HALO[i].x, &HALO[i].y, &HALO[i].z,
		 &HALO[i].vx, &HALO[i].vy, &HALO[i].vz,&HALO[i].lsat);
	}
      fclose(fp);
      fprintf(stderr,"popsim> done reading halo data [%d].\n",NHALO);

      // lets create a list of random numbers
      fprintf(stderr,"popsim> creating random numbers.\n",NHALO);
      for(i=0;i<NRANDOM;++i)
	{
	  UNIFORM_RANDOM[i] = drand48();
	  GAUSSIAN_RANDOM[i] = gasdev(&IDUM3);
	}
      // each task gets its own counter
      for(i=0;i<100;++i)
	IRAN_CURRENT[i] = (int)(drand48()*100);
      fprintf(stderr,"popsim> done with randoms.\n",NHALO);
      
      return;
    }

  /* Put this in a global so that we know which HOD
   * to use.
   */
  fprintf(stderr,"popsim> starting population for imag=%d, blue=%d\n",imag, blue_flag);

  imag_offset = 17;
  imag_mult = 1;
  if(STELLAR_MASS){
    imag_offset = 90;
    imag_mult = 5;
  }
  /* We'll do simple linear interpolation for the HOD
   */
  if(blue_flag)
    sprintf(fname,"mock_blue_M%d.dat",imag*imag_mult+imag_offset);
  else
    sprintf(fname,"mock_red_M%d.dat",imag*imag_mult+imag_offset);
  outf = fopen(fname,"w");

  //srand48_r (iseed, &drand_buf);

  // fit the high-mass satellite occupation function, force slope=1
  istart = 130;
  if(imag>=2)istart = 135;  
  iend = 140;
  if(imag>=2)iend = 145;
  bfit = 0;
  if(imag==0) {
    istart = 120;
    iend = 130;
  }
  if(blue_flag)
    {
      for(i=istart;i<=iend;++i)
	bfit += nsatb[imag][i] - i/10.0; 
      bfit = bfit/(iend-istart+1);
      for(i=iend;i<=160;++i)
	nsatb[imag][i] = 1*i/10.0 + bfit;      
    }
  else
    {
      for(i=istart;i<=iend;++i)
	bfit += nsatr[imag][i] - i/10.0; 
      bfit = bfit/(iend-istart+1);
      for(i=iend;i<=160;++i)
	nsatr[imag][i] = 1*i/10.0 + bfit;      
    }    
  
  for(i=0;i<NHALO;++i)
    {
      mass = HALO[i].mass;
      logm = log10(mass);
      ncen=N_cen(mass,imag,blue_flag);
      //drand48_r(&drand_buf, &r);
      r = tabulated_uniform_random(thisTask);
      //r = UNIFORM_RANDOM[IRAN_CURRENT[thisTask]];
      //IRAN_CURRENT[thisTask]++;
      //if(IRAN_CURRENT[thisTask]==NRANDOM)IRAN_CURRENT[thisTask]=0;
      if(r<ncen)
	{
	  fprintf(outf,"%.5f %.5f %.5f %f %f %f %d %f\n",
		  HALO[i].x, HALO[i].y, HALO[i].z,
		  HALO[i].vx, HALO[i].vy, HALO[i].vz,0,logm);
	}
      nsat = N_sat(mass,imag,blue_flag);      
      n1 = poisson_deviate(nsat,thisTask);
      for(j=1;j<=n1;++j)
	{
	  NFW_position(mass,xg,thisTask);
	  NFW_velocity(mass,vg,thisTask);
	  xh[0] = HALO[i].x; xh[1] = HALO[i].y; xh[2] = HALO[i].z;
	  boxwrap_galaxy(xh,xg);
	  if(isnan(xg[0]+xg[1]+xg[2]))continue;
	  fprintf(outf,"%.5f %.5f %.5f %f %f %f %d %f\n",
		  xg[0],xg[1],xg[2],
		  HALO[i].vx+vg[0],HALO[i].vy+vg[1],HALO[i].vz+vg[2],
		  1,logm);
	}	
    }
  fclose(outf);
  

}

void boxwrap_galaxy(float xh[], float xg[])
{
  int i;
  for(i=0;i<3;++i)
    {
      xg[i] += xh[i];
      if(xg[i]>BOX_SIZE)xg[i]-=BOX_SIZE;
      if(xg[i]<0)xg[i]+=BOX_SIZE;
      if(xg[i]>BOX_SIZE-BOX_EPSILON)xg[i] = BOX_SIZE-BOX_EPSILON;
    }
}

float N_cen(float m,int imag,int blue_flag)
{
  int im;
  float x0, y0, x1, y1, yp, logm;
  logm = log10(m)/0.1;
  im = (int)logm;
  x0 = im;
  x1 = im+1;
  if(blue_flag) {
    y0 = ncenb[imag][im];
    y1 = ncenb[imag][im+1];
  } else {
    y0 = ncenr[imag][im];
    y1 = ncenr[imag][im+1];
  }    
  yp = y0 + ((y1-y0)/(x1-x0)) * (logm - x0);
  //if(logm>124)
  // printf("CEN %e %f %d %f %f %f %f %f\n",m,logm,im,x0,x1,y0,y1,yp);
  if(yp<=-10)return 0;
  return pow(10.0,yp);
}
float N_sat(float m, int imag, int blue_flag)
{
  int im;
  float x0, y0, x1, y1, yp, logm;
  logm = log10(m)/0.1;
  im = (int)logm;
  x0 = im;
  x1 = im+1;
  if(blue_flag) {
    y0 = nsatb[imag][im];
    y1 = nsatb[imag][im+1];
  } else {
    y0 = nsatr[imag][im];
    y1 = nsatr[imag][im+1];
  }    
  yp = y0 + ((y1-y0)/(x1-x0)) * (logm - x0);
  if(yp<=-10)return 0;
  return pow(10.0,yp);
}

int poisson_deviate_old(float nave, int thisTask)
{
  struct drand48_data drand_buf;
  double p,pp;
  double r;
  int n;

  if(nave>50)
    {
      r = tabulated_uniform_random(thisTask);
      return (int)(r*sqrt(nave) + nave);
    }
  
  p=0;
  pp=1;

  while(p<pp)
    {
      r = tabulated_uniform_random(thisTask);
      if(nave<1)	
	n=(int)(r*20);
      else
	n=(int)(r*30*nave);
      p=poisson_prob(n,nave);
      pp = tabulated_uniform_random(thisTask);
      //drand48_r(&drand_buf, &pp);
    }
  return(n);
}

//===========================================================================
//=  Function to generate Poisson distributed random variables              =
//=    - Input:  Mean value of distribution                                 =
//=    - Output: Returns with Poisson distributed random variable           =
//===========================================================================
int poisson_deviate(float x, int thisTask)
{
  float r;
  int    poi_value;             // Computed Poisson value to be returned
  double t_sum;                 // Time sum value

  if(x>50)
    return (int)x;

  x = 1/x; //???
  // Loop to generate Poisson values using exponential distribution
  poi_value = 0;
  t_sum = 0.0;
  while(1)
  {
    r = tabulated_uniform_random(thisTask);
    t_sum = t_sum - x*log(r);
    //printf("POI %d %e %e %e %e\n",poi_value,x,r,log(r),t_sum);
    if (t_sum >= 1.0) break;
    poi_value++;
  }

  return(poi_value);
}


/* Poisson probability of n given n_average
 */
double poisson_prob(int n, float nave)
{
  int i;
  double fac=1;

  if(n>0)
    for(i=1;i<=n;++i)
      fac*=nave/i;

  return((double)(fac*exp(-nave)));
}

float tabulated_gaussian_random(int thisTask)
{
  float r;
  r = GAUSSIAN_RANDOM[IRAN_CURRENT[thisTask]];
  IRAN_CURRENT[thisTask]++;
  if(IRAN_CURRENT[thisTask]==NRANDOM)IRAN_CURRENT[thisTask]=0;
  return r;
}
float tabulated_uniform_random(int thisTask)
{
  float r;
  r = UNIFORM_RANDOM[IRAN_CURRENT[thisTask]];
  IRAN_CURRENT[thisTask]++;
  if(IRAN_CURRENT[thisTask]==NRANDOM)IRAN_CURRENT[thisTask]=0;
  return r;
}

/* Randomy generates a position away from the origin with 
 * a probability given by the NFW profile for a halo of the input
 * mass (and including the CVIR_FAC)
 */
float NFW_position(float mass, float x[], int thisTask)
{
  float r,pr,max_p,costheta,sintheta,phi1,signs,rvir,rs,cvir, mfac=1;
  double rr;
  cvir=halo_concentration(mass)*CVIR_FAC;
  rvir=pow(3*mass/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rs=rvir/cvir;
  max_p=NFW_density(rs,rs,1.0)*rs*rs*4.0*PI;

  for(;;) {
    r=tabulated_uniform_random(thisTask)*rvir;
    pr=NFW_density(r,rs,1.0)*r*r*4.0*PI/max_p;
    
    if(tabulated_uniform_random(thisTask)<=pr)
      {
	costheta=2.*(tabulated_uniform_random(thisTask)-.5);
	sintheta=sqrt(1.-costheta*costheta);
	signs=2.*(tabulated_uniform_random(thisTask)-.5);
	costheta=signs*costheta/fabs(signs);
	phi1=2.0*PI*tabulated_uniform_random(thisTask);
	
	x[0]=r*sintheta*cos(phi1);
	x[1]=r*sintheta*sin(phi1);
	x[2]=r*costheta;
	return r;
      }
  }
}

/* This is the NFW density profile
 */
float NFW_density(float r, float rs, float ps)
{
  return(ps*rs/(r*(1+r/rs)*(1+r/rs)));
}

/* This sets the velocity to be isotropic Gaussian.
 */
float NFW_velocity(float mass, float v[], int thisTask)
{
  static long IDUM2=-455;
  //static float fac = -1;
  float sigv,vbias=1,mfac=1;
  int i;
  float fac;

  fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)
    *3.09E19*sqrt(1+REDSHIFT);
  sigv=fac*pow(mass,1.0/3.0)/ROOT2;
  for(i=0;i<3;++i)
    v[i]=tabulated_gaussian_random(thisTask)*sigv;
  return(0);
}

float halo_concentration(float mass)
{
  return 10*pow(mass/1.0E14,-0.11);
}


