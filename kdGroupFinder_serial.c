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

/* Imported functions from numerical recipes 
 */
float qromo(float (*func)(float), float a, float b,
       float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);
float qtrap(float (*func)(float), float a, float b);

/* external functions
 */
float density2halo(float galaxy_density);
float density2host_halo(float galaxy_density);

/* Local functions
 */
float distance_redshift(float z);
float angular_separation(float a1, float d1, float a2, float d2);
void find_satellites(int icen, void *kd);
float radial_probability(float mass, float dr, float rad, float ang_rad);
float fluxlim_correction(float z);

/* Structure definition for galaxies.
 */
struct galaxy {
  float x,y,z;
  float ra, dec, redshift;
  float rco;
  float luminosity,
    magnitude,
    mstellar,
    psat,
    vmax;
  int igrp;
  int id;

  // halo properties  
  float mass,
    theta,
    rad,
    sigmav,
    mtot,
    nsat;
} *GAL;

// other global variables
int NGAL;
float BPROB=4;
float MINREDSHIFT;
float MAXREDSHIFT;

int main(int argc, char **argv)
{
  FILE *fp;
  char aa[1000];
  int i, ngrp, i1, *itmp, niter, MAX_ITER=5, j, *permanent_id, ngrp_prev;
  float frac_area, zmin, zmax, volume, *xtmp, nsat_tot;
  int fluxlim = 0, *flag;
  double galden;
  
  void *kd;


  
  if(argc<5)
    {
      fprintf(stderr,"kdGroupFinder inputfile zmin zmax frac_area [fluxlim] > out\n");
      exit(0);
    }
  fp = openfile(argv[1]);
  NGAL = filesize(fp);
  fprintf(stderr,"Allocating space for [%d] galaxies\n",NGAL);
  GAL = calloc(NGAL, sizeof(struct galaxy));
  flag = ivector(0,NGAL-1);
  
  MINREDSHIFT = zmin = atof(argv[2]);
  MAXREDSHIFT = zmax = atof(argv[3]);
  frac_area = atof(argv[4]);
  if(argc>5)
    fluxlim = atoi(argv[5]);
  fprintf(stderr,"input> %f %f %f %d\n",zmin, zmax, frac_area,fluxlim);

  // calculate the volume of the sample
  volume = 4./3.*PI*(pow(distance_redshift(zmax),3.0))*frac_area;
  volume = volume -  4./3.*PI*(pow(distance_redshift(zmin),3.0))*frac_area;

  for(i=0;i<=NGAL;++i) 
    {
      fscanf(fp,"%f %f %f %f",&GAL[i].ra,&GAL[i].dec,&GAL[i].redshift,&GAL[i].mstellar);
      GAL[i].ra *= PI/180.;
      GAL[i].dec *= PI/180.;
      GAL[i].id = i;
      GAL[i].rco = distance_redshift(GAL[i].redshift);
      // check if the stellar mass is in log
      if(GAL[i].mstellar<100)
	GAL[i].mstellar = pow(10.0,GAL[i].mstellar);      
      // check to see if we're doing a fluxlimited sample
      if(fluxlim)
	fscanf(fp,"%f",&GAL[i].vmax);
      else
	GAL[i].vmax = volume;
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in from [%s]\n",argv[1]);

  fprintf(stderr,"Volume= %e L_box= %f\n",volume, pow(volume, THIRD));
  fprintf(stderr,"Number density= %e\n",NGAL/volume);

  // first sort by stellar mass
  xtmp = vector(1,NGAL);
  itmp = ivector(1,NGAL);
  permanent_id = ivector(1,NGAL);
  for(i=1;i<=NGAL;++i)
    {
      xtmp[i] = -GAL[i-1].mstellar;
      itmp[i] = i-1;
    }
  fprintf(stderr,"sorting galaxies...\n");
  sort2(NGAL, xtmp, itmp);
  fprintf(stderr,"done sorting galaxies.\n");

  // do the inverse-abundance matching
  density2host_halo(0.01);
  fprintf(stderr,"Starting inverse-sham...\n");
  galden = 0;
  for(i1=1;i1<=NGAL;++i1)
    {      
      i= itmp[i1];
      galden += 1/GAL[i].vmax;
      GAL[i].mass = density2host_halo(galden);
      GAL[i].rad = pow(3*GAL[i].mass/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),THIRD);
      GAL[i].theta = GAL[i].rad/GAL[i].rco;
      GAL[i].sigmav = sqrt(BIG_G*GAL[i].mass/2.0/GAL[i].rad*(1+GAL[i].redshift));
      GAL[i].psat = 0;
    }
  fprintf(stderr,"Done inverse-sham.\n");

  // assume that NGAL=NGROUP at first
  ngrp = NGAL;

  // let's create a 3D KD tree
  fprintf(stderr,"Building KD-tree...\n");
  kd = kd_create(3);
  for(i = 1; i <= NGAL; ++i){
    j = itmp[i];
    GAL[j].x = GAL[j].rco * cos(GAL[j].ra) * cos(GAL[j].dec);
    GAL[j].y = GAL[j].rco * sin(GAL[j].ra) * cos(GAL[j].dec); 
    GAL[j].z = GAL[j].rco * sin(GAL[j].dec);
    permanent_id[j] = j;
    
    double pt[3] = {GAL[j].x, GAL[j].y, GAL[j].z};
    assert( kd_insert(kd, pt, (void*)&permanent_id[j]) == 0);
  }
  fprintf(stderr,"Done building KD-tree.\n");


  
  // now start the group-finding iteratin
  for(niter=1;niter<=MAX_ITER;++niter)
    {
      // first, reset the psat values
      for(j=0;j<NGAL;++j)
	{
	  GAL[j].psat = 0;
	  GAL[j].nsat = 0;
	  GAL[j].mtot = GAL[j].mstellar;
	  flag[j] = 1;
	}
      // find the satellites for each halo, in order of group mass
      ngrp_prev = ngrp;
      ngrp = 0;
      for(i1=1;i1<=ngrp_prev;++i1)
	{
	  i = itmp[i1];
	  flag[i] = 0;
	  find_satellites(itmp[i1],kd);
	  if(GAL[i].psat<0.5)
	    {
	      ngrp++;
	      xtmp[ngrp] = -GAL[i].mtot;
	      itmp[ngrp] = i;
	      if(fluxlim)
		xtmp[ngrp] *= fluxlim_correction(GAL[i].redshift);
	    }
	}
      // go back and check objects are newly-exposed centrals
      for(j=0;j<NGAL;++j)
	{
	  if(flag[j] && GAL[j].psat<0.5)
	    {
	      find_satellites(j,kd);
	      ngrp++;
	      xtmp[ngrp] = -GAL[i].mtot;
	      itmp[ngrp] = i;
	      if(fluxlim)
		xtmp[ngrp] *= fluxlim_correction(GAL[i].redshift);
	    }
	}
      
      // sort groups by their total stellar mass
      sort2(ngrp,xtmp,itmp);

      // reassign the halo masses
      nsat_tot = galden = 0;
      for(j=1;j<=ngrp;++j)
	{      
	  i= itmp[j];
	  galden += 1/GAL[i].vmax;
	  GAL[i].mass = density2host_halo(galden);
	  GAL[i].rad = pow(3*GAL[i].mass/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),THIRD);
	  GAL[i].theta = GAL[i].rad/GAL[i].rco;
	  GAL[i].sigmav = sqrt(BIG_G*GAL[i].mass/2.0/GAL[i].rad*(1+GAL[i].redshift));
	  nsat_tot += GAL[i].nsat;
	}
      fprintf(stderr,"iter %d ngroups=%d fsat=%f\n",niter,ngrp,nsat_tot/NGAL);
    }

  /* Output to disk the final results
   */
  for(i=0;i<NGAL;++i)
    {
      printf("%d %f %f %f %e %e %f %e %e %e %d\n",
	     i, GAL[i].ra*180/PI, GAL[i].dec*180/PI,GAL[i].redshift,
	     GAL[i].mstellar, GAL[i].vmax, GAL[i].psat, GAL[i].mass,
	     GAL[i].nsat, GAL[i].mtot, GAL[i].igrp);
    }
}

/* Distance-redshift relation
 */
float func_dr1(float z)
{
  return pow(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M),-0.5);
}
float distance_redshift(float z)
{
  float x;
  if(z<=0)return 0;
  x= c_on_H0*qromo(func_dr1,0.0,z,midpnt);
  return x;
}



/* Here is the main code to find satellites for a given central galaxy
 */
void find_satellites(int icen, void *kd) {
  int j, k;
  float dx, dy, dz, theta, prob_ang, vol_corr, prob_rad, grp_lum, p0, range;
  float cenDist;
  void *set;
  int *pch;
  double cen[3];
  double sat[3];

  // check if this galaxy has already been given to a group
  if(GAL[icen].psat>0.5)return;
  
  // Use the k-d tree kd to identify the nearest galaxies to the central.
  cen[0] = GAL[icen].x;
  cen[1] = GAL[icen].y;
  cen[2] = GAL[icen].z;
  
  // Nearest neighbour search should go out to about 4*sigma, the velocity dispersion of the SHAMed halo.
  // find all galaxies in 3D that are within 4sigma of the velocity dispersion
  range = 4*GAL[icen].sigmav/100.0*(1+GAL[icen].redshift)/
    sqrt(OMEGA_M*pow(1+GAL[icen].redshift,3.0) + 1-OMEGA_M);
  set = kd_nearest_range(kd, cen, range);
  
  // Set now contains the nearest neighbours within a distance range. Grab their info. 
  // Note that set will ALWAYS contain a node that is the same as the central galaxy (this is a quirk of the code--or me not using it properly--when computing nearest neighbour distances to a point that is already in the k-d tree). Make sure to reject this galaxy.
  
  while( !kd_res_end(set)) {
    
    // Get index value of the current neighbor      
    pch = (int*)kd_res_item(set, sat);
    j = *pch;
    kd_res_next(set);
    //printf("%d %d %f %f %f %f\n",j,icen,GAL[icen].x, GAL[j].x, range,sat[0]);
    
    // Skip if target galaxy is the same as the central (obviously).
    if(j == icen)continue;

    // skip if the object is more massive than the icen
    if(GAL[j].mstellar>=GAL[icen].mstellar)continue;
    
    // Skip if already assigned to a central.
    if(GAL[j].psat)continue;
    
    // check if the galaxy is outside the angular radius of the halo
    dz = fabs(GAL[icen].redshift - GAL[j].redshift)*SPEED_OF_LIGHT;
    theta = angular_separation(GAL[icen].ra,GAL[icen].dec,GAL[j].ra,GAL[j].dec);
    if(theta > GAL[icen].theta){
      continue;
    }
    
    // Now determine the probability of being a satellite
    //(both projected onto the sky, and along the line of sight).      
    prob_ang = radial_probability(GAL[icen].mass, theta, GAL[icen].rad, GAL[icen].theta);
    prob_rad = exp(-dz*dz/(2*GAL[icen].sigmav*GAL[icen].sigmav))
      *SPEED_OF_LIGHT/(RT2PI*GAL[icen].sigmav);
    
    // combine them into the total probability
    p0 = (1 - 1/(1 + prob_ang * prob_rad / BPROB));
    
    if(p0 < 0){
      printf("ZERO %e\n",p0);
      p0 = 0;
    }
    if(p0<0.5)continue;
    
    // this is considered a member of the group
    GAL[j].psat = p0;
    GAL[j].igrp = icen;
    GAL[icen].mtot += GAL[j].mstellar;
    GAL[icen].nsat++;
  }
  //exit(0);
  // Correct for boundary conditions
  dz = SPEED_OF_LIGHT* fabs(GAL[icen].redshift - MINREDSHIFT);
  vol_corr = 1-(0.5*erfc(dz/(ROOT2*GAL[icen].sigmav)));
  GAL[icen].nsat /= vol_corr;
  GAL[icen].mtot /= vol_corr;
  
  dz = SPEED_OF_LIGHT* fabs(GAL[icen].redshift - MAXREDSHIFT);
  vol_corr = 1-(0.5*erfc(dz/(ROOT2*GAL[j].sigmav)));
  GAL[icen].nsat /= vol_corr;
  GAL[icen].mtot /= vol_corr;
}

/* angular separation between two points in ra/dec
 */
float angular_separation(float a1, float d1, float a2, float d2)
{
  float cd1,cd2,sd1,sd2,ca1a2,sa1a2;

  return atan((sqrt(cos(d2)*cos(d2)*sin(a2-a1)*sin(a2-a1) + 
		    pow(cos(d1)*sin(d2) - sin(d1)*cos(d2)*cos(a2-a1),2.0)))/
	      (sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a2-a1)));
}



/* Probability assuming a projected NFW profile
 */
float radial_probability(float mass, float dr, float rad, float ang_rad)
{
  float c, x, rs, delta, f;

  dr = dr*rad/ang_rad;

  c = 10.0*pow(mass/1.0E+14,-0.11);
  rs = rad/c;
  x = dr/rs;

  if(x<1)
    f = 1/(x*x-1)*(1-log((1+sqrt(1-x*x))/x)/(sqrt(1-x*x)));
  if(x==1)
    f = 1.0/3.0;
  if(x>1)
    f = 1/(x*x-1)*(1-atan(sqrt(x*x-1))/sqrt(x*x-1));

  delta = DELTA_HALO/3.0*c*c*c/(log(1+c)-c/(1+c));

  return 1.0/c_on_H0*2*rs*delta*f;
}


/* This is calibrated from the MXXL BGS mock,
 * from ratio of luminosity density in redshift
 * bins relative to total 1/Vmax-weighted luminosity
 * density. (SLightly different than Yang et al).
 *
 * luminosity_correction.py
 */
float fluxlim_correction(float z)
{
  return 1;
  return pow(10.0,pow(z/0.4,4.0)*0.1);
}
