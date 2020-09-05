// Initialization //

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include <omp.h>
#include "nrutil.h"
#include "kdtree.h"
#include "groups.h"


struct galaxy *GAL;
int NGAL;
int OUTPUT=0;

/* Local functions
 */
float angular_separation(float a1, float d1, float a2, float d2);
void find_satellites(int icen, void *kd);
float radial_probability(float mass, float dr, float rad, float ang_rad);
float fluxlim_correction(float z);
void groupfind(void);

/* Variables for determining
 * if a galaxy is a satellite
 */
float BPROB=10;
float BPROB_RED = 5, BPROB_XRED=0;
float BPROB_BLUE = 15, BPROB_XBLUE=0;

/* Variables for weighting the
 * central galaxies of the blue galaxies
 */
float WCEN_MASS = 10.5,
  WCEN_SIG = 0.5,
  WCEN_MASSR = 10.5,
  WCEN_SIGR = 1.0,
  WCEN_NORM = 1;

float MINREDSHIFT;
float MAXREDSHIFT;
float FRAC_AREA;
int FLUXLIM, COLOR;
int ARGC;
char  **ARGV;

int main(int argc, char **argv)
{
  int i,istep,istart;
  double t0,t1;


  if(argc<5)
    {
      fprintf(stderr,"kdGroupFinder inputfile zmin zmax frac_area [fluxlim] [color] [wcenvalues 1-3] [pbvalues 1-4]> out\n");
      exit(0);
    }
  ARGC = argc;
  ARGV = argv;
  MINREDSHIFT = atof(argv[2]);
  MAXREDSHIFT = atof(argv[3]);
  FRAC_AREA = atof(argv[4]);
  if(argc>5)
    FLUXLIM = atoi(argv[5]);
  if(argc>6)
    COLOR = atoi(argv[6]);
  if(argc>7)
    {      
      WCEN_MASS = atof(argv[7]);
      WCEN_SIG = atof(argv[8]);
      WCEN_MASSR = atof(argv[9]);
      WCEN_SIGR = atof(argv[10]);
    }
  if(argc>11)
    {
      BPROB_RED = atof(argv[11]);
      BPROB_XRED = atof(argv[12]);
      BPROB_BLUE = atof(argv[13]);
      BPROB_XBLUE = atof(argv[14]);
    }
  fprintf(stderr,"input> %f %f %f %d\n",MINREDSHIFT, MAXREDSHIFT, FRAC_AREA, FLUXLIM);
  fprintf(stderr,"input> %f %f %f\n",WCEN_MASS, WCEN_SIG, WCEN_NORM);
  fprintf(stderr,"input> %f %f %f %f\n",BPROB_RED,BPROB_XRED,BPROB_BLUE,BPROB_XBLUE);

  OUTPUT = 1;
  groupfind();
  OUTPUT = 0;
  tabulate_hods();
  lsat_model();
  populate_simulation(-1,0);
  t0 = omp_get_wtime();
#pragma omp parallel private(i,istart,istep)
  {
    istart = omp_get_thread_num();
    istep = omp_get_num_threads();
    for(i=istart;i<10;i+=istep)
      {
	populate_simulation(i/2,i%2);
      }
  }
  t1 = omp_get_wtime();
  fprintf(stderr,"popsim> %.2f sec\n",t1-t0);
}

void groupfind()
{
  FILE *fp;
  char aa[1000];
  int i, i1, niter, MAX_ITER=5, j, ngrp_prev;
  float frac_area, zmin, zmax, nsat_tot;
  int fluxlim = 0, colors = 1;
  double galden, pt[3];
  long IDUM1 = -555;
  double t0,t1,t3,t4;

  static int *permanent_id, *itmp, *flag;
  static float volume, *xtmp, *lumshift;
  static void *kd;
  static int first_call=1, ngrp;

  if(first_call) {
    colors = COLOR;
    first_call = 0;
    fp = openfile(ARGV[1]);
    NGAL = filesize(fp);
    fprintf(stderr,"Allocating space for [%d] galaxies\n",NGAL);
    GAL = calloc(NGAL, sizeof(struct galaxy));
    flag = ivector(0,NGAL-1);
    
    fluxlim = FLUXLIM;
    zmin = MINREDSHIFT;
    zmax = MAXREDSHIFT;
    
    // calculate the volume of the sample
    volume = 4./3.*PI*(pow(distance_redshift(zmax),3.0))*FRAC_AREA;
    volume = volume -  4./3.*PI*(pow(distance_redshift(zmin),3.0))*FRAC_AREA;
    
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
	// check to see if we're using colors
	if(colors)
	  fscanf(fp,"%f",&GAL[i].color);

	fgets(aa,1000,fp);
      }
    fclose(fp);
    fprintf(stderr,"Done reading in from [%s]\n",ARGV[1]);
    
    fprintf(stderr,"Volume= %e L_box= %f\n",volume, pow(volume, THIRD));
    fprintf(stderr,"Number density= %e\n",NGAL/volume);
    
    // first sort by stellar mass
    xtmp = vector(1,NGAL);
    itmp = ivector(1,NGAL);
    permanent_id = ivector(1,NGAL);
    lumshift = vector(0,NGAL-1);
    for(i=1;i<=NGAL;++i)
      {
	// just for kicks, give each galaxy a random luminosity
	lumshift[i-1] = pow(10.0,gasdev(&IDUM1)*0.0);
	
	xtmp[i] = -(GAL[i-1].mstellar*lumshift[i-1]);
	itmp[i] = i-1;
      }
    fprintf(stderr,"sorting galaxies...\n");
    sort2(NGAL, xtmp, itmp);
    fprintf(stderr,"done sorting galaxies.\n");
    
    // do the inverse-abundance matching
    density2host_halo(0.01);
    fprintf(stderr,"Starting inverse-sham...\n");
    galden = 0;
    // reset the sham counters
    if(fluxlim)
      density2host_halo_zbins(-1);
    for(i1=1;i1<=NGAL;++i1)
      {      
	i= itmp[i1];
	galden += 1/GAL[i].vmax;
	if(fluxlim)
	  GAL[i].mass = density2host_halo_zbins(GAL[i].redshift);
	else
	  GAL[i].mass = density2host_halo_zbins(galden);
	GAL[i].rad = pow(3*GAL[i].mass/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),THIRD);
	GAL[i].theta = GAL[i].rad/GAL[i].rco;
	GAL[i].sigmav = sqrt(BIG_G*GAL[i].mass/2.0/GAL[i].rad*(1+GAL[i].redshift));
	GAL[i].psat = 0;
	j = i;
	GAL[j].x = GAL[j].rco * cos(GAL[j].ra) * cos(GAL[j].dec);
	GAL[j].y = GAL[j].rco * sin(GAL[j].ra) * cos(GAL[j].dec); 
	GAL[j].z = GAL[j].rco * sin(GAL[j].dec);
      }
    fprintf(stderr,"Done inverse-sham.\n");
    // assume that NGAL=NGROUP at first
    ngrp = NGAL;
  }
  
  // let's create a 3D KD tree
  fprintf(stderr,"Building KD-tree...\n");
  kd = kd_create(3);
  for(i = 1; i <= NGAL; ++i){
    j = i;
    permanent_id[j] = j;
    pt[0] = GAL[j].x;
    pt[1] = GAL[j].y;
    pt[2] = GAL[j].z;
    assert( kd_insert(kd, pt, (void*)&permanent_id[j]) == 0);
  }
  fprintf(stderr,"Done building KD-tree. %d\n",ngrp);



  // now start the group-finding iteratin
  for(niter=1;niter<=MAX_ITER;++niter)
    {
      t3 = omp_get_wtime();

      // first, reset the psat values
      for(j=0;j<NGAL;++j)
	{
	  GAL[j].psat = 0;
	  GAL[j].nsat = 0;
	  GAL[j].mtot = GAL[j].mstellar*lumshift[j];
	  if(GAL[j].color<0.8)
	    GAL[j].mtot *=
	      1.0/pow(10.0,0.5*(1+erf((log10(GAL[j].mstellar)-WCEN_MASS)/WCEN_SIG)));
	  else
	    GAL[j].mtot *=
		pow(10.0,0.5*(1+erf((log10(GAL[j].mstellar)-WCEN_MASSR)/WCEN_SIGR)));
	  flag[j] = 1;
	}
      // find the satellites for each halo, in order of group mass
      ngrp_prev = ngrp;
      ngrp = 0;
      t0 = omp_get_wtime();
      for(i1=1;i1<=ngrp_prev;++i1)
	{
	  i = itmp[i1];
	  flag[i] = 0;
	  find_satellites(itmp[i1],kd);
	}
      find_satellites(-1,kd);
      t1 = omp_get_wtime();
      for(i1=1;i1<=ngrp_prev;++i1)
	{      
	  i = itmp[i1];
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
	      xtmp[ngrp] = -GAL[j].mtot;
	      itmp[ngrp] = j;
	      if(fluxlim)
		xtmp[ngrp] *= fluxlim_correction(GAL[j].redshift);
	    }
	}
      
      // sort groups by their total stellar mass
      sort2(ngrp,xtmp,itmp);

      // reassign the halo masses
      nsat_tot = galden = 0;
      // reset the sham counters
      if(fluxlim)
	density2host_halo_zbins(-1);
      for(j=1;j<=ngrp;++j)
	{      
	  i= itmp[j];
	  galden += 1/GAL[i].vmax;
	  if(fluxlim)
	    GAL[i].mass = density2host_halo_zbins(GAL[i].redshift);
	  else
	    GAL[i].mass = density2host_halo(galden);
	  GAL[i].rad = pow(3*GAL[i].mass/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),THIRD);
	  GAL[i].theta = GAL[i].rad/GAL[i].rco;
	  GAL[i].sigmav = sqrt(BIG_G*GAL[i].mass/2.0/GAL[i].rad*(1+GAL[i].redshift));
	  nsat_tot += GAL[i].nsat;
	}
      t4 = omp_get_wtime();
      //for the satellites, set their host halo mass
      for(j=0;j<NGAL;++j)
	if(GAL[j].psat>0.5)
	  GAL[j].mass = GAL[GAL[j].igrp].mass;
      fprintf(stderr,"iter %d ngroups=%d fsat=%f (kdtime=%.2f %.2f)\n",niter,ngrp,nsat_tot/NGAL,
	      t1-t0,t4-t3);
    }

  /* Output to disk the final results
   */
  if(OUTPUT)
    {
      for(i=0;i<NGAL;++i)
	{
	  printf("%d %f %f %f %e %e %f %e %e %e %d\n",
		 i, GAL[i].ra*180/PI, GAL[i].dec*180/PI,GAL[i].redshift,
		 GAL[i].mstellar, GAL[i].vmax, GAL[i].psat, GAL[i].mass,
		 GAL[i].nsat, GAL[i].mtot, GAL[i].igrp);
	}
      fflush(stdout);
    }
  /* let's free up the memory of the kdtree
   */
  kd_free(kd);
  
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
void find_satellites(int ii, void *kd) {
  int i,nThreads;
  static int icen_set[24];
  static int cnt=0;

  if(GAL[ii].psat>0.5)return;
  
#pragma omp parallel shared(nThreads)
  {
#pragma omp  single 
  {
    nThreads = omp_get_num_threads();  }
  }
  if(ii>=0)
    {
      icen_set[cnt] = ii;
      cnt++;
      if(cnt<nThreads)return;
    }
  
  //printf("resetting count: %d %d\n",cnt,nThreads);
  //for(i=0;i<cnt;++i)
  //  printf("cnt[%d]= %d\n",i,icen_set[i]);
  cnt = 0;

#pragma omp parallel
  {
    int j, k, thisTask, icen;
  float dx, dy, dz, theta, prob_ang, vol_corr, prob_rad, grp_lum, p0, range;
  float cenDist, bprob;
  void *set;
  int *pch;
  double cen[3];
  double sat[3];
  int local_cnt = 0;

  

  
  thisTask = omp_get_thread_num();
  icen = icen_set[thisTask];
  //printf("ICEN %d %d\n",thisTask,icen);
  //fflush(stdout);
  // check if this galaxy has already been given to a group
  if(GAL[icen].psat>0.5)goto END_SAT;

  cen[0] = GAL[icen].x;
  cen[1] = GAL[icen].y;
  cen[2] = GAL[icen].z;

  range = 4*GAL[icen].sigmav/100.0*(1+GAL[icen].redshift)/
    sqrt(OMEGA_M*pow(1+GAL[icen].redshift,3.0) + 1-OMEGA_M);
  set = kd_nearest_range(kd, cen, range);
  
  // Set now contains the nearest neighbours within a distance range. Grab their info. 
  
  while( !kd_res_end(set)) {

    local_cnt++;
    // Get index value of the current neighbor      
    pch = (int*)kd_res_item(set, sat);
    j = *pch;
    kd_res_next(set);
    /*
    if(thisTask==0){
      printf("JJ %d %d %f %f %f %f\n",j,icen,GAL[icen].x, GAL[j].x, range,sat[0]);
      fflush(stdout);
    }
    */
    
    // Skip if target galaxy is the same as the central (obviously).
    if(j == icen)continue;

    // skip if the object is more massive than the icen
    if(GAL[j].mstellar>=GAL[icen].mstellar)continue;
    
    // Skip if already assigned to a central.
    if(GAL[j].psat>0.5)continue;
    
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

    // set the background level
    if(GAL[j].color>0.8) 
      bprob = BPROB_RED + (log10(GAL[j].mstellar)-9.5)*BPROB_XRED;
    else
      bprob = BPROB_BLUE + (log10(GAL[j].mstellar)-9.5)*BPROB_XBLUE;
    
    // combine them into the total probability
    p0 = (1 - 1/(1 + prob_ang * prob_rad / bprob));
    
    if(p0 < 0){
      printf("ZERO %e\n",p0);
      p0 = 0;
    }
    if(p0<0.5)continue;
    if(GAL[j].igrp)continue;
    
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
  
  END_SAT:
  //printf("NCNT %d %d %d\n",thisTask,icen,local_cnt);
  
#pragma omp barrier
  dz =0;
  }
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
