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

/* Local functions
 */
float mean_galaxy_separation(float z);

/* Local globals
 */
int TEST_ID;

void fof(void *kd)
{
  int i,ii,j, iprev, n;
  float plink, zlink, galsep, bp, bz;
  float proj_sep, rad_sep, theta;
  double pt[3], sat[3];
  int ngrp = 0;
  int *grpmem;
  int *grpid;
  void *set;
  int *pch;

  grpmem = calloc(NGAL,sizeof(float));
  grpid = calloc(NGAL,sizeof(float));

  /* initialize the fields
   */
  for(i=0;i<NGAL;++i)
    {
      GAL[i].igrp = -1;
      GAL[i].next = -1;
      GAL[i].mtot = 0;
    }
			
  
  /* Take these linking lengths from
   * Berlind et al 2006 (now modified)
   */
  bp = 0.1;
  bz = 1.0;

  /* Increasing for fun
   */
  //bp = 0.15;
  //bz = 1.5;

  // what is local linking lengths?
  galsep = mean_galaxy_separation(GAL[i].redshift);
  plink = bp*galsep;
  zlink = bz*galsep;
  
  //fprintf(stderr,"0.000");
  for(i=0;i<NGAL;++i)
    {
      i = TEST_ID;
      
      //printf("here %d\n",i);
      // if linked, skip
      if(GAL[i].igrp>=0)continue;
      if(i%1000==1)
	fprintf(stderr,"\b\b\b\b\b%.3f",i*1.0/NGAL); 

      // set up new group
      iprev = i;
      grpmem[ngrp] = 1;
      grpid[ngrp] = i;
      ngrp++;
      GAL[i].igrp = i;
      GAL[i].mtot = GAL[i].mstellar;
      
      // what is local linking lengths?
      galsep = mean_galaxy_separation(GAL[i].redshift);
      plink = bp*galsep;
      zlink = bz*galsep;
      
      // find all neighbors of this galaxy
      pt[0] = GAL[i].x;
      pt[1] = GAL[i].y;
      pt[2] = GAL[i].z;

      set = kd_nearest_range(kd, pt, zlink*2);
      
      while( !kd_res_end(set)) {
	
	// Get index value of the current neighbor      
	pch = (int*)kd_res_item(set, sat);
	j = *pch;
	kd_res_next(set);
	if(j==i)continue;
	if(GAL[j].igrp!=-1)continue;
	theta = angular_separation(GAL[i].ra,GAL[i].dec,
				   GAL[j].ra,GAL[j].dec);
	proj_sep = theta*(GAL[i].rco+GAL[j].rco)/2;
	rad_sep = fabs(GAL[i].rco-GAL[j].rco);
	//if(i==1)printf("X %d %f %f\n",j,proj_sep/plink, rad_sep/zlink);    
	if(proj_sep<=plink && rad_sep<=zlink)
	  {
	    GAL[j].igrp = i;
	    GAL[iprev].next = j;
	    GAL[i].mtot += GAL[j].mstellar;
	    grpmem[ngrp-1]++;
	    iprev = j;
	    //printf("%d %d\n",i,grpmem[ngrp-1]);
	  }
      }
      
      /* now go through the list of all neighbors 
       * to find their neighbors
       */
      ii = GAL[i].next;
      //printf("%d\n",GAL[ii].next);
      while(ii>=0) {

	pt[0] = GAL[ii].x;
	pt[1] = GAL[ii].y;
	pt[2] = GAL[ii].z;
	
	set = kd_nearest_range(kd, pt, zlink*2);
	
	while( !kd_res_end(set)) {
	  
	  // Get index value of the current neighbor      
	  pch = (int*)kd_res_item(set, sat);
	  j = *pch;
	  kd_res_next(set);
	  if(j==ii)continue;
	  if(GAL[j].igrp>=0)continue;
	  theta = angular_separation(GAL[ii].ra,GAL[ii].dec,
				     GAL[j].ra,GAL[j].dec);
	  proj_sep = theta*(GAL[ii].rco+GAL[j].rco)/2;
	  rad_sep = fabs(GAL[ii].rco-GAL[j].rco);
	  if(proj_sep<=plink && rad_sep<=zlink)
	    {
	      GAL[j].igrp = i;
	      GAL[iprev].next = j;
	      GAL[i].mtot += GAL[j].mstellar;
	      grpmem[ngrp-1]++;
	      iprev = j;
	      //printf("%d %d\n",i,grpmem[ngrp-1]);
	    }
	}
	ii = GAL[ii].next;
      }
      GAL[i].nsat = grpmem[ngrp-1];
      return;
      
      /* print out first group
       */
      if(i==1){
	
	j = i;
	n = 1;
	do {	  
	  printf("%d %d %e %e %e\n",n,j,
		 GAL[j].ra*180/PI, GAL[j].dec*180/PI, GAL[j].mstellar);
	  j = GAL[j].next;	  
	  n++;
	} while(j!=-1);
	GAL[i].nsat = n;
	return;
      }
    }

  /* Let's output the groups
   */
  for(i=0;i<ngrp;++i)
    {
      j = grpid[i];
      printf("%d %d %d %e\n",i,j,grpmem[i],GAL[j].mtot);
    }
  exit(0);
}

float mean_galaxy_separation(float z)
{
  static float rg=-1;
  // simple if volume-limited samples
  if(!FLUXLIM)
    {
      if(rg<0)
	rg = pow(GALAXY_DENSITY,-THIRD);
      return rg;
    }
}

void test_fof(void *kd)
{
  FILE *fp;
  int n, i, j, cnt, cnt2, flag,nhalo, icen, flag2, ntrue, nfof, j1;
  float x1, *mass, theta, zlink, plink, galsep, rad_sep, proj_sep;
  long long int ii, *upid;
  int *listid;

  TEST_ID = 1;
  fof(kd);

  // allocate memory
  listid = calloc(NGAL,sizeof(int));
  upid = calloc(NGAL,sizeof(long long int));
  mass = calloc(NGAL,sizeof(float));
  
  // for each halo, assign the true members to the halo and do the testing procedure

  // first read in the true halo information
  fp = openfile("sham_haloinfo_test6.rdz");
  if(filesize(fp)!=NGAL) {
    fprintf(stderr,"ERROR: file size mismatch with haloinfo file %d %d\n",
	   NGAL, filesize(fp));
  }
  for(i=0;i<NGAL;++i)
    fscanf(fp,"%f %ld %ld %ld %d",&mass[i],&ii,&upid[i],&ii,&listid[i]);
  fclose(fp);

  // lets find the rproj and dz of known groups
  ii = 1;
  galsep = mean_galaxy_separation(0);
  zlink = 0.8*galsep;
  plink = 0.14*galsep;
  for(ii=0;ii<NGAL;++ii)
    {
      if(mass[ii]<12.5)continue;
      if(upid[ii]!=-1)continue;

      GAL[ii].mass = pow(10.0,mass[ii]);
      GAL[ii].rad = pow(3*GAL[ii].mass/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),THIRD);
      GAL[ii].theta = GAL[ii].rad/GAL[ii].rco;
      GAL[ii].sigmav = sqrt(BIG_G*GAL[ii].mass/2.0/GAL[ii].rad*(1+GAL[ii].redshift));
      
      ntrue = nfof = 1;
      TEST_ID = ii;
      fof(kd);
      flag = 0;
      for(i=0;i<NGAL;++i)
	{      
	  /*
	  theta = angular_separation(GAL[ii].ra,GAL[ii].dec,
				     GAL[i].ra,GAL[i].dec);
	  proj_sep = theta*(GAL[ii].rco+GAL[i].rco)/2;
	  rad_sep = fabs(GAL[ii].rco-GAL[i].rco);
	  */
	  if(GAL[i].igrp==ii) {
	    nfof++;
	    if(GAL[i].mstellar>GAL[ii].mstellar)flag = 1;
	  }
	  if(listid[i]==ii)ntrue++;
	  //printf("MOO %d %f %f %d\n",i,proj_sep/plink, rad_sep/zlink, GAL[i].igrp);      
	}
      if(nfof<5)continue;
      if(ntrue==0)continue;
      // get the center
      //j=group_center(ii,kd);
      j = 0;
      j1 = iterative_center(ii);
      //ntrue++;
      printf("%f %f %d %d %d %d %d %f %.0f\n",nfof*1./ntrue,
	     GAL[ii].nsat/ntrue, nfof,ii, j,j1, flag ,mass[ii], GAL[ii].nsat);
      fflush(stdout);
    }
  exit(0);
}
