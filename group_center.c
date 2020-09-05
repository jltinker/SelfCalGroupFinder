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

// global for testing luminosty thresholds
float LUM_LIMIT=0;

int group_center(int icen0, void *kd)
{
  float range, rhalo, theta;
  int j, i, n, icen, nmem, nn;
  void *set;
  int *pch;
  double cen[3];
  double sat[3];
  int list[10000], toplist[10];
  int NUM_TEST=2; //how many galaxies to test?

  static int *llist, *vlist, flag=1, *id;
  static float *lumtot, *vsat, *lum;

  // allocate arrays for sorting members
  if(flag)
    {
      llist = ivector(1,NUM_TEST*10);
      vlist = ivector(1,NUM_TEST*10);
      lumtot = vector(1,NUM_TEST*10);
      vsat = vector(1,NUM_TEST*10);
      id = ivector(1,1000);
      lum = vector(1,1000);
      flag = 0;
    }
  
  // allocate space for the memebers (with buffer)
  
  
  // find all members of the group
  nmem = 0;
  for(i=0;i<NGAL;++i)
    {
      if(GAL[i].igrp==icen0 && GAL[i].mstellar>LUM_LIMIT)
	{
	  nmem++;	  
	  id[nmem] = i;
	  lum[nmem] = -GAL[i].mstellar;
	}
    }

  if(nmem==0)return -1;
  // sort to get brightest galaxies
  sort2(nmem, lum, id);

  if(nmem<3)return id[1];
  //if(GAL[icen0].mass>1.0E14)NUM_TEST=10;
  if(nmem<NUM_TEST)NUM_TEST=nmem;

  
  
  // make range only 2sigma_v
  range = 2*GAL[icen0].sigmav/100.0*(1+GAL[icen0].redshift)/
    sqrt(OMEGA_M*pow(1+GAL[icen0].redshift,3.0) + 1-OMEGA_M);

  // get the halo radius
  rhalo = GAL[icen0].theta;

  //printf("search %f %f %f %e\n",range,rhalo,GAL[icen0].sigmav,
  //	 GAL[icen0].rco);
  
  for(i=1;i<=NUM_TEST;++i)
    {
      icen = id[i];
      cen[0] = GAL[icen].x;
      cen[1] = GAL[icen].y;
      cen[2] = GAL[icen].z;

      set = kd_nearest_range(kd, cen, range);
      lumtot[i] = 0;
      vsat[i] = 0;
      llist[i] = icen;
      vlist[i] = icen;
      n= 0;
      nn= 0;
      while( !kd_res_end(set)) {
	
	// Get index value of the current neighbor      
	pch = (int*)kd_res_item(set, sat);
	j = *pch;
	kd_res_next(set);
	theta = angular_separation(GAL[icen].ra,GAL[icen].dec,
				   GAL[j].ra,GAL[j].dec);
	if(theta > rhalo){
	  continue;
	}
	if(j==icen)continue;
	n++;
	lumtot[i] -= pow(GAL[j].mstellar,0.25)/pow(theta + rhalo/10,0.25);
	//vsat[i] += (GAL[j].redshift-GAL[icen].redshift)*SPEED_OF_LIGHT;
      }
      //vsat[i] /= n;
      //vsat[i] = fabs(vsat[i]);
      //printf(" %d %e %e %d %d %d\n",i,lum[i],lumtot[i], n, icen, icen0);
  }

  // who has biggest lumtot?
  sort2(NUM_TEST,lumtot,llist);
  // who has smallest vsat?
  //sort2(NUM_TEST,vsat,vlist);

  // return the one with highest lumtot
  return(llist[1]);
}


void test_centering(void *kd)
{
  FILE *fp;
  int n, i, j, cnt, cnt2, flag,nhalo, icen, flag2;
  float x1, *mass;
  long long int ii, *upid;
  int *listid;


  LUM_LIMIT=-0.4*(-17-4.65);
  fprintf(stderr,"Testing centering... %f\n",LUM_LIMIT);
  LUM_LIMIT = pow(10.0,LUM_LIMIT);
  
  // allocate memory
  listid = calloc(NGAL,sizeof(int));
  upid = calloc(NGAL,sizeof(long long int));
  mass = calloc(NGAL,sizeof(float));
  
  // for each halo, assign the true members to the halo and do the testing procedure

  // first read in the true halo information
  fp = openfile("sham_haloinfo_test5.rdz");
  if(filesize(fp)!=NGAL) {
    fprintf(stderr,"ERROR: file size mismatch with haloinfo file %d %d\n",
	   NGAL, filesize(fp));
  }
  for(i=0;i<NGAL;++i)
    fscanf(fp,"%f %ld %ld %ld %d",&mass[i],&ii,&upid[i],&ii,&listid[i]);
  fclose(fp);

  // go through. send any halo with M>1e13
  cnt = cnt2= nhalo = 0;
  for(i=0;i<NGAL;++i)
    GAL[i].igrp = -1;
  
  for(i=0;i<NGAL;++i)
    {
      //printf("%ld %e\n",upid[i],mass[i]);
      if(upid[i]>=0)continue;
      if(mass[i]<12.5)continue;
      nhalo++;
      //first, set the central
      GAL[i].mass = pow(10.0,mass[i]);
      GAL[i].rad = pow(3*GAL[i].mass/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),THIRD);
      GAL[i].theta = GAL[i].rad/GAL[i].rco;
      GAL[i].sigmav = sqrt(BIG_G*GAL[i].mass/2.0/GAL[i].rad*(1+GAL[i].redshift));
      GAL[i].igrp = i;
      GAL[i].nsat = 0;
      //printf("%e %e %e %e\n",mass[i],GAL[i].ra, GAL[i].dec, GAL[i].redshift);
      flag = 0;
      for(j=0;j<NGAL;++j)
	if(listid[j]==i) {
	  GAL[j].igrp=i;
	  GAL[i].nsat++;
	  // let's check the frequency of Lsat>Lcen
	  if(GAL[j].mstellar>GAL[i].mstellar)flag = 1;
	  // printf("%e %e %e %e\n",mass[i],GAL[j].ra, GAL[j].dec, GAL[j].redshift);
	}      
      if(flag)cnt++;
      // test centering...
      //printf("flag = %d\n",flag);
      icen = group_center(i,kd);
      flag2 = 0;
      if(icen!=i){ cnt2++; flag2=1; }
      // let's check the results
      //if(nhalo==3)exit(0);
      // print out for plots
      printf("CEN %e %d %d %d\n",mass[i],flag, flag2,i);
    }
  printf("%f\n",cnt*1.0/nhalo);
  printf("%f\n",cnt2*1.0/nhalo);
  exit(0);
}

int iterative_center(int ii)
{
  int i,j,n,ntot,i1,i2,ncur, NUMLUM=2, iii, ncnt, ibrightest;
  int skipflag[10000], localid[1000], lumrank[1000], zrank[1000];
  float ra, dec, w, zz, dz, POW;
  static float *theta, *lum, redshift[1000];
  static int first=1, *id, *galid;

  POW = 2;
  if(first)
    {
      lum = vector(1,10000);
      theta = vector(1,10000);
      id = ivector(1,10000);
      galid = ivector(0,10000);
      first = 0;
    }
  
  ntot = GAL[ii].nsat;
  n = GAL[ii].nsat;
  j = ii;

  if(ntot>10)NUMLUM=3;
  if(ntot>20)NUMLUM=5;

  //NUMLUM = 2;
  
  j = ii;
  for(i=0;i<n;++i)
    {
      localid[i] = i;
      galid[i] = j;
      skipflag[i] = 0;
      lum[i] = -GAL[j].mstellar;
      redshift[i] = GAL[j].redshift;
      j = GAL[j].next;
    }

  sort2(n,&lum[-1],&localid[-1]);
  for(i=0;i<n;++i)
    lumrank[localid[i]] = i+1;
  ibrightest = localid[0];
  
  j = ii;
  for(i=0;i<n;++i)
    {
      localid[i] = i;
      redshift[i] = GAL[j].redshift;
      j = GAL[j].next;
    }

  sort2(n,&redshift[-1],&localid[-1]);
  for(i=0;i<n;++i)
    zrank[localid[i]] = i+1;

  /*
  for(i=0;i<n;++i)
    printf("%d %d %d %d %e %d %e\n",i,galid[i],localid[i],lumrank[i],GAL[galid[i]].mstellar,
	   zrank[i],GAL[galid[i]].redshift);
  printf("%d %d %d %e\n",ibrightest, zrank[ibrightest], lumrank[0]-1, galid[lumrank[0]-1], GAL[galid[lumrank[0]-1]].mstellar);
  exit(0);
  */
  //printf("ZRANK %d %d %d\n",n, zrank[0], zrank[ibrightest]);

  // try this: if zrank is <=0.2 take second-brightest
  //if(zrank[ibrightest]*1.0/n<0.2)return galid[lumrank[1]-1];
  //return galid[lumrank[0]-1];
  
  while(n>2) {
    // get center of light for all non-skip memebers
    ra = dec = w = 0;
    j = ii;
    ncnt = 0;
    for(i=0;i<ntot;++i)
      {
	if(skipflag[i])continue;
	if(lumrank[i]<=NUMLUM) { ncnt++; iii=j; }
	ra = ra + pow(GAL[j].mstellar,POW)*GAL[j].ra;
	dec = dec + pow(GAL[j].mstellar,POW)*GAL[j].ra;
	zz = zz + pow(GAL[j].mstellar,POW)*GAL[j].rco;
	w = w + pow(GAL[j].mstellar,POW);
	j = GAL[j].next;
      }
    ra = ra/w;
    dec = dec/w;
    zz = zz/w;

    /* if we only have 1 galaxy within the
     * luminosity threshold, return that galaxy.
     */
    if(ncnt==1)
      return iii;//galid[iii];
    
    // get distance to the center
    j = ii;
    ncur = 0;
    for(i=0;i<ntot;++i)
      {
	if(skipflag[i])continue;
	ncur++;
	id[ncur] = i;
	theta[ncur] = angular_separation(ra, dec, GAL[j].ra, GAL[j].dec)*zz;
	dz = (GAL[j].rco-zz)*0;
	theta[ncur] = theta[ncur]*theta[ncur] + dz*dz;
	j = GAL[j].next;
      }
    
    // sort to get closest;
    sort2(ncur, theta, id);
    // remove the last
    skipflag[id[ncur]]=1;
    n--;
  }

  // now that we only have 2 galaxies, take brightest
  i1 = id[1];
  i2 = id[2];
  if(GAL[galid[11]].mstellar>GAL[galid[i2]].mstellar)return galid[i1];
  return galid[i2];
						       
}
