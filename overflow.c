void find_satellites(int icen, void *kd) {
  int j, k, ii, neighbor[100000], i, cnt;
  float dx, dy, dz, theta, prob_ang, vol_corr, prob_rad, grp_lum, p0, range;
  float cenDist, bprob, mtot, nsat;
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
  // get the list of all neighbors to farm out to multithreads
  cnt = 0;
  while( !kd_res_end(set)) {
    pch = (int*)kd_res_item(set, sat);
    neighbor[cnt] = *pch;
    kd_res_next(set);
    cnt++;
  }

  
  mtot = nsat = 0;
#pragma omp parallel private(ii,j,dz,theta,prob_ang,prob_rad,bprob,p0)
  {
    mtot=0;nsat=0;
#pragma omp for  reduction(+:mtot) reduction(+:nsat)
    for(ii=0;ii<cnt;++ii) {
    j = neighbor[ii];
    
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
    if(p0>0.5){
      // this is considered a member of the group
      GAL[j].psat = p0;
      GAL[j].igrp = icen;
      mtot += GAL[j].mstellar;
      nsat ++;
    }
    //GAL[icen].mtot += GAL[j].mstellar;
    //GAL[icen].nsat++;
    }
  }
  GAL[icen].mtot += mtot;
  GAL[icen].nsat = nsat;

  
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
