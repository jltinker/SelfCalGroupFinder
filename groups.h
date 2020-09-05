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


/* Structure definition for galaxies.
 */
extern struct galaxy {
  float x,y,z;
  float ra, dec, redshift;
  float rco;
  float luminosity,
    magnitude,
    mstellar,
    psat,
    color,
    propx,
    propx2,
    weight,
    vmax;
  int igrp;
  int id;
  int listid;
  int next;
  int grp_rank;

  // halo properties  
  float mass,
    theta,
    rad,
    sigmav,
    mtot,
    nsat;
} *GAL;

// other global variables
extern int NGAL;

/* Structure for the halos in the simulations
 */
extern struct halo {
  float x,y,z,vx,vy,vz,mass,lsat;
  int nsat;
} *HALO;
extern int NHALO;

extern int OUTPUT;
extern int FLUXLIM;
extern int STELLAR_MASS;
extern int SECOND_PARAMETER;
extern float FRAC_AREA;
extern float MAXREDSHIFT;
extern float MINREDSHIFT;
extern float GALAXY_DENSITY;

/* Imported functions from numerical recipes 
 */
float qromo(float (*func)(float), float a, float b,
       float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);
float qtrap(float (*func)(float), float a, float b);
float gasdev(long *idum);

/* other functions shared by multiple files
 */
float distance_redshift(float z);
float density2host_halo_zbins(float z);
float density2host_halo_zbins3(float z, float vmax);
float density2halo(float galaxy_density);
float density2host_halo(float galaxy_density);
void lsat_model(void);
int search(int n, float *x, float val);
void test_centering(void *kd);
int group_center(int icen0, void *kd);
float angular_separation(float a1, float d1, float a2, float d2);
void test_fof(void *kd);


