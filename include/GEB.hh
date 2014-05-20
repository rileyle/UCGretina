#ifndef GEB_h
#define GEB_h 1

#define GEB_TYPE_DECOMP         1
#define GEB_TYPE_S800PHYSDATA   9
#define GEB_TYPE_G4SIM         11
#define S800PHYSDATA_TYPETAG  0xABCD1234

#define MAX_SEGS 8              /* max. number of segments to take in events */
#define MAX_INTPTS (2*MAX_SEGS) /* max. number of interaction points */
#ifdef HIGHMULT
#define MAX_SIM_GAMMAS 40       /* max. simulated gammas per event */
#else
#define MAX_SIM_GAMMAS 10       /* max. simulated gammas per event */
#endif

typedef struct gebData
{
  int type;                     /* type of data following */
  int length;
  long long int timestamp;
} GEBDATA;

typedef union {
  GEBDATA gd;
  unsigned int header[4];
} GEBHEADER;

//interaction points
typedef struct ip{
  float x, y, z;
  float e;
  int seg;
  float seg_ener;
} IP;

typedef struct crys_ips_abcd5678 {
  int type;          /* defined as abcd5678 */
  int crystal_id;
  int num;           /* # of int pts from decomp, or # of nets on decomp error */
  float tot_e;       /* dnl corrected */
  int core_e[4];     /* 4 raw core energies from FPGA filter (no shift) */
  long long int timestamp;
  long long trig_time;    /* not yet impl */
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  float prestep;    /* avg trace value before step */
  float poststep;   /* avg trace value following step */
  int pad;          /* non-0 on decomp error, value gives error type */
  IP ips[MAX_INTPTS];
} CRYS_IPS;

typedef struct g4sim_emitted_gamma{
  float e;
  float x, y, z;
  float phi, theta;
  float beta;
} EG;

typedef struct g4sim_abcd1234 {
  int type;          /* defined as abcd1234 */
  int num;           /* # of emitted gammas */
  int full;          /* is full energy */
  EG gammas[MAX_SIM_GAMMAS];
} G4SIM_EGS;

typedef struct S800_physicsdata {
  int32_t type;    /* defined abcd1234 for indicating this version */
  float crdc1_x;   /* Crdc x/y positions in mm */
  float crdc1_y;
  float crdc2_x;
  float crdc2_y;
  float ic_sum;    /* ion chamber energy loss         */
  float tof_xfp;   /* TOF scintillator after A1900    */
  float tof_obj;   /* TOF scintillator in object box  */
  float rf;        /* Cyclotron RF for TOF            */ 
  int32_t trigger; /* Trigger register bit pattern    */
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /* from here corrected values extracted from data above */ 
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  float ic_de;
  /* TOF values with TOF correction applied (from afp/crdc x) */
  float tof_xfpe1;
  float tof_obje1;
  float tof_rfe1;
  /* Trajectory information at target position calculated from 
     a map and afp/bfp/xfp/yfp. New map and you need to re-calc */
  float ata; /* dispersive angle        */
  float bta; /* non-dispersive angle    */
  float dta; /* dT/T T:kinetic energy   */
  float yta; /* non-dispersive position */
} S800_PHYSICSDATA;

#endif
