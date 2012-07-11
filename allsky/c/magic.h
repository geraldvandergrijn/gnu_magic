/* float interpolate(float ddeg,float ylat,float xlon); */ 
/*const int naod=10,nssa=3,ngg=2,nh2o=18,no3=3 */
#define daod 10
#define dssa 3
#define dgg 2
#define dh2o 18
#define do3 3
#define xd 360
#define yd 180     /* only needed for cavesal !!*/
#define SIZE 120

float tttt[360][180];

/* struct and function to read the aerosol climatology -->*/
void raeroclim(int dx, int dy,int month, char aeroclim[]);

struct aerosol {
  float **aod; /** the aerosol optical depth **/
  float **ssa; /** the single scattering albedo **/
  float **gg;  /** the assymetry parameter **/ 
  float *lat;  /* the latitude */
  float *lon; /* the longitude */ 
} aero ;


/* struct and function to read the water vapour climatology -->*/
void rhclim(int dx, int dy, int month, char a[]);

struct absorber {
  float **val; /** matrix pointer to water vapour */
  float *lat;  /* vector pointer to the latitude */
  float *lon; /* vector pointer to the longitude */
} h2o ;

/* struct to ease the implementation of a climatology chosen by the user */ 

void read_something(int x, int y, int m, char a[], int n);

struct dummy {
  float **val; /** matrix pointer to water vapour */
  float *lat;  /* vector pointer to the latitude */
  float *lon; /* vector pointer to the longitude */
}some;

/* struct and function to read the cave surface albedo 
   ---> new for version 0x3 21.05.2008 */
/* float csal[xd][yd]; */

struct cavesal
{
  float csal[xd][yd];    /* the respective surface albedo */
  float dval[xd][yd];  /* d-value */
  float u0[xd][yd];  /* cosine of sZA */
} csal;

struct sarbsal
{
  float **val; /* the SARB sal value */ 
  float *lat; /* the respective latitude */
  float *lon; /* the longitude */  
  unsigned char **lu; /* the land use value */
} sal;

int rcavesal(void);
int rsal(int dx, int dy, int month,char a[]);
int rsarbsal(char a[], char b[], int xsdim, int ysdim);

float normsal(float dval, float u0, float sal);
float corsal(short int igbp, float sza);
/* <--- new for version 0x3 */


/* struct and function to read the clear sky luk up table -->*/

void rclut(int ngg, int nssa, int naod, int nh2o, int no3, char clut[]);

struct aerolut{
  float Im[dgg][dssa][daod];
  float gtau[dgg][dssa][daod];
  float ag[dgg][dssa][daod];
  float btau[dgg][dssa][daod];
  float ab[dgg][dssa][daod];
  float aod[daod];
  float ssa[dssa];
  float gg[dgg];  
};
struct aerolut alut;

struct h2olut
{float  x[dh2o];
  float yg[dh2o];
  float yb[dh2o];
} hlut;

struct o3lut
{float x[do3];
  float yg[do3];
  float yb[do3]; 
}olut;


float interpolate2(float lat,float lon,int xdim,int ydim ,float **var,float x[],float y[]); /* the new interpolation function */
float lint(float xa[],float ya[], int n, float x);
int borders(float xa[], int n,float x);
