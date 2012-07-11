/* ----------------------------------------------------------------------------------------------------------- 

definitions of functions for the calculation of sun earth geometry


Richard Müller 15.02.2009
               included definitions for the new module geoloc 

------------------------------------------------------------------------------------------------------------- */

#define RECLENG (256)                 /* Recordlaenge fuer XPIF und PIF-Format */
#define PI (3.141592654)
#define pi (3.141592654)
#define deg2rad pi/180.0   /* conversion factors from degree to radiance and vice versa */
#define rad2deg 180.0/pi
#define h 42164.0                     /* Height of geostationary Satellite */
#define R_E (6378.155)                /* Erdradius am Aequator in km (DLR)  */
#define R_P (6356.5358)               /* Erdradius am Pol in km (DLR)  */
#define D_SAT (42164)                 /* Distanz Erdmittelpunkt-Satellit in km */
#define RE2H2 -1737122034.796         /* r_e * r_e - h * h */
#define REP2 1.0068137630899          /* r_e * r_e / r_p * r_p */
#define RPE2 0.9932323500             /* r_p * r_p / r_e * r_e */
#define EPSI2 0.00676765              /* (r_e * r_e - r_p * r_p)/(r_e * r_e)*/

struct bildinfo
  { 
   /* VCS Parameters are read from XPIF-HEADER  
    - do not change order of structure elements ! */
    unsigned char  nrbytes;    /* Number of Bytes per Pixel */
    unsigned short nrrecords;  /* Number of Records per Line 1REC=256byte */ 
    unsigned long  nrlines;    /* Number of Lines */
    unsigned long  nrcolumns;  /* Number of Columns with Information */
    unsigned char  add_header; /* Number of additional header Records */ 
    unsigned char  pif_comp;   /* Compression function? 0=no compression */
    unsigned short pif_ahx;    /* Number of additional header records following image data*/ 
    unsigned char  pif_nb;    /* Number of bits per Pixel */
    unsigned char  nav_func;   /* Navigation function code */
    long line_off;             /* Line offset */
    long col_off;              /* Column offset */
    long nav_lres;             /* Line resolution */
    long nav_cres;             /* Column resolution */
    char nav_data[95];         /* specific navigation data, set to 0 if unused */
    unsigned char dat_f;       /* greylevel function */
    char dat_nam[16];          /* data_name ascii */
    char dat_enh[16];          /* default enhancement function */
    char dat_sp[29];           /* spare set to zero */
    short dat_toff;            /* temperature offset */
    char src_obj[16];          /* image source object ascii */
    short src_year;            /* year UT of nominal image time */
    short doy;                 /* day of year of nominal image time  */
    short aq_time;             /* aquisition time = hour*100+minute(UT) of NIT */
    char src_sp[34];           /* spare set to zero */
    char src_version[8];       /* File format version  */
    /* UNI-Oldenburg parameters, not directly given in XPIF-HEADER */
    unsigned long nrcolbig;    /* Number of Columns with black frame */
   
    float fl_space_count; /* Added by CaP */ 
  };

float aperture=18.0,satpos=0.0 ;
  /*the aperture angle of the satellite, needed for the correct calculation of lat and lon 
    if the cloud albedo interface is used 
    the longitudinal satellite position Meteosat is usally at lon=0 ,
    both quantities are overwritten by the input of the config
    file, initialisation is done in order to be able to test older verions..*/
  
int geoloc(int line, int column, int nrlines, double *lat_ptr, double *lon_ptr);


/* Laenge und Breite aus HeaderInformationen nach VCS-Konvention PIF bzw XPIF
   Achtung: gezaehlt werden Line und Column ab 1 (Fortran-Konvention) */


double equation_of_time(int doy, double *gamma, double *f, double *dec);
/*   Equation of Time in Minuten:
     Abweichung der Sonne um 12:00 Mittags von der Suedposition
     siehe: Iqbal, An Introduction to Solar Radiation  */

double cos_zenit(double lat,double lon,double dec,double lon_sun);
		/* Cosinus des Sonnenzenitwinkels */

/* RM May 2008: added functions for the caclulation of the solar azimuth and scatter angle */

double cos_azimuth(double lat,double lon, double dec,double sza,double lon_sun);
 /* cosine of the solar azimuth */

double scatt_angle(double mcossza, double mcosazi, double malpha, double mabeta);
/* the scattering angle alpha=viewing zenith sat, beta=azimuth sat */


double azimuth(double lat,double lon, double dec,double sza,double lon_sun);


void xpifinfo(unsigned char *headline,struct bildinfo *temp);

/* liest aus dem XPIF-Header alle noetigen Informationen 
   und gibt sie als Struktur zurueck zur weiteren Verarbeitung */


int fregrid(float lonrev[], float latrev[], int xrdim, int yrdim, int xdim, int ydim, float **var1, float **var2);


int fregrid2(float lonrev[],float latrev[],int xrdim,int yrdim,int xdim,int ydim,float **var1,float **var2,char cdlfile[]);

int fcdloutCS(float lonrev[], float latrev[], int xrdim, int yrdim,float **var1, float **var2, char cdlfile[]);
