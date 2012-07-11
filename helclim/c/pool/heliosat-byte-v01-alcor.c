/******/

/* V1x22 based on University of Oldenburg Heliosat cloudindex scheme 
   albedo_and_cloudindex_shadow_msg.c
   MODIFICATIONS:
   RM 10.01.06: added the calculation and output of the surface reflectance in %  
   RM May  2008:   added automatic rho_max calculation.
   RM June 2008:   added automatic calculation of sigma_ground
                   bugfix: generation of reflectance maps now includes month and year in the name
                   changed loops, in order to be consistent with Eumetsat version loops start at k=0 
                   instead of k=1 !!! 
                   minor changes-> see comments within the code 
   RM July 2008:   -impelemented new scatter angle correction for the calculation of rhomax
                   -calculation of satazi and satzen now in main routine instead of msg_rgeo...
                   -dpos can be provided, it replaces sig_ground in the function call
                    sig_ground. Combined with automatic rhomax calculation also sig_ground is
                    calculated automatically, hence  has no longer to be provided within 
                    the heliosat call, instead dpos has to be provided, dpos is the position 
                    of Meteosat relativ to lon=0, West is negative
                    e.g. Meteosat is at 0 degree, hence dpos=0, Meteosat-8 was on 3.4 W, 
                    hence dpos=-3.4. 
                    By intorducin dpos there has been a need to modify the equation for satazi 
                    calculation. The equation has been extended to
                    cover all longitudes and not only Meteosat position at lon=0.
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    It is important to provide Heliosat during the call the longitude 
                    position of Meteosat using argv parameter -s
                   -Added functionality that the reference region changes with dpos in order
                    to get appropriate reference region also fo Meteosat East.
   RM Jan/Feb 2008:-impplemted new routines and modified others for the calculation of the
                    geolocation and GMT in order to be able to use MFG and MSG.
                    The original code was only able to work with MSG.
                   -addded functions and commands to read a byte matrix (MFG) instead of 
                    a short int matrix (MSG). 
                   -implemented comands to transform the MFG image byte information of to 
                    short integer in order to have always the same output format of the cloud 
                    index files.
   RM Mar      2009: modified while(cloudy  && cloud_free_quantity) loop as the condition lead to
                     endless loops on some machines, e.g SUN, even if it make no sense to use other 
                     machines and not LINUX machines some companies do. 

   RM Jul      2009: -fixed GMT bug, RM assumed Nort/South detector scannes noth and south of Earth 
                      in ||  but it scans two subsequent lines in ||, hence bug fix of Carlo has 
                      been implemented 
                     -replaced in program call -z Pramater sig_shadow, as it is not needed,  by MFG.
                      If -z is 1, it is assumed that the input data comes from MFG, dark offset is then
                      set to 4.5 otherwise (MSG) to 51.
                     -geoloc (lat,lon calculation) is used only for MFG, for MSG geomsg is used.
                      todo: In the next version it is foreseen to use geoloc for both satellites again,
                      but as a prerequisite geoloc has to be modified (e.g. provide aperture angle in
                      function call) 
   COMPILATION:      
   compile with  gcc heliosat-vXX.c -lm -m32 -o DIR/myheliosat.exe
*/


#include <stdio.h>             /* Ein- und Ausgabedeklarationen            */
#include <stdlib.h>            /* Umwandlung von Zahlen, Speicherverwaltung*/
#include <stdbool.h>
#include <string.h>            /* Verarbeitung von Zeichenketten           */
#include <ctype.h>             /* Fuer die Funktion tolower()              */
#include <math.h>              /* mathematische Funktionen                 */

#include "./satgeo-32b.h"         /* Geolocation und alles fuer MSG           */

#define CLOUDI_MAX (1.2)       /* binaer Speichern des Cloud Index         */
#define CLOUDI_MIN (-0.2)
#define NO_VALUE   (9999)      /* falls Sensorausfall */
#define BYTE_MAX (1023)        /* binaer Speichern des Cloud Index         */
#define BYTE_MIN (0)
#define HDL (256)              /* Headleng=RECLENG=256 für XPIF */
                               /* Headleng=260 für X260 Oldenburg */

#define maxi(A,B) ((A)>(B) ? (A) : (B))  /* Kernighan/Ritchie:          */ 
#define mini(A,B) ((A)<(B) ? (A) : (B))  /*  Programmieren in C, ANSI C */

#define BILDANZ (32)           /* Zahl der Tage eines Monats + GrundalbedoBild */

double slope, y_intercept;     /* affine transformation of cloud index to short int values */
             
void exit_print_usage(void);
void exit_no_list(char *name);
void exit_no_write(char *name);
void exit_no_read(char *name);
void heapsort(unsigned long n, short int ra[]);
bool file_exists(const char *filename);

short int **smatrix(long nrl, long nrh, long ncl, long nch);
/* allocate a short int matrix with subscript range m[nrl..nrh][ncl..nch] 
   A. Hammer 2005 - concept of Numerical Recipes in C */
void free_smatrix(short int **m, long nrl, long nrh, long ncl, long nch);
/* free a short int matrix allocated by smatrix()
   A. Hammer 2005 - concept of Numerical Recipes in C  */
 
void read_xpif_header(short int xpif_head_length,char *image_name,unsigned char *headline);

short int **read_xpif_to_matrix(short int xpif_head_length,char *image_name,struct bildinfo *imageinfo);

/* RM FEB 2009, added stuff to read in the MFG data ,bmatrix and read_xpif_Bmatrix */
/* RM FEB 2009 -----> */
typedef unsigned char byte;
byte **bmatrix(long nrl, long nrh, long ncl, long nch);
/* allocate a  byte matrix with subscript range m[nrl..nrh][ncl..nch] 
  concept of Numerical Recipes in C */
void free_bmatrix(byte **m, long nrl, long nrh, long ncl, long nch);
byte **read_xpif_to_Bmatrix(short int xpif_head_length,char *image_name,struct bildinfo *imageinfo);
/* <----- RM FEB 2009 */
short int groundreflectivity_and_cloudindex(short int *reflectivities, long zahl, 
   short int sigma_g, short int sigma_s, short int rhomax);
/* calculates ground reflectivity and cloud index 
   A. Hammer 2005 */
short int cloudindex(short int reflectivity,short int ground,short int rhomax,short int shadow);
   
/****************** M A I N ************************************************
 *  
 *
 *    input : Zeitreihe von Satelittenbildern vis006 vis008 oder HRV (xpif)
 *    output: 1 Bodenalbedobild  
 *
 *
 ****************** M A I N ***********************************************/
 

int main(int argc, char *argv[])
{ 
  struct bildinfo geoinfo;        /* alles zum Bildbearbeiten aus dem Header*/ 
  struct bildattribut{
      char name[128];          /* Struktur, um fuer jedes Bild die geo-  */
      int doy;                 /* zu verarbeiten                         */
      double time;
      double gamma; 
      double f;
      double dec;
      double EOT;
      double lon_sun;  
      unsigned char headline[HDL];  /* Recordlaenge fuer "echte" Groesse    */ 
      short int **imagedata;  /* matrix; msg data are assigned to imagedata */
      byte **Bimagedata; /* Matrix, mfg data are assigned to Bimagedata */
  };
  struct bildattribut bildfolge[BILDANZ]; /* Zahl der Listenelemente: Monat = BILDANZ */

  int rho_max=0;           /* maximale Reflektivitaet bei Met-8 Bildern   */
                             /* Cosinuskorrektur                            */
  double GMT;                /* Abscannzeitpunkt                            */
  int quiet=0;               /* Bildschirminformationen ?                   */
  int accm;                  /* parameter for MFG GMT calculation           */
  int lin, col;              /* Zaehlvariablen fuer Zeile und Spalte        */
  int n=0;                   /* Anzahl der verarbeiteten  Bilder            */
  int decide;                /* ausserhalb der Erde oder nicht              */
  double lat,lon;            /* latitude and longitude in radians !         */
  double corr;               /* normalisierte Radianz                       */
  double coszen;             /* Cosinus des Zenitwinkels                    */
  /* ---> SELFCALIBRATION: declarations */
  /* RM May 2008: added file pointer frhomax,fstatrho */
  FILE *out_image,*list,*f,*frhomax,*fstatrho;
  /* RM May 2008: fper and kend1 is needed for the caclulation of the percentile */
  short int *fper = NULL;  
  int kend1=0;
  /* RM May 2008: added averaged angles used for the caclculation of the scatter angle scata */
  float msza=0;  /* average of the solar zenith angle over time and space within the calibration region */
  float mazi=0;  /* average of the solar azimuth angle over time and space within the calibration region */
  float msatazi=0; /* average of the sat azimuth angle .... */
  float msatzen=0; /* average of the satellite zenith angle */
  float scata=0;   /* the mean scattering angle */
  double satazi, satzen;     /* geometry of MSG, satellite zenith and azimuth  */
  float r2d=180.0/3.1416;
  float esd;      /* variable for the geometric correction of rho_max */
  float dpos=0.0; /* dpos is the position of Meteosat relativ to lon=0, West is negative, East positive */
  /* <--- SELFCALIBRATION declarations */

                              /* Dateinamen von Bildern                      */ 
  char out_imagefile[160],groundpathname[128],
       cloudpathname[128],bild[128],bildpathname[128];            
  short int **out;           /* Matrix fuer die Ausgabewerte       	    */
  short int reflectivity_vector[BILDANZ]; 
  short int min, hour;

  char *cp;                         /* temporaere Werte                     */
  char zeitname[4];                 /* temporaerer ASCII-String der Zeit    */
         
  long zahl,anzahl;
  double cloudindex ;
  double db_hours_per_line; /* added in order to caclulate the timeoffset */
  char *channel="MFGVIS"; /* added, currently hard wired, shoul be provided via heliosat call once the thermal irradiance is
			included */
  short int MFG=1; /* added, currently hard wired, should be provided via heliosat call, once the routine is successfully tested for 
                   MFG and MSG !! */
  int i,k,laenge;
  char quite='\0';
  short int sig_ground, sig_shadow=0; /* width of distribution of cloudfree pixels 
                                        and shaded pixels */
  /* sig_shadow is zero, the shadow module will not be used, it is too slow and probably not needed if 
  the fuzzy logic Heliosat version is used*/

  short int(temp);
  float doffset=0; /* the dark/space offset of the instrument */

  slope=(BYTE_MAX-BYTE_MIN)/(CLOUDI_MAX-CLOUDI_MIN);
  y_intercept=BYTE_MAX-CLOUDI_MAX*slope;                   /* DEFINE */
  *out_imagefile = '\0';
  for(i = 0;i <= BILDANZ;i++){*bildfolge[i].name = '\0';}
  

/*** optionale Argumente von der command line lesen... ***/

  if (argc<4)             /* wenn nicht mindestens zwei opt.Parameter.. */
    exit_print_usage();   /* .. dann schliesse das Programm */

  for (i=1; i<argc; i++)  /* fuer alle Argumente..                          */
  {
    cp=argv[i];
    if (*cp!='-' && *cp!='/') /* pruefe ob Argument mit '-' oder '/' beginnt */
    {
      fprintf(stderr,"bodenalbedo: no option specified (%s)\n",cp);
      exit_print_usage();
    }
    cp++;
    if (*(cp+1))            /* pruefe ob an uebernaechster Stelle ein Blank  */
    {
      fprintf(stderr,"bodenalbedo: missing blank. Usage: -%c <arg>\n",*cp);
      exit(1);
    }
    *cp=(char) tolower((int) *cp); /* es wird nicht zwischen Klein- und      *
				    * Grossschreibung unterschieden          */

    if (*cp=='l') list=fopen(argv[++i],"rt");
    else if (*cp=='b') strcpy(bildpathname,argv[++i]);
    else if (*cp=='c') strcpy(cloudpathname,argv[++i]);
    else if (*cp=='g') strcpy(groundpathname,argv[++i]);
    /* ---> SELFCALIBRATION dpos instead of sig_ground !!! */
    else if (*cp=='s') dpos=atof(argv[++i]); /*dpos instead of sig_ground, see  program header for details */
    else if (*cp=='z') MFG=atoi(argv[++i]); 
    /* MFG instead of sig_shadow MFG=1: MFG data expected MFG=0: MSG data */ 
    else if (*cp=='q') quite='q';
    else {
      fprintf(stderr,"bodenalbedo: unknown option: %s\n",cp);
      exit_print_usage();
    }
  }
  if (dpos != 0.0) 
  {printf("Meteosat position not at longitude=0 ! Please check -s in the heliosat call & set it to 0 if Met is at lon=0 \n");}

  if(MFG==1)
  {printf("assume MFG image data as input !! \n please check if you realy use MFG: -z has to be 0 for MSG \n");}
  else{
      printf("assume MSG image data as input !! \n please check if you realy use MSG: -z has to be 1 for MFG \n");
  }


 /* oeffnen aller Dateien, die in Liste stehen, Header auswerten, Matrizen zuweisen */
 zahl=0;
 while(!feof(list) && anzahl < BILDANZ )
       /* CaP - inserting the limit of BILDANZ to avoid coredumps*/
       /* ||  anzahl < BILDANZ )        */
 { 
   anzahl=zahl;
   fscanf(list,"%s\r",bildfolge[anzahl].name);
   strcpy(bild,bildpathname);
   strcat(bild,"/");
   strcat(bild,bildfolge[anzahl].name);
   
   if(!quite)printf("%d %s\n",anzahl,bild);
   
   read_xpif_header(HDL,bild,&bildfolge[anzahl].headline[0]); /* store header for output XPIF-images */
   
   
    if(HDL==256)
    {
     bildfolge[anzahl].headline[11]=(unsigned char) 0;          /* no additional headers wanted !      */
    }
    if (MFG==1)
    { 
      bildfolge[anzahl].Bimagedata=read_xpif_to_Bmatrix(HDL,bild,&geoinfo);
      bildfolge[anzahl].imagedata=smatrix(0,geoinfo.nrlines-1,0,geoinfo.nrcolumns-1);
      /* RM 17.07.09: dark offset and snanning time is different for MFG and MSG */
      db_hours_per_line=(25./geoinfo.nrlines/60.);
      doffset=4.5*4;
      for (lin = 0;lin < geoinfo.nrlines;lin++) {  /* change type and rotate the image, afterwards free bmatrix */
	  for(col = 0;col < geoinfo.nrcolumns;col++){ 
	      /* temp=(short int)bildfolge[anzahl].Bimagedata[geoinfo.nrlines-lin-1][geoinfo.nrcolumns-col-1]; */
	      /* correct if old openMTP converter is used */
	  bildfolge[anzahl].imagedata[lin][col]=4*(short int)bildfolge[anzahl].Bimagedata[geoinfo.nrlines-lin-1][geoinfo.nrcolumns-col-1]; 
	  /* correct if new openMTP to XPIF converter is used */
          /* bildfolge[anzahl].imagedata[lin][col]=4*(short int)bildfolge[anzahl].Bimagedata[lin][col]; */

	  /* multiplication with factor 4 needed in order to have only one BYTE MAX definition */  
           /* (short int)bildfolge[anzahl].Bimagedata[lin][col]; */
	      /* (short int)bildfolge[anzahl].Bimagedata[geoinfo.nrlines-lin-1][geoinfo.nrcolumns-col-1]; */
            /*       printf("i,j,temp= %d %d %d",lin, col, temp); */
	  }
      }
      free_bmatrix(bildfolge[anzahl].Bimagedata,0,geoinfo.nrlines-1,0,geoinfo.nrcolumns-1); 
	/* bildfolge[anzahl].Bimagedata=read_xpif_to_Bmatrix(HDL,bild,&geoinfo); */ 
    } 
    else{
	bildfolge[anzahl].imagedata=read_xpif_to_matrix(HDL,bild,&geoinfo);
	db_hours_per_line=(12./geoinfo.nrlines/60.);
	doffset=51;
    }
    printf("Iam here 1 %d %d \n",geoinfo.nrlines, geoinfo.nrcolumns); 
    bildfolge[anzahl].doy=geoinfo.doy;
    hour=geoinfo.aq_time /100;
    min=fmod(geoinfo.aq_time,100);
    bildfolge[anzahl].time=hour+min/60.0;
    bildfolge[anzahl].EOT=equation_of_time(bildfolge[anzahl].doy,&bildfolge[anzahl].gamma,&bildfolge[anzahl].f,&bildfolge[anzahl].dec);
    zahl++; 
  } 
 n=anzahl;    
 printf ("n= %d \n",n);
 /*  printf("doy is %04d %04d \n",geoinfo.doy, geoinfo.aq_time );   
    printf("nrs is %d %d %d %d \n",geoinfo.nrlines,geoinfo.nrcolbig,geoinfo.nrcolumns, geoinfo.nrlines );
    printf("nrs is %d %d %d %d \n",geoinfo.line_off,geoinfo.col_off,geoinfo.nav_lres, geoinfo.nav_cres );  */
  
 /* orig out=smatrix(0,geoinfo.nrlines-1,0,geoinfo.nrcolbig-1); RM 24.02->*/
  out=smatrix(0,geoinfo.nrlines-1,0,geoinfo.nrcolumns-1);
 
 /***  Normalize images   ---------------------------------------- ***/
 
 /* allocate memory for fper 42610*/
/* fper = malloc(45000 * (n+1) * sizeof(short int));
 if(fper == NULL)
 {
     fprintf(stderr, "not enough memory for array fper, needed for automatic calibration \n");
     exit;
     } */
  /* CaP - allocating only 2 bytes for now. Using realloc later in the loop below (line 719).
The former memory allocation (previous line) causes severe memory errors and 
program crashes */
 
  /* --> SELFCALIBRATION: initial allocation of field needed for the percentiles */
  fper = malloc(sizeof(short int));
  /* <-- */
  printf("Iam here 2 \n");
  accm=0;
 for (lin = 0;lin < geoinfo.nrlines;lin++)   /* fuer alle Zeilen */
 { 
/* 12.10.09 included new lines for the calculation of GMT. Now MFG and MSg can be processed */
 /* old, orig for MSG GMT=bildfolge[anzahl].time+timeoffset(geoinfo.nav_lres,geoinfo.line_off,lin,15); */
     /* GMT : Abscannzeitpunkt dieser Zeile bestimmen */
      /* 17.07.09 fixed GMT bug, accm needed to consider that two subsequent lines are scanned in || */
     accm=accm+abs((lin-1)%2); 
     GMT=bildfolge[anzahl].time+(double)lin*db_hours_per_line; 
     /* printf("now here vv \n");*/
    if ((strcmp(channel,"MFGVIS")==0) && MFG==1)
     {  /* CaP 29.04.2009 - in the vis imagery adjacent lines are scanned in couples */
	 GMT=bildfolge[anzahl].time+(double)(accm-1)*db_hours_per_line;
	 /*  printf("GMT,accm= %f %d \n",GMT ,accm); */
     }
     /* orig  GMT=bildfolge[anzahl].time+timeoffset(geoinfo.nav_lres,geoinfo.line_off,lin,15); */  
     /* GMT=bildfolge[anzahl].time+timeoffset(geoinfo.nav_lres,geoinfo.line_off,lin,15);  */

     for (k = 0 ;k <= n;k++)  /* fuer alle Bilder */
     {
	 /* fuer Abscannzeit dieser Zeile und jeden Tag den Sonnenstand bestimmen */
	 bildfolge[k].lon_sun=(12.0-GMT-(bildfolge[k].EOT/60.0))*15*(PI/180.0);
     }
     /* printf("ere re re %d \n",k); */
     for(col = 0;col < geoinfo.nrcolumns;col++)  /* fuer alle Spalten */
     { 
         if(MFG==1){
	     decide=geoloc(lin,col,geoinfo.nrlines,&lat,&lon);
         }
         else{
	decide=geomsg(geoinfo.nav_cres,geoinfo.nav_lres,geoinfo.col_off,geoinfo.line_off,lin,col,&lat,&lon);
	 }
         /* if (decide == 0){
         printf("the geoloc i,j,lat,lon %d %d %f %f \n", lin, col,lat*180/pi,lon*180/pi);
         } */
	     /* orig decide=geomsg(geoinfo.nav_cres,geoinfo.nav_lres,geoinfo.col_off,geoinfo.line_off,lin,col,&lat,&lon); */

	 /* ---> SELFCALIBRATION: calculation of additional angles */
         /* RM Jul 2008: added equations for the calculation of satellite zenith and azimuth 
		  !! it is assumed that the satellite is at lat=lon=0
		  h=42164.0 height of geostationary Satellite */
	 satazi = atan(sin(lon-dpos/r2d)/tan(lat)); /* rm 072008 modified formular, covers all positions of Meteosat now */
	 /* old satzen = fabs(atan(sin(lat) * h/(cos(lat) * h - R_E))); */
    satzen=90/r2d-atan((cos(dpos/r2d-lon)*cos(lat)-0.1512)/sqrt(1-pow(cos(dpos/r2d-lon),2)*pow(cos(lat),2)));     /* <--- SELFCALIBRATION: calculation of additional angles */         
	 if (decide==0)
	 {
	     for(k = 0 ;k <= n;k++)  /* fuer alle Bilder */
	     {
		 corr=0;
		 coszen=cos_zenit(lat,lon,bildfolge[k].dec,bildfolge[k].lon_sun);
		 if(coszen>0.02)
		     corr = (int)((bildfolge[k].imagedata[lin][col]-doffset)/coszen); /* orig 51 instead of doffset */
		 corr = mini(corr,BYTE_MAX);
		 corr = maxi(corr,BYTE_MIN); 
		 bildfolge[k].imagedata[lin][col] = corr;
		 /* RM 20.11.07:  automatic calculation of rho_max starts, sort vector content according 
		    to the value, in ascending order using heapsort, Numerical Recipies in C 
		    the 95 percentil * 1.1 is than the desired rhomax value */
		 /* ---> SELFCALIBRATION: Calculate the percentile and additional angles
                         needed for the calculation of the scatter angle */
	  if(lat > -0.92 && lat < -0.8 && lon  < (-0.01 + dpos/r2d) && lon > (-0.262+dpos/r2d)) /* -0.436*/     
		/* RM 072008: added +dpos in order to shift the reference region for Meteosat-East
		dpos=position of Meteosat relative to lon=0 in degree, west is negative */  
               {
		
		fper[kend1]=(short int)bildfolge[k].imagedata[lin][col]; /* assign values within region to fper */
		kend1=kend1+1;
		/* calculate the mean angles as basis for the caclulation of the scatter angle, 
		   needed for the geometric correction of rhomax */
		fper = realloc(fper,(kend1+1)*sizeof(short int));
		msza=((kend1-1)*msza+acos(coszen))/kend1; 
		mazi=((kend1-1)*mazi+azimuth(lat,lon,bildfolge[k].dec,acos(coszen),bildfolge[k].lon_sun))/kend1;
		msatazi=((kend1-1)*msatazi+satazi)/kend1;
		msatzen=((kend1-1)*msatzen+satzen)/kend1;
	        }
            /* <--- end of SELFCALIBRATION: Calculate the percentile and additional angles
                         needed for the calculation of the scatter angle */
	     }
	 }
         if (decide==1){
	     out[lin][col] = (int) 0;
	     for(k=1;k<=n;k++){
		 /*  if(geoinfo.nrlines<200){doffset=doffset+ bildfolge[15].imagedata[lin][col];} */
		     bildfolge[k].imagedata[lin][col] = 0;
	     }
	 }
     }
 }
 printf("doffset= %f \n", doffset);
 printf("now I am here \n"); 

 /* ---> SELFCALIBRATION: sorting the fper field and assignment of the percentiles */
 heapsort(kend1, fper); /* order the counts */
 int p1=(int)(0.1*(kend1+1));
 int p3=(int)(0.3*(kend1+1));
 int p5=(int)(0.5*(kend1+1));
 int p7=(int)(0.7*(kend1+1));
 int p8=(int)(0.8*(kend1+1));
 int p9=(int)(0.9*(kend1+1));
 int pm=(int)(0.95*(kend1+1));
 /* <--- SELFCALIBRATION: sorting the fper field and assignment of the percentiles */
 /* ---> SELFCALIBRATION: empirical correction of anisotropy for rho_max */
 /* RM 28.05.08 empirical correction for scatter angle effect, so far formula only proven for 3 - 33 degree,
    RM 30.07.08 implemented new formula, correction works from 3-70 degree (scata), but probably only for 
             the reference regions and not everywhere, above scattering angle of 32 the uncertainty of
             the geometric corrections is significnat higher */
 scata=scatt_angle(msza, mazi, msatzen, msatazi);
 /* added new correction for rho_max, seems that the dominant effect comes from the
   sza dependent change of the cloud albedo, the clouds seems to reflect pretty well
   isotropic for a given SZA and rel AZA lower than 60 degree
   the correction formula results from a linear fit (gnuplot) applied to a 
   2.5 year time series for all slots with SZA>60 */
   if (msza*r2d > 80.0 && msza*r2d <0.0) esd=1.0;  /* gcor =1 for SZA>80 and <0*/
   else
   {esd=1+0.0017*(45.0-msza*r2d); }

 /* <--- SELFCALIBRATION: empirical correction of anisotropic cloud reflection for rho_max */
 /* ---> SELFCALIBRATION: Open and write the output files for information  related to  the 
         self calibration approach */
  fstatrho=fopen("rhostat.out","a") ;
 fprintf(fstatrho, "# kend1= %d doy= %d hour= %d min = %d \n",kend1,geoinfo.doy, hour,min);
 fprintf(fstatrho, " gcor= %f msza= %f mazi= %f satazi= %f satzen= %f scattangle= %f rhomax %f pm= %d \n", esd, msza*r2d, (mazi)*r2d, msatazi*r2d, msatzen*r2d, scata*r2d, fper[pm]*esd*1.1, fper[pm]);
 fprintf(fstatrho, "# ESD= %f rho(10,30,50,70,80,90,95)= %d  %d %d %d %d %d  %d rhomax %f \n",esd,fper[p1],fper[p3],fper[p5],fper[p7], fper[p8],fper[p9],fper[pm], fper[pm]*esd*1.1);
   if (hour==13 && min==0)
   {
       frhomax=fopen("rhomax.out","w");  /* RM 072008 changed r+ to w, r+ does not work properly on my machine */
       rho_max=fper[pm]*esd*1.1;
       fprintf(frhomax, "%d \n", rho_max);
       fprintf(fstatrho,"## The noon rho_may is %d \n",rho_max );
       fclose(frhomax); 
   }
   else
   {
       /* if((frhomax=fopen("rhomax.out","r")) == NULL) */
       if(file_exists("rhomax.out")) 
       {
           frhomax=fopen("rhomax.out","r");
	   fscanf(frhomax, "%d", &rho_max);	   
	   fprintf(fstatrho,"# Use noon rhomax from rhomax.out %d \n",rho_max );
	   printf("Use noon rhomax from rhomax,out %d \n",rho_max );
	   fclose(frhomax); 
       }
       else
       {
	   fprintf(fstatrho,"# WARNING Cannot open rhomax file, hence instead of the noon rhomax that of the actual slot is used,which might lead to higher uncertainty.\n");
	   rho_max=(int)fper[pm]*esd*1.1;
	   fprintf(fstatrho,"# Use rhomax of actual slot %d \n",rho_max );
	   /*fclose(frhomax);*/
       }  
   }
   sig_ground=(short int)(0.035 * rho_max);
   fprintf(fstatrho,"#  sig_ground =%d \n \n",sig_ground );
   fclose(fstatrho);
   free(fper);
   /* <--- SELFCALIBRATION: Open and write the output files for information  related to  the 
      self calibration approach */  
   printf("Finished normalizing images and automatic calibration \n");
 /*
  *  Grundalbedo and Cloud Index  ----------------------------------------
  */

 for (lin = 0;lin < geoinfo.nrlines;lin++)   /* fuer alle Pixel (Zeilen ... */
 {
   for(col = 0; col < geoinfo.nrcolumns; col++)  /* ... und Spalten)            */
   { 
     for(k = 0; k<=n; k++)
     {
       /* copy all values to a vector, from these values ground albedo will de determined 
          reflectivity_vector contains all reflectivity values for one location */
        reflectivity_vector[k]=bildfolge[k].imagedata[lin][col]/bildfolge[k].f;
	/* /f correct for earth sun-distance **/
     } 
	 /* calculation of ground albedo: 
	    out[lin][col]contains ground albedo
	    reflectivity_vector contains cloud index values */
     out[lin][col]=groundreflectivity_and_cloudindex(&reflectivity_vector[0], n, sig_ground, sig_shadow, rho_max);
     for(k = 0; k<=n; k++)
     { 
       /* copy all cloud index values from vector back to matrix */
       bildfolge[k].imagedata[lin][col]=reflectivity_vector[k];

     }   
   }
  }
 
/***********************************************************************
  *  Ausgeben aller Bilder und Freisetzung des Speicherplatzes...
 ***********************************************************************/ 
  
  /* Bodenalbedonamen zusammensetzen  */
  /* RM changed the naming of the reflectance files, mothly means for every slot are calculated and has to be saved,
     with the old naming they are  overriden */ 
     strcpy(out_imagefile,groundpathname);
     strcat(out_imagefile,"/");
     sprintf(zeitname,"%04d",geoinfo.aq_time);
     strncat(out_imagefile,bildfolge[1].name,6); /* added string for year and month */
     /*strncat(out_imagefile,zeitname,4);*/
     strcat(out_imagefile,zeitname);
     /* printf("out_imagefile = %s bildfolge.name= %s\n",out_imagefile, bildfolge[1].name); */ 
     if(HDL==256)strcat(out_imagefile,"hhmm.REF"); 
     if(HDL==260)strcat(out_imagefile,"hhmm.REF.X260"); 
  /* Grundalbedo speichern */
     
   printf("nrbytes= %d",geoinfo.nrbytes) ;
  if (!(out_image=fopen(out_imagefile,"wb"))) 
     exit_no_write(out_imagefile);
    
  if (MFG==1) 
  {geoinfo.nrbytes=geoinfo.nrbytes*2;}   /* dangerous, only for testing */

   fwrite(&bildfolge[1].headline[0],HDL,1,out_image);    

   /* orig fwrite(&out[0][0],geoinfo.nrcolbig*geoinfo.nrlines*(long)geoinfo.nrbytes,1,out_image);  
   fclose(out_image);
   free_smatrix(out,0,geoinfo.nrlines-1,0,geoinfo.nrcolbig-1); */
   fwrite(&out[0][0],geoinfo.nrcolumns*geoinfo.nrlines*(long)geoinfo.nrbytes,1,out_image);  
   fclose(out_image);
   free_smatrix(out,0,geoinfo.nrlines-1,0,geoinfo.nrcolumns-1);
  
  for(k=0;k<=n;k++)             
  {  
     strcpy(out_imagefile,cloudpathname);
     strcat(out_imagefile,"/");
     strncat(out_imagefile,bildfolge[k].name,strlen(bildfolge[k].name)-4);
     if(HDL==256)strcat(out_imagefile,"CI.XPIF");
     if(HDL==260)strcat(out_imagefile,"CI.X260"); 
     if (!(out_image=fopen(out_imagefile,"wb"))) 
         exit_no_write(out_imagefile);
     fwrite(&bildfolge[k].headline[0],HDL,1,out_image);   
     /*  orig
	 fwrite(&bildfolge[k].imagedata[0][0],geoinfo.nrcolbig*geoinfo.nrlines*(long)geoinfo.nrbytes,1,out_image);  */
   
     /* memcpy(void *dest, const void *src, size_t n);*/
     bildfolge[k].Bimagedata=bmatrix(0,geoinfo.nrlines-1,0,geoinfo.nrcolumns-1);
     for (lin = 0;lin < geoinfo.nrlines;lin++) {  
	 for(col = 0;col < geoinfo.nrcolumns;col++){ 
	     bildfolge[k].Bimagedata[lin][col]=(byte)(bildfolge[k].imagedata[lin][col]/4);
	 }
     }
	 fwrite(&bildfolge[k].Bimagedata[0][0],geoinfo.nrcolumns*geoinfo.nrlines,1,out_image);
	 free_bmatrix(bildfolge[k].Bimagedata,0,geoinfo.nrlines-1,0,geoinfo.nrcolumns-1);
	 fclose(out_image);
	 free_smatrix(bildfolge[k].imagedata,0,geoinfo.nrlines-1,0,geoinfo.nrcolumns-1); /* orig geoinfo.nrcolbig-1 */
  }

  return(0);
}


/***********************************************************************
 *  exit_*
 *
 *    different error messages
 *
 ***********************************************************************/
void exit_print_usage(void)
{
  fprintf(stderr,"\n  Iterative calculation of ground reflectivity and cloud index images \n");
  fprintf(stderr,"\n  with cloud shadow algorithm \n");
  fprintf(stderr,"\n  usage: albedo_and_cloudindex_shadow_msg -l <listname> [other options]\n");
  fprintf(stderr,"\n  options: -l <listname>       : list of image names\n");
  fprintf(stderr,"\n           -b <imagepathname>  : where do I find the images?\n");
  fprintf(stderr,"\n           -g <groundpathname> : where to put ground reflectivity image?\n");
  fprintf(stderr,"\n           -c <cloudpathname>  : where to put cloud index images?\n");
  /* RM 072008 replaced sig_ground bny dpos */
  fprintf(stderr,"\n           -s <dpos>   : default 0, position of Met relativ to lon=0 in degree, West is negative \n");
  fprintf(stderr,"\n            -z <sigma_shadow>   : default 5, no shadow 0?\n");
  

  exit(1);
}

void exit_no_list(char *name)
{
  fprintf(stderr,"groundalbedo_and_cloudindex: No valid list with images in %s\n",name);
  exit(1);
}

void exit_no_write(char *name)
{
  fprintf(stderr,"groundalbedo_and_cloudindex: No write access to %s\n",name);
  exit(1);
}

void exit_no_read(char *name)
{
  fprintf(stderr,"groundalbedo_and_cloudindex: No read access to %s\n",name);
  exit(1);
}


/* -------------------------------------------------------------------------------------------------- */
short int groundreflectivity_and_cloudindex(short int *reflectivities, long zahl, short int sigma_g, short int sigma_s, short int rhomax)
{
  /* calculates ground reflectivity and cloud index

     the resulting ground albedo contains no shadows
     shaded pixels get an appropriate cloud index value for overcast situations

     input:  reflectivities organized in a vector (20-30 values)
     output: one ground reflectivity value
             cloud index values organized in a vector


     Dr. Annette Hammer, University of Oldenburg, July 2005
     based on algorithm of Rolf Kuhlemann, University of Oldenburg
  */
  
  short int k,kk=0; /* added kk in order to modify while loop */
  short int cloudy;                /* counts cloudy values   */
  short int cloud_free_quantity;   /* counts cloudfree values */
  double average;                  /* average of cloudfree values == ground reflectivity*/
  short int shadow_count,shadow_det,
        dist,p,aveold,low,distold;
  short int reflec[BILDANZ], q[BILDANZ];


  cloud_free_quantity = shadow_det = zahl;

  shadow_count = 0 ;
  for(k = 0; k<=zahl; k++)q[k]=0;
  
  do                       /* Iteratives Ausfiltern von Schattenpixeln */
    {
      for(k = 0; k<=zahl; k++)
	{
	  reflec[k]=reflectivities[k];
	}
      shadow_det=dist=p=0;
      average=0.0;
      aveold=BYTE_MAX;         /* Nur fuer Detektionsschleife hochgesetzt */
      do                          /*  Iteratives Ausfiltern bewoelkter Pixel */
	{
	  average = 0.0;
	  low =BYTE_MAX ;      /* Nur fuer Detektionsschleife hochgesetzt */
	  for(k = 0; k<=zahl; k++)
	    {
		if( k != (q[k]-1)) 
              /* q[k]-> q[k]-1 RM changed in order to be consistent with loop from 0 instead of 1 */
	      {
		average+=reflec[k];
		  {
		    if(low > reflec[k] && reflec[k] != 0)
		      {
			low = reflec[k];
			p = k;   /* Tag des dunkelsten Pixels eines Slots */
		      }
	           }
               }
	    }
         average = average/(cloud_free_quantity - shadow_count);   /* berechne den Mittelwert */
	 cloudy = cloud_free_quantity = 0;
	 for(k = 0; k<=zahl; k++)                   /* Pixel bewoelkt ? */
	    {
	      if(reflec[k] > (average + sigma_g))
		{
		 reflec[k]  = 0;   /* bewoelkte Pixel bei Mittelung nicht beruecksichtigten */
		 cloudy++;
		}
	      else
	        {
		  if(reflec[k]) ++cloud_free_quantity;
		  /* unbewoelkt und kein Sensorausfall,
				     und nicht vorher ausgefiltert */
		}
	    }

	  if (sigma_s>0)
	  { 
	   distold=dist;
	   dist=aveold - (int) average;    /* Schrittweitenvergleich + Sigma-Schatten */
	   if(dist > 0 && distold > 0 && dist > (distold + sigma_s) && average)
	    {
	      shadow_count++;                /* Zaehler, um jeweils fuer richtige Detektion,   */
	      shadow_det++;                  /* richtigen Schleifenabbruch und richtige        */
	      q[p] = p-1;                      /* Statistik zu sorgen.                           */
	      /* p-> p-1 RM changed to be consistent with loop from 0 instead of 1 */
		                             /* q[k] belegt dunkelstes Pixel und ist */
              cloud_free_quantity=zahl;      /* noetig, um gleiche Pixel wie Schatten mit      */
	      cloudy = 0;                    /* korrigiertem Cloudindex zu berechnen.          */
            }
	   aveold= (short int) average;
	  }
          kk=kk+1; /* RM 02,2009 have been added in order to be able to modify while lop */
        }while(kk <= zahl); /* the orig while condition leads to endless loops on some machines !! */
      /* orig  while(cloudy  && cloud_free_quantity); */ /* solange bewoelkte Pixel aus gefiltert werden
	                                          und Mittelwertbildung */
      
    } while(shadow_det != 0 );
    
    for(k = 0; k<=zahl; k++)                   /* cloudindex bestimmen */
       {
	reflectivities[k]=cloudindex(reflectivities[k],(short int)average,rhomax,q[k]);
       }
   return((short int) average);       /* verbleibender Mittelwert wird als Grundalbedo zurückgegeben */
}

/*------------------------------------------------------------------*/


short int cloudindex(short int reflectivity,short int ground,short int rhomax,short int shadow)
{
 double n,relation_k;
 
 if(reflectivity==0)
   {
    return(NO_VALUE);
   }
 else
 {
 if(shadow==0) /* Heliosat method as usual */
   {
   if(rhomax!=ground) 
     {n=(reflectivity-ground+0.0)/(rhomax-ground);}
   else 
     {n=CLOUDI_MAX;}
   }
 else         /* if there was a cloud shadow */
   {
    if(ground>0) 
       {
        relation_k=reflectivity/ground;  /* ist ein clearsky index ! */
	n= 1-relation_k;
	if (n > 0.80)
	    {
	     n = (3.6667 - sqrt(13.4447 + 6.6668 * (relation_k - 2.0667))) / 3.3334 ; 
	     /* cloud index aus clearsky index */
	     /* Umkehrabbildung zu clearsky index aus cloud index*/
	    }
       }
    else {n=CLOUDI_MAX;}
   }
   
   n=maxi(n,CLOUDI_MIN);
   n=mini(n,CLOUDI_MAX);
      
   return ((short int)(slope*n+y_intercept));
  }
}


/*------------------------------------------------------------------*/

#define NR_END 1
#define FREE_ARG char*

short int **smatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a short int matrix with subscript range m[nrl..nrh][ncl..nch] 
   A. Hammer 2005 - concept of Numerical Recipes in C */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        short int **m;
        
        /* allocate pointers to rows */
        m=(short int **) malloc((size_t)((nrow+NR_END)*sizeof(short int*)));
        if (!m) {
	         fprintf(stderr,"allocation failure 1 in smatrix()");
		 exit(1);
		 }
        m += NR_END;
        m -= nrl;


        /* allocate rows and set pointers to them */
        m[nrl]=(short int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(short int)));
        if (!m[nrl]) {
	               fprintf(stderr,"allocation failure 2 in smatrix()");
		        exit(1);
		      }
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

/*--------------------------------------------------------------------------*/ 

void free_smatrix(short int **m, long nrl, long nrh, long ncl, long nch)
/* free a short int matrix allocated by smatrix()
   A. Hammer 2005 - concept of Numerical Recipes in C  */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

/*--------------------------------------------------------------------------*/ 

short int **read_xpif_to_matrix(short int xpif_head_length,char *image_name,struct bildinfo *imageinfo)
{ 
unsigned char headline[HDL];        /* die Headerzeile */ 
 short int **imagedata;
 short int nrcolcorr;
FILE *in_image;                         /* FILEname des Eingabebildes */

 if (!(in_image=fopen(image_name,"rb")))   /* das Bild muss zum Lesen zu oeffnen sein            */
     exit_no_read(image_name); 
 /*
   *  Auswertung des Headers und Einlesen des Bildes in eine Matrix..............................
  */

 if(xpif_head_length==256)
 {  
  fread(&headline[0],xpif_head_length,1,in_image);
  xpifinfo(&headline[0],imageinfo); 
 }
 else if(xpif_head_length==260)
 {
   fread(imageinfo,xpif_head_length,1,in_image);
  (*imageinfo).nrcolbig=(*imageinfo).nrcolumns;

 }
 else
  {
  fprintf(stderr,"Crash: The given headerlength is not 256 or 260!\n");
  exit(1);
  }  
  fseek(in_image, (long) (((*imageinfo).add_header)*RECLENG+xpif_head_length), SEEK_SET);
  imagedata=smatrix(0,(*imageinfo).nrlines-1,0,(*imageinfo).nrcolbig-1);
  fread(&imagedata[0][0],(*imageinfo).nrcolbig*(*imageinfo).nrlines*sizeof(short),1,in_image); 
  fclose (in_image);
  return imagedata;

}
/*--------------------------------------------------------------------------*/ 

void read_xpif_header(short int xpif_head_length,char *image_name,unsigned char *headline)
{ 
FILE *in_image;                         /* FILEname des Eingabebildes */

if (!(in_image=fopen(image_name,"rb")))   /* das Bild muss zum Lesen zu oeffnen sein            */
     exit_no_read(image_name); 
 /*
   *  Einlesen des Headers .............................
  */

 if(xpif_head_length==256 || xpif_head_length==260)
 {
  fread(&headline[0],xpif_head_length,1,in_image);
 }
  else
  {
  fprintf(stderr,"Crash: The given headerlength is not 256 or 260!\n");
  exit(1);
  } 
 fclose(in_image);
}

/* RM FEB 2009, added stuff to read in the MFg data ,bmatrix and read_xpif_Bmatrix */
/* RM FEB 2009 -----> */

byte **bmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a byte matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    byte **m;

    /* allocate pointers to rows */
    m=(byte **) malloc((size_t) ((nrow+NR_END)*sizeof(byte*)));
    if (!m) {	         
	fprintf(stderr,"allocation failure 1 in smatrix()");
	exit(1);
    }
    m += NR_END;
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(byte *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(byte)));
    if (!m[nrl]) 
    {fprintf(stderr,"allocation failure 2 in bmatrix()");
	exit(1);
    }
    m[nrl] += NR_END;
    m[nrl] -= nrl;

    for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

void free_bmatrix(byte **m, long nrl, long nrh, long ncl, long nch)
/* free a byte matrix allocated by bmatrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

unsigned char **read_xpif_to_Bmatrix(short int xpif_head_length,char *image_name,struct bildinfo *imageinfo)
{ 
unsigned char headline[HDL];        /* the header line */ 
byte **imagedata;                  /* alternatively unsigned char */
 short int nrcolcorr;
FILE *in_image;                         /* filename of the file containing the image data */

 if (!(in_image=fopen(image_name,"rb")))   /* check if the file can be open for reading            */
     exit_no_read(image_name); 
 /*
   evaluation of the XPIF header and assignment of the image data to a byte matrix */ 
 

 if(xpif_head_length==256)
 {  
  fread(&headline[0],xpif_head_length,1,in_image);
  xpifinfo(&headline[0],imageinfo); 
 }
 else if(xpif_head_length==260)
 {
   fread(imageinfo,xpif_head_length,1,in_image);
  (*imageinfo).nrcolbig=(*imageinfo).nrcolumns;

 }
 else
  {
  fprintf(stderr,"Crash: The given headerlength is not 256 or 260!\n");
  exit(1);
  }  
    fseek(in_image, (long) (((*imageinfo).add_header)*RECLENG+xpif_head_length), SEEK_SET);
  imagedata=bmatrix(0,(*imageinfo).nrlines-1,0,(*imageinfo).nrcolumns-1);
  fread(&imagedata[0][0],(*imageinfo).nrcolumns*(*imageinfo).nrlines*sizeof(char),1,in_image); /* orig nrcolbig */
  fclose (in_image);
  return imagedata;

}
/* <----- RM FEB 2009 */

/*--------------------------------------------------------------------------*/ 



/* RM May 2008 implemented function to sort an vector into ascending numerical order , 
  from Numerical Recipes in C, Secon Edition, Cambridge University Press */
void heapsort(unsigned long n, short int ra[])
{
	unsigned long i, ir, j, l;
	short int rra;
	
	if(n < 2) return;
	l=(n>>1)+1;
	ir=n;
	
	for(;;)
	{
		if(l>1)
		{
			rra = ra[--l];
		}
		else
		{
			rra = ra[ir];
			ra[ir] = ra[1];
			if(--ir == 1)
			{
				ra[1] = rra;
				break;
			}
		}
	
	i = l;
	j = l+l;
	while(j <= ir)
	{
		if(j < ir && ra[j] < ra[j+1])
		j++;
		
		if(rra < ra[j])
		{
			ra[i] = ra[j];
			i = j;
			j <<= 1;
		}
		else
		j = (ir+1);
	}
	ra[i] = rra;
	}
}

bool file_exists(const char *filename)
{
   FILE *f = fopen(filename, "r");

   if (!f) return false;

   fclose(f);

   return true;
}

/*bool file_exists(const char *filename)
{
    if (FILE *file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    return false;
    }*/

#include "./satgeo-32b.c"
