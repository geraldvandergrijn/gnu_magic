#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include "magic.h"     /* the magic.h declarations */
#include "nr_magic.h" /* definition of byte matrix and vector */
#include "nrutil.h"   /* math definitions from NR */
#include "magic-geo.h" /*the definition of variables needed for geolocation and geometry calculation */
/* includes function regrid2, which transform irregular 
 satellite grid to regular user defined grid */
/* todo: arrays start at index 0, interpolation assumes start at 1 ?
borders and lint are already redifened for arrays 
starting at 0, are there ohters ? */
/* global variables */
/* compile with: gcc-3.4  -lmathlib magic-clear-s1.c -o tmp.exe
or...   gcc-3.3  -lm magic-v0x71.c -o tmp.exe
 */
float sza,coszen; /* solar zenith angle, cosine of sza */
int main(void)
{
  /** general info
 programer: R. Mueller. 
 PURPOSE:
  A: Calculates the clear sky global and direct irradiance, stand alone  
  B: Calculates the all sky global and direct irradiance.
     pre-requisite for option B is the availibility of 
     cloud index/cloud albedo or cloud transmission information.
     Cloud albedo information can also be derived from cloud optical depth
     using the cloud index fomrula.

 NEEDED INPUT: 
   - A aerosol climatology (default KINNE or GADS/OPAC),
   - A water vapour climatology (default NCEP).
   - optional, an ozone climatology, if not available default values
     will be used (is currently not implemented). 
   ALL clear sky climatologies are in ASCII, the used climatologies
   can be replaced by climatologies the user selected. The dimension
   is defined within the config file. However the ASCII format of the
   climatologies has to be identical, see the respective clims for details. 
   - optional a binary file (byte or short int with or without header)
     containing information of cloud
   - parameter/variables for the steering of the program, 
     all provided in the configuration file  magic-config.inp            
 OUTPUT: A file conatining the irradiance data on either 
 a "regular grid" in ASCII format or optinally on Meteosat 
 projection in binary format.
 VERSION: 0.1

 /* to do add funcition to calculate SZA for given date & time*/ 

  /* ACKNOLEDGMENT: Thanks to the libRadtran team. LibRadtran has been used
 to derive the parameterisations used in this code !
**/
  /* the program starts with the declaration of the needed variables 
   and parameters */
  int xadim, yadim, xhdim, yhdim; /*dimension of aerosol and h2o input */ 
  int naod=10, nssa=3, ngg=2, nh2o=18, no3=3; /* lut dimension, do not change */
  int xdim,ydim,region,cloudcalc; /* dimension of cloud file, cloudcalc flag */
  float lonb,lone,latb,late,dlon,dlat; /* longitude, latitude information for region*/

  int xodim, yodim, nsal, xndim, yndim,dummy,xrdim,yrdim;
  int i,j,k,ii,jj, month, year, doy, hour, minute,gads; /* doy=day of year !*/
  int nsza=6, kh, kl,khs,kls,col,lin,coloffset=0,linoffset=0,regrid=0 ;
  int rlin,rcol;
  float w,ws,stime;
  char aeroclim[100], cross, hclim[100],clut[100],salclim[100],luclim[100];
  /* names of input files */  
  char cdummy[100];
  /* char aeroclim[100], header[120], cross, hclim[100],clut[100]; */
  char outfile[100],nfile[100],nlonfile[100],nlatfile[100],nszafile[100],cdlfile[100];   
  /* name of input / output files */
  float  aod, ssa, H2O, O3, rgrid,ozone; /* atmospheric parameter ,rgrid=regrid flag */
  float Gtmp[2],Gmlb,Btmp[2],Bmlb,nmin,nmax; /* the solar irradiance quantities */
    double salb=0.1,xlon,ylat,xplon; /* latitude, longitude of site */
  int xsdim=2160,ysdim=1080; /*hard wired dimensions for SARB/CERES SAL */
  int bsize,klat,klon; /* bsize, size of one element of cloud array */
  float ni,ki,GMT;
  int headleng=0; 
  FILE *outGcld,*outGclear;

  struct SSIC{
   float **g;  /* the clear sky solar surface irradiance */
   float **b; /* the clear sky beam (direct) irradiance */
  }clear;

  struct image{
    FILE *f ;    /* pointer to the image containing the cloud index */
    char header[120];    /* header */
    byte **n; /* the cloud index as byte */ 
    short int **in; /*the cloud index as short int */
    float **g; /* the all sky solar surface irradiance */
    float **b; /* the all sky beam irradiance */
  };
  struct image cloud;   /* cloud[n] read n images*/
  int read_n;
 

  float *lonreg,*latreg;

  /* using double instead of float leads to zero values in the printf 
      command, using %lf instead of %f necessary in this case*/
 

  FILE *configf; 
  configf = fopen("magic-config.inp","r");  
    if (configf == NULL){
      printf("ERROR: config file seems to be not available \n");
    }
  printf("start reading config file");
  fscanf (configf, "%s \n", &cdummy);
  fscanf(configf, "%d %d %d %d %d", &year, &month, &doy, &hour, &minute);
  /* read the dimensions and the name of the aorosol climatology */
  fscanf (configf, "%s \n", &cdummy);
  fscanf(configf, "%d %d %s", &xadim, &yadim, &aeroclim);
  fscanf (configf, "%s \n", &cdummy);
 /* read the diemensions and the name of the water vapour climatology */
  fscanf(configf, "%d %d", &xhdim, &yhdim);
  fscanf(configf, "%s %d", &hclim);
  fscanf (configf, "%s \n", &cdummy);
  /* read the name of the clear sky look-up table */
  fscanf(configf, "%s ", &clut);
  /* read the name of the surface albedo and the landuse climatology */
  fscanf(configf, "%s %s", &salclim, &luclim);
  fscanf(configf, "%f", &ozone ); 
 /* read the name of the output file for binary and cdl ASCII output */
  fscanf (configf, "%s \n", &cdummy);
  fscanf(configf, "%s", &outfile);
  fscanf(configf, "%s", &cdlfile);
  /* read the dimension of the all sky G and beam output */
  fscanf (configf, "%s \n", &cdummy);
  fscanf(configf, "%d", &cloudcalc);
  fscanf(configf, "%d %d", &bsize, &headleng); /* size in byte and headleng */
  fscanf(configf, "%d %d", &xndim, &yndim);   /* dimension of cloud index */
  fscanf(configf, "%s", &nfile); /* the cloud index file */
  fscanf(configf, "%f %f %f", &satpos,  &aperture, &stime); /* satellite position & aperture angle,
  and time satellite needs to scan 1 line */ 
  fscanf(configf, "%f %f", &nmin,  &nmax); /* dimension of the cloud index file */
  fscanf(configf, "%d", &regrid); /* regrid in case of clouds ?*/
  fscanf (configf, "%s \n", &cdummy);
  fscanf(configf, "%f %f %f %f %f %f", &lonb,&lone,&latb,&late,&dlon,&dlat);
  /* define the region */ 
  fscanf(configf, "%s", &nlonfile);
  fscanf(configf, "%s", &nlatfile); /* the associated lat,lon files */
  fscanf(configf, "%s", &nszafile); /* the associated sza file */
  fclose(configf);
  printf("reading config file finished");  
  

  FILE *logfile;
  logfile = fopen("logfile.out","w");

/* fgets(hclim,sizeof(hclim),configf); */
/* read the dimension of the aerosol climatology */
  fprintf(logfile, "year,month,day of year,hour,minute: \n %d %d %d %d %d \n", year,month,doy,hour,minute);
  fprintf(logfile, "the dimension of aerosols and h2o climatology \n");
  fprintf(logfile, "%d %d %d %d \n", xadim, yadim, xhdim, yhdim);
  fprintf(logfile, "the name of the aerosol, h2o  and surface albedo and landused climatology \n");
  fprintf(logfile,"%s\n %s\n %s\n %s\n", aeroclim,hclim,salclim,luclim);
  fprintf(logfile, "the dimension of the LUT parameter, ngg,nssa,naod,nh2o,no3: \n");
  fprintf(logfile, "%d %d %d \n %d %d \n", ngg,nssa,naod,nh2o,no3);   
  if (ngg > dgg || nssa >dssa || naod > daod || nh2o > dh2o || no3 > do3)
    {fprintf(logfile,"one of the dimension parameter is larger than the arrays define in the LUT strucuture, see magic.h, this leads to segmentation fault ! \n");
    fprintf(logfile,"the correct values in the config files are:\n 2 3 10 \n 18 3");}
  if (ngg != dgg || nssa != dssa != naod > daod != nh2o > dh2o != no3 > do3)
    {fprintf(logfile,"The dimension parameter does not match the definition of the original LUT,please be aware what you are doing, see definition of struct in magic.h \n");
    fprintf(logfile,"the correct values in the config files are:\n 2 3 10 \n 18 3");
    }
  fprintf(logfile, "the name of the look-up table: %s \n", clut);
  /* printf("%s\n", aeroclim);
     printf("%s\n", hclim); */
  fprintf(logfile, "the names of the output files: \n %s %s \n", outfile, cdlfile);   
  fprintf(logfile, "calculation for cloudy sky 0=false, 1=true: %d \n", cloudcalc );
  fprintf(logfile, "if cloud true:output for regrided area  0=false, 1=true: %d \n", regrid );
  fprintf(logfile, "The defined sub-region, longitude:start/end latitude: start/end, x,y-resolution \n");
  fprintf(logfile, "%f %f %f %f %f %f \n", lonb,lone,latb,late,dlon,dlat);
  fprintf(logfile, "the following parameter are only needed if the cloud option is used \n");
  fprintf(logfile, "the file containing the cloud information: \n %s \n", nfile);
  fprintf(logfile, "the dimension of the cloud file xdim,ydim= %d %d \n", xndim,yndim);
    fprintf(logfile, "The headerlength is %d and the size of the array elements is %d \n",headleng, bsize); 
  fprintf(logfile,"aperture = %f satpos= %f scanning time for 1 line %f\n", aperture,satpos,stime);
  fprintf(logfile,"the heliosat cloud index ranges from -0.2 to 1.2 \n");
  fprintf(logfile,"the used range is nmin, nmax %f %f", nmin, nmax);
  fprintf(logfile, "currently not in use \n %s\n %s\n %s\n", nlonfile,nlatfile,nszafile);
  fclose(logfile);
   
  /*  struct geometry
  { float nlon[xndim][yndim];
    float nlat[xndim][yndim];
    float nsza[xndim][yndim];
  } geo;
  leads to core dump  */
  
  rhclim(xhdim,yhdim,month,hclim); /* function reads the h2o information */ 

  raeroclim(xadim,yadim,month, aeroclim); /* function reads the aerosol information */

  rclut(ngg,nssa,naod,nh2o,no3,clut); /* function reads the LUT */

  rsarbsal(salclim,luclim,xsdim,ysdim);  /* function reads the surface albedo */

  /* read_something(int xdim, int ydim, int month, char clim[], int nr);*/

   /* read the cloud index, !! image[lin][col] */
  if (cloudcalc != 0){ 
    printf("read the n file, xndim= %d yndim= %d \n", xndim,yndim);
    cloud.f = fopen(nfile,"rb");  
    if (cloud.f == NULL){
      printf("ERROR: cloud index file seems to be not available %s \n",cloud.f);
    }
    else{
      cloud.in=smatrix(0,yndim-1,0,xndim-1);
      switch (bsize)
	{
	case 1: {
	  cloud.n=bmatrix(0,xndim-1,0,yndim-1);
	  fseek(cloud.f, (long) headleng, SEEK_SET);
	  read_n=fread(&cloud.n[0][0],bsize,xndim*yndim,cloud.f);
	  /*	  &cloud.in[0][0]=4*(short int)&cloud.n[0][0]; might work ?*/
	  for (lin=0;lin<yndim;lin++){
	    for (col=0;col<xndim;col++){
	      cloud.in[lin][col]=4*(short int)&cloud.n[lin][col];
	    }
	  }
	  if (read_n*bsize != xndim*yndim*bsize){ 
	    printf("problem reading n-image, xdim ydim does not match size of image %d %d \n",                       read_n*bsize, xndim*yndim*bsize);
	    
	    /*  int a[10][5];
		int b[10][5];
		  memcpy(b[0], a[0], 10 * 5 * sizeof(int));*/
	    
	  }
	  free_bmatrix(cloud.n,0,yndim-1,0,xndim-1); 
	  break;
	}
	case 2: { /* lines to read SHORT INT matrix -->*/
	  /* read_n = fread(&cloud.n[0][0],xndim*yndim,1,cloud.f);*/
	  /* <---*/
	  fseek(cloud.f, (long) headleng, SEEK_SET);
	  read_n=fread(&cloud.in[0][0],bsize,xndim*yndim,cloud.f);
	  if (read_n*bsize != xndim*yndim*bsize){   
	    printf("problems reading image, xdim ydim does not match size of image %d %d \n",                        read_n*bsize, xndim*yndim*bsize);
	  } 
	  break;
	}
	default: {
	  printf("type (binary size) of cloud information currently not supported,                                  only 1byte and 2 byte types (e.g. short int) are supported");
	    }
	  /*exit(8); */
	} 	  

      xdim=xndim;
      ydim=yndim; 
      cloud.g=matrix(0,ydim-1,0,xdim-1); /* allocate matrix for all sky G */
      cloud.b=matrix(0,ydim-1,0,xdim-1); /* allocate matrix for all sky beam */
      /* allocate the matrix for the clear sky solar surface irradiance and beam irradiance */
      clear.g=matrix(0,ydim-1,0,xdim-1);
      clear.b=matrix(0,ydim-1,0,xdim-1);
    } 
    fclose(cloud.f);
  } 
  else{
    /* allocate the matrix for the clear sky solar surface irradiance and beam irradiance,
       using defined region in the input file to calculate xdim,ydim */
    xdim=(int)(fabs((lonb-lone)/dlon)+1); 
    ydim=(int)(fabs((latb-late)/dlat)+1);
    printf("only clear sky for defined region, xdim=%d, ydim=%d \n",xdim,ydim);
    clear.g=matrix(0,ydim-1,0,xdim-1);
    clear.b=matrix(0,ydim-1,0,xdim-1);
  }
  
  
  /* EVERYTHING SHOULD BE PREPARED  START WITH THE CALCULATION NOW!!! */
  
  
  for (lin=0;lin<ydim;lin++){  /* lin=line (lat,yvalues)*/
    for (col=0;col<xdim;col++){  /* col=coloumn (lon,x-values) */
      /* 10.0+l*0.05; */
       if (cloudcalc != 0){
         rlin=linoffset+lin;  
  /* it is assumed that lin, col starts at ?? upper left corner 
     in the cloud image file */
         rcol=coloffset+col;
	 geoloc(rlin,rcol,yndim,&ylat,&xlon);
         /* printf("i,j,lon,lat= %d %d %f %f \n",l,m,xlon,ylat); */
	 GMT=hour+(double)minute/60.+(((float)rlin/(float)xdim)*((float)stime/60.)); 
	 /* if (GMT > hour+30 || GMT < hour-30){printf("GMT= %f \n");}*/
	 /* +(12./xdim/60.) has to be changed add sanning of sat */
	 /* minute in 60 or 100 units expected ?? */
       }
       else{
	 xlon=lonb+col*dlon;
	 ylat=latb+lin*dlat;
	 GMT=(float)hour+(double)minute/60.0;  
	 /* +dlon has to be changed add scanning of sat */
	 /* minute in 60 or 100 units expected ?? */
       }
       /*  if (cloud !=0 && (xndim =< xdim || yndim =< ydim ))
	   {printf("Error if Magic is used for cloudy sky xndim,yndim has to be greater than xdim,ydim");} */
       if(ylat < -89.5 || ylat > 89.5 || xlon == 999){
	 clear.g[lin][col]=-1;
       }
       else{ 
	 double gamma,f,dec,lon_sun,EOT;
	 EOT=equation_of_time(doy,&gamma,&f,&dec);
	 lon_sun=(12.0-GMT-(EOT/60.0))*15*(PI/180.0);
	 coszen=cos_zenit(ylat*deg2rad,xlon*deg2rad,dec,lon_sun);
         sza=acos(coszen); 
	 /* printf("geo= %f %f %f %f %f %f \n",ylat,xlon, EOT, GMT, lon_sun, sza) ;*/
	 /* GMT=time (hout+minute) 
	    bildfolge[k].lon_sun=(12.0-GMT-(bildfolge[k].EOT/60.0))*15*(PI/180.0); */
	 
	 /*  printf("ylat,xlon= %f %f \n", ylat,xlon); */
	 /* aerosol and  h2o climatology expected 
	    to be in equidistant lat, lon grid, max resolution currently 1 deg*/
	 aod=interpolate2(xlon, ylat,xadim,yadim,aero.aod,aero.lon,aero.lat);
	 ssa=interpolate2(xlon, ylat,xadim,yadim,aero.ssa,aero.lon,aero.lat); 
	 H2O=interpolate2(xlon, ylat,xhdim,yhdim,h2o.val,h2o.lon,h2o.lat);
	 /*printf("lin,col,H2O= %d %d %f \n", lin,col,H2O);*/
	 
     	 if (xlon < 0  && sal.lon[xsdim-1] > sal.lon[0]){  
	   /* currently needed, magic assumes input to range from
	      -180 - 180, the SARB SAL goes from 0-360 in 1/6 degree steps */
	   xplon=360.0+xlon;
	   /* printf("xplon2 = %f \n",xplon) ; */
	   if (xplon >= (360.0 - 0.08333)){xplon=0;} 
	   if (xplon > sal.lon[xsdim-1] && xplon <  360.0-0.08333)
           {xplon=sal.lon[xsdim-1];
	   /* printf("xplon3 = %f %f \n",xplon,sal.lon[xsdim-1]) ;*/
           } 
	   /* printf("xplon3 = %f %f \n",xplon,sal.lon[xsdim]) ;*/
         }
	 else 
	   {xplon=xlon;
	   sal.lon[0]=0.0;} 
     /* only for sal from SARB, wich goes from 0-360 degree instead
          of the standard 
     for security implement such checks also for the others & replace
       */ 

     salb=interpolate2(xplon,ylat,xsdim,ysdim,sal.val,sal.lon,sal.lat);  
     /* the albedo map shows snow south-west of the Alps, close to the meditereanen sea
        which is unlikely to be real, correct for it */
     if(salb > 0.3 && xplon > 6.0 && xplon < 18.0 && ylat > 41.0 && ylat < 45.0 )
       {salb=0.2;} 
     /* get the interpolated SAL value for a SZA of 60 degree */
     klat=borders(sal.lat,ysdim,ylat);
     klon=borders(sal.lon,xsdim,xplon);
  /* get the landuse of the nearest neighbour */
     /* PT  printf ("salb for correction = %f", salb); */
     if(salb<0){printf("salb before cor %f %f %f %f \n", xplon,xlon,ylat,salb);} 
     salb=salb*corsal(sal.lu[klon][klat],sza);
     if (salb >= 0.9){salb=0.9;} /* albedo map and landuse map are not 100% consistent
	    leading to values above 1 for regions with snow in Albedo but tundra in landuse */
     if(salb<0){printf("salb after cor %f %f %f \n", xplon,ylat,salb);}
  /* correct the surface albedo to the actual SZA */
     /* PT printf("interpol lon= %f lat= %f aod %f ssa %f h20 %f sal %f\n", xlon, ylat, aod, ssa, H2O,salb); */
    /*printf("H2O=  %f\n",H2O) ; */
    /* now start to caclulate the solar irradiance */
    /* 1st step, interpolate G for the aerosol state of the given
       lat,lon */
    /*  1a  find the borders of aod, ssa in the mlb lut x-value vector*/
    /* for testing */
    
    kh=borders(alut.aod,naod,aod);   
    kl=kh-1;
    w=(alut.aod[kh]-aod)/(alut.aod[kh]-alut.aod[kl]);
    /* w is the wheighting factor for var[kl] and 1-w for var[kh] !! */
    if (w < 0 || w > 1){
      printf("problem in function lint w aod is not between 1 an 0 !!!");
      printf("kh= %d %d w= %f\n",kh,kl,w);
      printf("aodlut= %f %f aod= %f\n",alut.aod[kh],alut.aod[kl],aod);
      }
    /* PT printf("kh= %d alut_kh= %f alut_kl= %f w= %f \n",kh,alut.aod[kh],alut.aod[kl],w); */

    if(ssa < 0.7){ssa=0.7;}
    khs=borders(alut.ssa,nssa,ssa);   
    kls=khs-1;
    /* printf("khs= %d %d \n",khs,kls); */
    ws=(alut.ssa[khs]-ssa)/(alut.ssa[khs]-alut.ssa[kls]); 
    /* ws is the wheighting factor for var[kls] and 1-ws for var[khs] !! */
    if (ws < 0 || ws > 1){
      printf("problem in borders w ssa is not between 1 an 0 !!!");
       printf("kh= %d %d w= %f\n",kh,kl,ws);
      printf("ssalut= %f %f ssa= %f\n",alut.ssa[kh],alut.ssa[kl],ssa);
    } 


    /* afterwards interpolate to G using MLB function for the two 
       border aod values, currently for fixed gg  */
     for (i=kl;i<=kh;i++)
       { /*  Gtmp[i]=lint(alut.ssa,nssa,ssa,i);*/
	 /* printf("the pow %f", pow(3,2));*/ /* pow(2,4) = 2**4 */
     Gtmp[i-kl]=ws*alut.Im[1][kls][i]*exp(alut.gtau[1][kls][i]/
     (pow(cos(sza),alut.ag[1][kls][i])))*cos(sza)+(1-ws)*alut.Im[1][khs][i]
      *exp(alut.gtau[1][khs][i]/(pow(cos(sza),alut.ag[1][khs][i])))*cos(sza);  
     /* the direct irradiane = beam */
     Btmp[i-kl]=ws*alut.Im[1][kls][i]*exp(alut.btau[1][kls][i]/
     (pow(cos(sza),alut.ab[1][kls][i])))*cos(sza)+(1-ws)*alut.Im[1][khs][i]
      *exp(alut.btau[1][khs][i]/(pow(cos(sza),alut.ab[1][khs][i])))*cos(sza); 

     /* printf("i-kl=%d \n", i-kl) ; */
   /* printf("1st terms %f \n",alut.ag[1][kls][i]);
     printf("1st term  %d %d %d %f \n",1,kls,i,pow(cos(sza),alut.ag[1][kls][i]));
     printf("the Gi is %f \n", Gtmp[kh-i+1] );  */
     }
     for (i=0;i<=1;i++)
       {
	 Gmlb=w*Gtmp[0]+(1-w)*Gtmp[1];
	 Bmlb=w*Btmp[0]+(1-w)*Btmp[1];
       }
     /* PT printf("that is GMLB= %f w=%f 
	Gtmp1 & 2 %f %f\n",Gmlb,w,Gtmp[0],Gtmp[1]); */
     /* Gmlb=lint(alut.aod,naod,aod); */
     /*  Correct the global irradiance if values are 
       not identical to assumed standart atmosphere,
       15mm H20, SALB=0.2, O3=345 DU */

    /* first the O3 and H20 corretion */
    /* O3 not implemented so far, however works like H20 */
     /* if(Gmlb < 0){printf("Gmlb<0 %d %d %f \n", l,m,Gmlb);} ok*/
     /* multiplation with f takes care for the sun earth distance !!*/
	Gmlb=Gmlb*f+lint(hlut.x,hlut.yg,nh2o,H2O)*pow(cos(sza),0.9);
        Bmlb=Bmlb*f+lint(hlut.x,hlut.yb,nh2o,H2O)*pow(cos(sza),0.9);
	/* printf("G after h2o cor %f \n",Gmlb);  */ 

    /* finally, the albedo correction, the RTM relationship 
       is not linear but can be approximated as linear */
     Gmlb=Gmlb*(0.98+0.1*salb);
     /*if(Gmlb !> -0.01 && !<1368){Gmlb=-9.9;} */
     if(isnan((double)Gmlb)){Gmlb=0.0;}  /* to avoid nans */
     if(isnan((double)Bmlb)){Bmlb=0.0;}

     clear.g[lin][col]=Gmlb; 
     /*  if(clear.g[l][m] < 0){printf("<0 %d %d %f %f\n", l,m,clear.g[l][m],salb);} */
     if(coszen>0.07)
     	clear.b[lin][col]=Bmlb/coszen;
     else
     	clear.b[lin][col]=Bmlb/0.07;

     /* clear.b[lin][col]=ssa; */

     /* the surface albedo does not effect the beam, 
	hence no correction for Bmlb !! */
     /*     fprintf(outf,"final G %f %f %f\n",xlon,ylat,Gmlb); 
       if (fclose(outf) != NULL){
       printf("ERROR: output file can not be closed \n");
  } */
     
     if(cloudcalc != 0 && clear.g[lin][col] != -1)
       {  int drange;
	 /* printf("hello here comes the ne stuff"); */
         switch (bsize){
         case 1: drange=256;
	   break;
	 case 2: drange=1024;
           break;
	 default: printf("not supported")  ;
         }
	 ni=nmin+(nmax-nmin)*cloud.in[lin][col]/drange;
         /* if e.g. COD is scaled from 0-1 using the 95 percentile of COD, 
             this equation has to be used */
         
	 /*245 todo have to be changed to a variable*/
	 /* printf("ni= %f \n",ni); */
	 if(ni<=-0.2)         ki= 1.2; /* todo if loop -0.2<ni to ni>-0.2 */
	 if(-0.2<ni && ni<=0.8)ki= 1-ni;
	 if(0.8<ni && ni<=1.1) ki= 2.0667-3.6667*ni+1.6667*ni*ni; /* overcast */
	 if(1.1<ni)           ki= 0.05;
	 cloud.g[lin][col]=ki*clear.g[lin][col];    
         if(ni<=0.6){ 
           if(ki > 1){ki=1;}  /* we assume no increase of beam due to 3d cloud effect*/
	   cloud.b[lin][col]=clear.b[lin][col]*pow((ki-0.38*(1-ki)),2.5);
         }
         else{
	   cloud.b[lin][col]=0.0;
         }
           /* *dsun*dirfrac*(kindex-0.38*(1-kindex))**2.5*/
	   /* todo dsun has to be defined, dirfrac not needed  */
	 /* printf("the final cloud results %d %d %f %f %f \n",l,m,clear.g[l][m],cloud.g[l][m],ni);*/
      }
     
       }
     }
   }
   
/* free the fields that are no longer needed */
   if(cloudcalc != 0){
   free_smatrix(cloud.in,0,yndim-1,0,xndim-1);
   }
   free_matrix(h2o.val,0,xhdim-1,0,yhdim-1);
   free_vector(h2o.lat,0,yhdim-1);
   free_vector(h2o.lon,0,xhdim-1);
   free_matrix(aero.aod,0,xadim-1,0,yadim-1);
   free_matrix(aero.ssa,0,xadim-1,0,yadim-1);
   free_matrix(aero.gg,0,xadim-1,0,yadim-1);
   free_vector(aero.lat,0,yadim-1);
   free_vector(aero.lon,0,xadim-1);
   free_matrix(sal.val,0,xsdim-1,0,ysdim-1);
   free_vector(sal.lat,0,ysdim-1);
   free_vector(sal.lon,0,xsdim-1);

   
   /* printf("before tcdlout \n");
      tcdlout();
      printf("after tcdlout \n"); 
   */
    printf("for regrid the cdl file %s \n", cdlfile);
    if (regrid != 0 && cloudcalc != 0)
     {
       /* define the regular grid, irradiance data will be regridded to the
          regular grid and .. */
       xrdim=(int)(fabs((lonb-lone)/dlon))+6; 
       yrdim=(int)(fabs((latb-late)/dlat))+6;
       printf("xrdim= %d, yrdim=%d \n", xrdim, yrdim);
       lonreg=vector(0,xrdim-1);
       latreg=vector(0,yrdim-1);
       
       for (i=0;i<xrdim;i++){
	 lonreg[i]=lonb+(i-2)*dlon;
       }
       for (j=0;j<yrdim;j++){ 
	 latreg[j]=latb+(j-2)*dlat;
	 /* printf("latreg %f \n", latreg[j]); */
       }
       printf("latb= %f latreg[0]= %f, latreg[yrdim-1]=%f \n", latb, latreg[0], latreg[yrdim-1]);
       printf("start to regrid, regrid=%d \n",regrid);
       printf("cdlfile= %s \n", cdlfile);
       fregrid2(lonreg,latreg,xrdim,yrdim,xndim,yndim,cloud.g,cloud.b,cdlfile);
       /* fregrid(lonreg,latreg,xrdim,yrdim,xdim,ydim,clear.g,clear.b); */
       printf("after regrid \n");
       free_vector(lonreg,0,xrdim-1);
       free_vector(latreg,0,yrdim-1);
     }
    if(cloudcalc == 0)
      {
       printf("write clear sky output xdim= %d, ydim=%d \n", xdim, ydim);
       lonreg=vector(0,xdim-1);
       latreg=vector(0,ydim-1);
       for (i=0;i<xdim;i++){
	 lonreg[i]=lonb+i*dlon;
       }
       for (j=0;j<ydim;j++){ 
	 latreg[j]=latb+j*dlat;
	 /* printf("latreg %f \n", latreg[j]); */
       }
       printf("before cdlout \n");
	fcdloutCS(lonreg,latreg,xdim,ydim,clear.g,clear.b,cdlfile);
       printf("after ...\n");
	free_vector(lonreg,0,xdim-1);
	free_vector(latreg,0,ydim-1);
      }


   /*
     for (i=1;i<ydim;i++){
     for (j=1;j<xdim;j++){
     if (nfile != 0 && cloudcalc != 0){
     printf("the final cloud results %d %d %f %f %f \n",i,j,clear.g[i][j],cloud.g[i][j],ni)
     }
     else{
     printf("the final clear results %d %d %f  \n",i,j,clear.g[i][j])
     }
     }
     }
   */

   int nrbytes=4;
   if (nfile != 0 && cloudcalc != 0){
     printf("write binary clod file \n");
     outGcld = fopen(outfile,"w");
     fwrite(&cloud.g[0][0],nrbytes,xdim*ydim ,outGcld);
     fclose(outGcld);
  }
   else{
     printf("write binary clear sky file \n"); 
     outGclear = fopen(outfile,"w");
     fwrite(&clear.g[0][0],nrbytes,xdim*ydim ,outGclear);
     printf("finishes \n");
     fclose(outGclear);
   }

    /* free the memory of the remaining dynamic arrays that have been allocated */   
    /* -> the start */
    free_matrix(clear.g,0,ydim-1,0,xdim-1);
    free_matrix(clear.b,0,ydim-1,0,xdim-1);
    if (nfile != 0 && cloudcalc != 0){
     free_matrix(cloud.g,0,ydim-1,0,xdim-1);
     free_matrix(cloud.b,0,ydim-1,0,xdim-1);
    }

  /* -> the end */


}


 /*   function:         borders

     PURPOSE:        searching for the location of x and provide the 
                     borders of x within the discret x-array xa(n)

     CALLING MODE:   function call

     PARAMETER DESCRIPTION:
     INTEGER  klo, khi, k   -> counters
     REAL xa(n)   -> array containing the discrete x values   
     REAL x       -> value for which the borders xa(i) and
                    xa(i-1) are scanned
     INTEGER n    -> Dimension of the array xa() 
     INTEGER i    -> i=khi upper counter of the border
                     the interval

     PROGRAMMER:    based on numerical recipes, modifications 
                    Richard Mueller

     VERSION: 0.9         
     DATE:    27.07.2008 */


 int borders(float xa[], int n,float x){  
   /* independent on increase or decreas of X(i) !! */ 
     int  i, klo, khi, k;
     klo=0 ; /* has to be 1 and n for arrays starting at i=1 */
     khi=n-1 ;
     while (khi-klo > 1){
       k=(khi+klo)/2;
       if (xa[klo] < xa[khi])  
	 {
	   if (xa[k] > x)
	     {khi=k;}
	   else
	     {klo=k;}
	 }
       else
	 { 
	   if (xa[k] < x)
	     {khi=k;}
	   else
	       {klo=k;}
	   
	 }
       
     }
     if ((khi-1) != klo) 
       { printf("WARNING problems in borders subroutine khi-1 NEQ klo");}
     /* only for testing ------------------------------------->
	if ((k .NE. khi) .AND. (k .NE. klo )) then
	print*, 'WARNING problems in borders k is not eq klo or khi'
	print*, 'k,khi,klo', k,khi,klo 
	endif */
     i=khi;
     return i;    
   }
   

int oborders(float xa[], int n,float x)
{  
  int  i, klo, khi, k;
  klo=0 ; /* has to be 1 and n for arrays starting at i=1 */
  khi=n-1 ;
  while (khi-klo > 1){
    k=(khi+klo)/2;
     if (xa[k] > x)
       {khi=k;}
     else
       {klo=k;}
     
  }
  if ((khi-1) != klo) 
    { printf("WARNING problems in borders subroutine khi-1 NEQ klo");}
/* only for testing ------------------------------------->
  if ((k .NE. khi) .AND. (k .NE. klo )) then
   print*, 'WARNING problems in borders k is not eq klo or khi'
   print*, 'k,khi,klo', k,khi,klo 
 endif */
  i=khi;
  return i;    
  }



float lint(float xa[],float ya[], int n, float x)
{  /* performs a linear interpolation */
  int  i, klo, khi, k;
  float dG,w;
  klo=0 ;
  khi=n-1 ;
  while (khi-klo > 1){
    k=(khi+klo)/2;
    if (xa[k] > x)
      {khi=k;}
    else
      {klo=k;} 
  }
  if ((khi-1) != klo) 
    { printf("WARNING problems in borders subroutine khi-1 NEQ klo");}
  /* only for testing ------------------------------------->
     if ((k .NE. khi) .AND. (k .NE. klo )) then
     print*, 'WARNING problems in borders k is not eq klo or khi'
     print*, 'k,khi,klo', k,khi,klo 
 endif */
  w=(xa[khi]-x)/(xa[khi]-xa[klo]);
   if (w >= 0 && w <= 1){     
     /*   printf("that is dgh20 %f %d %f %d %f \n",x, khi,xa[khi], klo,xa[klo]); 
	  printf("that is dgh20 %f %d %f %d %f \n",w, khi,ya[khi], klo,ya[klo]); */
   }
   else
   {printf("problem in function lint w is not between 1 an 0 !!!");  dG=-1000;}
  dG=w*ya[klo]+(1-w)*ya[khi];
  return dG;    
}

float interpolate2(float lon,float lat,int xdim, int ydim,float **var, float x[], float y[])
     /* 2 dimensional interpolation */
     /* parameter:
     wlatl,wlonl,wlath,wlatl    -> the weighting factors    
     latl,lonl,lath,latl        -> the associated lat and lon borders
     lat,lon                    -> latitude, longitude for which the weighting
                                        factors are calculated 
     programmed by Richard Mueller - Free University of Palatinum
 
     version: 0.2
     revision history:
     30.07.2007 version 0.1, does not use lat,lon fields but assumes
      specific start end of lat,lon gird and caluclates lat(i),lon(i)
      with start + k*ddeg (ddeg=resolution of grid)
     30.04.2008 version 0.2  
     uses lat lon fields of climatologies, more flexibel now 
     works on increasing or decreasing grids !!*/
{  float wlonh,wlonl,wlath,wlatl;
   /*,ddeg;*/
   float lonh,lonl,lath,latl;
   /*,lat,lon; */
   float weighting,dummi,varint ;
   int i,j,nlonl,nlonh,nlath,nlatl;         
   /*   printf("xdim= %d, ydim= %d \n",xdim,ydim);
      for(i=0;i<xdim;i++){
      for(j=0;j<ydim;j++){
         printf("h20 climi %f %f %f \n", x[i],y[j],var[i][j]) ;
	 }} */
   /* the borders and weighting for longitute */ 
     nlonh=borders(x,xdim,lon);
     nlonl=nlonh-1;
     wlonl=1-(lon-x[nlonl])/(x[nlonh]-x[nlonl]);
     wlonh=1-wlonl;
   /* the borders and weighting for latitude */ 
     nlath=borders(y,ydim,lat);
     nlatl=nlath-1;
     wlatl=1-(lat-y[nlatl])/(y[nlath]-y[nlatl]);
     wlath=1-wlatl;			   
     /* print commands for the interpolation test */
     /* printf("nlonl,...%d %d %d %d \n",nlonl,nlonh,nlatl,nlath);
     printf("lon= %f lat= %f, xl= %f xh= %f yl= %f yh= %f varl= %f varh=%f varix= %f \n"
     ,lon,lat,x[nlonl],x[nlonh],y[nlatl],y[nlath],var[nlonl][nlatl],var[nlonh][nlath], var[76][12]);
     printf("wlonl... %f %f %f %f \n",wlonl,wlonh,wlatl,wlath );*/
     /*   if ((khi-) 
    { printf("WARNING problems in borders subroutine khi-1 NEQ klo");
     hier Warnung für ws einbauen !!!*/
    dummi=wlonl*x[nlonl]+wlonh*x[nlonh]; /* dummy, only for testing the weighting factors */
   /*printf("dummilon= %f %f\n", dummi);*/
       dummi=wlatl*x[nlatl]+wlath*x[nlath];  

       /* the same for lat....   */
   /* printf("dummilat= %f \n", dummi); */
   varint=wlath*wlonl*var[nlonl][nlath]+wlatl*wlonl*var[nlonl][nlatl]
     +wlath*wlonh*var[nlonh][nlath]+wlatl*wlonh*var[nlonh][nlatl];
  

   /* printf("wlonlat= %f %f %f %f\n",wlonh,wlonl,wlath,wlatl);
   printf("nlonlat= %d %d %d %d\n",nlonh,nlonl,nlath,nlatl);
   printf("varint= %f %f \n",varint,var[70][28] );*/
     return varint;
}


void read_something(int xdim, int ydim, int month, char clim[], int nr)
{ /* dummy function to read somethig */
  int i,j,k;
  float tmp[13];
  char header[120];
  /*mad thing, when hclimf is read after aerosol an error message occur,
   file not found  */
  FILE *climf;   /* open the file containing the h2o climatology */
  climf = fopen(clim,"r"); 
  if (climf == NULL){
    printf("ERROR: user specific climatology seems to be not available \n");
  }
  else
    {some.lat=vector(0,ydim-1);
     some.lon=vector(0,xdim-1);
     some.val=matrix(0,xdim-1,0,ydim-1);
}
  /* read the data, the file can have nr columns, but
     only one is processed, hence dummy values are used as buffer */
  /* no header assumed 
    fscanf(hclimf,"%120s \n", &header); 
     printf("header %120s \n", header); 
     printf("position after header: %ld\n", ftell(hclimf));
     printf("Hallo %d \n", xhdim); */
  for(i=0;i<xdim;i++)
    { for(j=0;j<ydim;j++)
      {   /* warning, no fscanf of i,j if i,j are the counter !! */
	  fscanf(climf, "%f %f", &some.lon[i],&some.lat[j]);
	  for(k=1;k<=nr;k++){
	    fscanf(climf, "%f", &tmp[k]);
	    if (k==month){
	     some.val[i][j]=tmp[k];
	    }
	  }     
	  /*   printf("h20 clim %f %f %f \n", h2o.lon[i],h2o.lat[j],h2o.val[i][j]) ; */
	}
    }
  if (fclose(climf) != NULL){
    printf("ERROR: user specific climatology file can not be closed \n");
  }



  /*#include <stdio.h> */
  printf("finished program read the something\n") ; 
}

void rhclim(int xhdim, int yhdim, int month, char hclim[])
{  /* reads the water vapour climatology */
  int i,j,k;
  float tmp[13];
  char header[120];
  /*mad thing, when hclimf is read after aerosol an error message occur,
   file not found  */
  FILE *hclimf;   /* open the file containing the h2o climatology */
  hclimf = fopen(hclim,"r"); 
  /*hclimf = fopen("/home/richard/archive/software/4nsolis/climatologies/rhl1-ltmm-ncep.dat","r");*/
  if (hclimf == NULL){
    printf("ERROR: H2O climatology seems to be not available \n");
  }
  else
    {h2o.lat=vector(0,yhdim-1);
     h2o.lon=vector(0,xhdim-1);
     h2o.val=matrix(0,xhdim-1,0,yhdim-1);
}
  /* read the H20 data, the H20 file contains 12 month, but
     only one is processed, hence dummy values are used as buffer */
     fscanf(hclimf,"%120s \n", &header); /* scan the header */
     printf("header %120s \n", header); 
     printf("position after header: %ld\n", ftell(hclimf));
   printf("Hallo %d \n", xhdim);
  for(i=0;i<xhdim;i++)
    { for(j=0;j<yhdim;j++)
      {   /* warning, no fscanf of i,j if i,j are the counter !! */
	  fscanf(hclimf, "%f %f", &h2o.lon[i],&h2o.lat[j]);
	  for(k=1;k<=12;k++){
	    fscanf(hclimf, "%f", &tmp[k]);
	    if (k==month){
	      h2o.val[i][j]=tmp[k];
	    }
	  }     
	  /*   printf("h20 clim %f %f %f \n", h2o.lon[i],h2o.lat[j],h2o.val[i][j]) ; */
	}
    }
  if (fclose(hclimf) != NULL){
    printf("ERROR: H2O climatology file can not be closed \n");
  }
}


void raeroclim(int xadim, int yadim, int month, char aeroclim[]){ 
  /* reads the aerosol climatology */
float tmp[13];
 char header[120];
  int i,j,k;
  FILE *aeroclimf;   /* open the file containing the aerosol climatology */
  aeroclimf = fopen(aeroclim,"r");
  if (aeroclimf == NULL){
    printf("ERROR: Aerosol climatology seems to be not available \n");
  }

  /* read the aerosol data, the aerosol file contains 12 month, but
     only one os processed, hence dummy values are used as buffer */
  aero.lat=vector(0,yadim-1);
  aero.lon=vector(0,xadim-1);
  aero.ssa=matrix(0,xadim-1,0,yadim-1);
  aero.aod=matrix(0,xadim-1,0,yadim-1);
  aero.gg=matrix(0,xadim-1,0,yadim-1); 
  fscanf(aeroclimf,"\n %120s \n", &header); /* scan the header */
  printf("header %120s \n", header); 
  printf("position after header: %ld\n", ftell(aeroclimf));
  for(i=0;i<xadim;i++)
    {
      for(j=0;j<yadim;j++)
        /* tha assymetry parameter is set to a constant value 
	   of 0.7*/
	{ 
          aero.gg[i][j]=0.7;  
	  fscanf(aeroclimf, "%f %f", &aero.lon[i],&aero.lat[j]);
	  for(k=1;k<=12;k++){
	    fscanf(aeroclimf, "%f", &tmp[k]);
	    if (k==month){
	      /* aero.aod[i][j]=tmp[k];*/
	      aero.aod[i][j]=tmp[k];  
	    }
	  }    
	  for(k=1;k<=12;k++){
	    tmp[k]=-9999.0;
	    fscanf(aeroclimf, "%f", &tmp[k]);
	    if (k==month){
	      aero.ssa[i][j]=tmp[k];
	    }
	  }  
	  /* printf("%d %d  %f %f ", i,j, aero.lon[i],aero.lat[j]);  
	     printf("%f %f %f\n", aero.aod[i][j],aero.ssa[i][j],aero.gg[i][j] );     */
	}
    }
  if (fclose(aeroclimf) != NULL){
    printf("ERROR: aerosol climatology file can not be closed \n");
  }
  /*  fclose(aeroclimf); */
}


void rclut(int ngg, int nssa , int naod, int nh2o, int no3, char clut[]){
  /* reads the clear sky look up table */
  char header[100];
  int i,j,k;
  FILE* clutf; /* open the file containing the look-up tables */
  char line[100];
  clutf = fopen(clut,"r"); 
  /*hclimf = fopen("/home /richard/archive/software/4nsolis/climatologies/rhl1-ltmm-ncep.dat","r");*/
  if (clutf == NULL){
    printf("ERROR: Clear sky LUT seems not to be available \n");
  }
  /* read the H20 data, the H20 file contains 12 month, but
     only one is processed, hence dummy values are used as buffer */
  /*fgets(header,sizeof(header),fclut); */ /* scan the header */    
  fgets(header,sizeof(header),clutf);
  printf("header %70s \n", header); 
  printf("in LUT position after header: %ld\n", ftell(clutf));
  for(i=0;i<ngg;i++){ 
    for(j=0;j<nssa;j++){  
      fscanf(clutf," %f %f",&alut.ssa[j],&alut.gg[i] );
      /* fscanf(clutf,"%s %f %s %f",&bla1,&alut[i][j].ssa,&bla2,alut[i][j].gg);*/
      printf("ssa,gg: %f %f \n", alut.ssa[j], alut.gg[i]);
      /* warning, no fscanf of i,j if i,j are the counter !! */
      printf("BEFORE k loop \n");
      /*printf("yes \n"); */
      for(k=0;k<naod;k++){
	/* printf("in k loop \n"); */ 
	fscanf(clutf, "%f %f %f %f %f %f \r", &alut.aod[k], &alut.Im[i][j][k],
	       &alut.gtau[i][j][k],&alut.ag[i][j][k],
	       &alut.btau[i][j][k],&alut.ab[i][j][k]);
	printf("the alut %d %d %d %f %f \n",i,j,k,alut.aod[k],alut.Im[i][j][k]);       
	printf("%f %f ", alut.gtau[i][j][k],alut.ag[i][j][k]);
	   printf("%f %f \n", alut.btau[i][j][k],alut.ab[i][j][k]);                       
      }     
      /*printf("%f %f %f \n", h2o.lon[i],h2o.lat[j],h2o.val[i][j]);*/
    }
  }
  fgets(header,sizeof(header),clutf);
  printf("line %70s \n", header); 
  /*printf("position after header: %ld\n", ftell(clutf)); */
  fgets(line,sizeof(line),clutf);
  printf("line %70s \n", line); 
  for(k=0;k<nh2o;k++){
    fscanf(clutf, "%f %f %f \r", &hlut.x[k], &hlut.yg[k], &hlut.yb[k]);
    printf("the h20 params %d %f %f %f \n",k,hlut.x[k], hlut.yg[k],hlut.yb[k]);                    
  }     
  fgets(line,sizeof(line),clutf);
  printf("line %70s \n", line); 
  for(k=0;k<no3;k++){
    fscanf(clutf, "%f %f %f \r",&olut.x[k], &olut.yg[k], &olut.yb[k]);
    printf("the O3 params %f %f \n",olut.yg[k], &olut.yg[k]);                              
  }    
  if (fclose(clutf) != NULL){
    printf("ERROR: file containing clear lut can not be closed \n");
  }  
}




float normsal(float dval, float u0, float sal)
{
  /*
    u0=consine of the solar zenith angle
    sala= the surface albedo corrected to AM2, SZA=60 degree 
        relative to the SZA of the CAVE value at AM 2, SZA=60 degree */
  float sala;
  float corsza;
  corsza=(1.0+dval)/(1.0+2.0*dval*u0);
  sala=sal * corsza;
  return sala;
}

float corsal(short int igbp, float sza)
{
  float corsza;
  float dval[20];
  dval[0]=0.4;      /* Evergreen Needle Forrest, IGBP class 1 */
  dval[1]=0.44;     /* Evergreen Broad Forrest, IGPB class 2 */
  dval[2]=0.32;     /* Deciduos Needle Forrest, IGPB class 3 */
  dval[3]=0.39;     /* Deciduos Braod Forrest, IGPB class 4 */  
  dval[4]=0.22;     /* Mixed Forest, class 5 */
  dval[5]=0.28;     /* Closed Shrubs, class 6 */
  dval[6]=0.40;     /* Open Shrubs, class 7 */
  dval[7]=0.47;     /* Woody savanna, class 8 */
  dval[8]=0.53;     /* Savanna, class 9 */
  dval[9]=0.53;     /* Grassland, class 10 */
  dval[10]=0.35;    /* Wetland, class 11 */
  dval[11]=0.41;    /* Cropland, class 12 */
  dval[12]=0.10;    /* Urban, class 13 */
  dval[13]=0.40;    /* Crop Mosaic, class 14 */
  dval[14]=0.10;    /* Anatartic Snow, class 15 */
  dval[15]=0.40;    /* Barren/Desert, class 16 */
  dval[16]=0.41;    /* Ocean Water, class 17 */
  dval[17]=0.58;    /* Tundra, class 18 */
  dval[18]=0.10;    /* Fresh Snow, class 19 */
  dval[19]=0.10;    /* Sea Ice, class 20 */
  corsza=(1.0+dval[igbp-1])/(1.0+2.0*dval[igbp-1]*cos(sza));
  return corsza;
}
 
int rsarbsal(char salclim[], char luclim[], int xsdim, int ysdim)
{ /* reads the SARB/CERES surface albedo*/
  int i,j,k;
  int salp;
  float sza=1;
  float dummy=0.0;
  char dum1[xsdim*ysdim],dum2[xsdim*ysdim];  
  /* to do provide xdim,ydim with function call, instead of using fixed values */ 
  float res;
  FILE *salclimf; 
  FILE *luclimf;
  int rsize; 
 /* open the file containing the sal climatology */
  salclimf = fopen(salclim,"r");
  if (salclimf == NULL)
    {
    printf("ERROR: SAL climatology seems to be not available \n");
    }
  luclimf = fopen(luclim,"r");
  if (luclimf == NULL)
    {
    printf("ERROR: Land Use climatology seems to be not available \n");
    }

  /* read the LARC NASA albedo map from the SARB (Surface and Atmospheric radiation) 
     working group, climatology, values normalised to 60 degree */
  /* 10 minute, (1/6 of a degree) resoultion on an equal angle map covering the globe */
  /* first test , ASCII will be changed !! */
  sal.lat=vector(0,ysdim-1);  /*regular grid, hence lat/lon have not to be a matrix */
  sal.lon=vector(0,xsdim-1);
  sal.val=matrix(0,xsdim-1,0,ysdim-1);
  sal.lu=bmatrix(0,xsdim-1,0,ysdim-1);

  k=xsdim*ysdim;
  rsize=fread((char*)&dum1,1,k,salclimf);
  if (rsize != k) {printf("problem reading albedo map, size do not match \n");}
  /* read the IGBP landuse information provided by SARB in order to
     be able to calculate SAL for varying SZA */ 
  rsize=fread((char*)&dum2,1,k,luclimf);
  if (rsize != k) {printf("problem reading landuse map, size do not match \n");}

 for(j=0;j<ysdim;j++) /* alt xsdim, here fixed */
    {
    for(i=0;i<xsdim;i++)  /*alt ysdim, here fixed */
	{
       /* the command for the ascii file... fscanf(salclimf, "%d", &salp); */
	  k=i+j*xsdim;
          /* for transformation not needed is considered in main ii=1080+i
	     if (i > 1080) {ii=i-1080;} */
	  /*for ascii	sal.val[i][j]=salp/100.0; */
          sal.val[i][j]=dum1[k]/100.0; 
          sal.lu[i][j]=dum2[k];
       }
    /*  printf("%d %d %f %d \n", i,j, sal.val[i][j],dum1[k]); */  
    }    

  /* printf("%d %d  %f %f ", i,j, sal.lon[i],sal.lat[j]);  
     printf("%f %f %f\n", sal.aod[i][j],sal.ssa[i][j],sal.gg[i][j] );     */
  res=1/6;
  for(i=0;i<xsdim;i++)
    {
      dummy=(float)i/6;
      sal.lon[i]=dummy;
      /* if (dummy >= 180.0 ){
	 sal.lon[i]=-360+dummy;  not needed any more, shift to main prgram 
      } */
      /*  printf("i %d lon %f dummy %f \n",i,sal.lon[i],dummy); */
    }  
  dummy=0.0;
  for(j=0;j<ysdim;j++)
    {
      dummy=(float)j*1/6;
      sal.lat[j]=90.0-(float)(1/12)-dummy; 
      /* printf("j %d lat %f dummy %f \n",j,sal.lat[j],dummy); */
    }  
  
  /* the correction is now done in the main program !! */
  /*for(i=0;i<xsdim;i++){
    for(j=0;j<ysdim;j++){
      printf("bef lon %f lat %f sal %f lu %d \n",sal.lon[i],sal.lat[j],sal.val[i][j],sal.lu[i][j]); 
    }
  }*/


  if (fclose(salclimf) != NULL){
    printf("ERROR: sal climatology file can not be closed \n");
  }
  /*  fclose(salclimf); */


}

#include "nr_magic.c"
#include "nrutil.c"
#include "magic-geo.c"
