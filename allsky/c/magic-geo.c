/* this c program compiles differnt modules for the calculation of 
geometry parameters, geolocation, scattering angle calucaltion and
time functions */


/* MODULE geoloc:
   the module geoloc is based on the Eumetsat routine refgeo in fortran,
   it has been rewritte in C andmodified by R.W. Mueller 
   !!!! works currrently only for satellite positioned at 0 degree longitude 
 
Purpose:  the program converts line and column of Meteosat images to latitude and
          longitude

Input parameter:  line and column of the pixel, and numnber of lines of the image (nrlines)  
                  line measured from southern end of frame !
                  column measured from eastern end o frame !
                  in float in order to enable sub-pixel accuracy.
                  Integer values correspond to the middle of the pixel.

Output parameter: lati     -latitude of the pixel (degrees North from Equator) 
                  long 180.0/pi
#define h 4216g     -longitude of the pixel (degrees East from Greenwich) 
                  disk  -flag, set to 1 (true) if the pixel is on the disk.
                                  set to 0 (false) if the pixel is in space. 


Version: 0.1, based on refgeo from Eumetsat.*/




int geoloc(int line, int column, int nrlines, double *lat_ptr, double *lon_ptr)
{
#define altitude 42164.0   /* distance from earth centre to satellite */
#define req 6378.140      /* radius of the earth at the equator and pole */
#define rpol 6356.755 
#define oblate 1.0/298.257 /* earth oblateness  */
/* #define pi 3.141592653   */



    double cenlat,aline,asamp,tanal,tanas,det,k,nrcols;  /* nrcols = number of columns */
    double step,x,y,z,a,b,c,p,q,r;
    short int disk; /* pixel is on the disk (1) or not (0) */
    float subres=2.4;

/* step is the radiometer step as seen by the spacecraft in degrees. The image (whole disk) represents
an 18x18 degree field of view divided up on an equi-angulat basis. This has to be divided by the line/col 
numbers of the image. For the whole disk number of lines and columns are equal, step 
reduces to 18/nrlines */
     
	/* step=18.0/nrlines; orig */
    step=aperture/nrlines;
    /* alternatively atan step  = atan(subres/(altitude-req)) */
    /* step  = atan(subres/(altitude-req))*rad2deg; */
    /* tan (step) * (altituide-req) = sub-satellite resolution in km */  
/* !! line 

/* convert line/column values to angular offsets from center points */
/* the image is assumed to be quadratic with nrlines=nrcols */
    nrcols=nrlines;
    asamp=-(column - 0.5 * nrcols+0.5) * step;  /* orig if measured from eastern end of frame !!*/
   /* tan (step) * (altituide-req) = sub-satellite resolution in km */
    aline=(line - 0.5 * nrlines+0.5) * step;
/* convert degree to radians as the trigonometric functions of C aspects radians */
    asamp=asamp * deg2rad;
    aline=aline * deg2rad;
/* Calculate tangents of angle offsets */
    tanal=tan(aline);
    tanas=tan(asamp);
/* calculate the components of the vector from the spacecraft in the viewing direction of the respective pixel */
    p=-1;
    q=tanas;
    r=tanal * sqrt (1.+ q*q);
/*    printf("pow,(),q,r= %f %f %f %f \n",pow((r * req/rpol),2),(req/rpol),q,r); */

/* the location of the pixel/point on the earth can be identified by solving a quadratic equation for the intersection 
 between the earth´s surface and the viewing line from the spacecraft.
 if there exist no real roots then there is no intersection, otherwise the required root is the one nearer to the
 spacecraft */

    a=q * q + pow((r * req/rpol),2) + p*p;
    b=2. * altitude * p ;
    c=altitude * altitude - req * req;
 
/* Calculate the determinant. If it is negative (no real roots of the quadratic equations) there is no intersection between the line
   of sight and the disk and the pixel is not on the disk */

    det= b * b - 4 * a* c;
    /*
    printf("p,q,r= %f %f %f %f \n",p,q,r,det);   
    printf("p,q,r= %f %f %f %f %f\n",a,b,c,det,column,line);  */
    if(det < 0.) 
    {
         
	disk=0;
        *lat_ptr=999;
        *lon_ptr=999;
           return(1); 
    }
   else
   {
       /* printf("disk =0 \n"); */
       disk=1;
       k=(-b -sqrt(det))/(2*a);
       x=altitude + k*p;
       y=k*q;
       z=k*r;
       *lon_ptr=-atan(y/x);   /* orig *lon_ptr=atan(y/x); if column is measured from eastern end of frame as for openMTP !!*/
       cenlat=atan(z*cos( *lon_ptr)/x);
       /* this is the geocentric latitude. It is converted to the geodetic (geographic) latitude */
       *lat_ptr=-atan( tan(cenlat)/pow((1.0-oblate),2) ); /* orig *lat_ptr=-atan(..) if lines are measured from southern end 
							     of frame ); as for openMTP !*/
       /* convert lat,lon back from radians to degree is not performed as Helisat assumes lat & lon in radians 
	*lat_ptr=*lat_ptr*rad2deg;
        *lon_ptr=*lon_ptr*rad2deg;
        */ 
	*lat_ptr=*lat_ptr*rad2deg;
        *lon_ptr=*lon_ptr*rad2deg+satpos;
 
   return(0); 
   } 
}



/*------------------------------------------------------------------------------------------------*/
double equation_of_time(int doy, double *gamma, double *f, double *dec)
{
  /* Equation of time in minutes: 
     Difference of sun position at 12 GMT and south position
     see: Iqbal, An Introduction to Solar Radiation  */

  double E_t;
  double cosgamma,cos2gamma,singamma,sin2gamma;
  
  *gamma = 2*PI* (doy - 1.0)/365.0;                    /* Tageswinkel in Rad */
  cosgamma=cos( *gamma);
  cos2gamma=cos(2* *gamma);
  singamma=sin( *gamma);
  sin2gamma=sin(2* *gamma);
  /* f gives the correction of irradiance for sun earth distance !! */

  *f = 1.00010 + 0.034221*cosgamma                  /* Erdbahn als Ellipse in Rad  */
    + 0.001280*singamma              
    + 0.000719*cos2gamma
    + 0.000077*sin2gamma;
  *dec = (0.006918 - 0.399912*cosgamma              /*  Deklination in Rad */
      + 0.070257*singamma
      - 0.006758*cos2gamma
      + 0.000907*sin2gamma
      - 0.002697*cos(3* *gamma)
      + 0.00148*sin(3* *gamma));
  E_t = (0.000075 + 0.001868*cosgamma               /* E_t in Minuten  */
      - 0.032077*singamma
      - 0.014615*cos2gamma
      - 0.04089*sin2gamma) * (4.0 *180.0/PI);
      
  return(E_t);
}

/*------------------------------------------------------------------------------------------------*/
/* no modifications needed for MFG */
double cos_zenit(double lat,double lon,double dec,double lon_sun)
		/* Cosine of the solar zenith angle */
{
double w;                                   /* hour angle */
double(cos_theta_zenith);

w = (lon -lon_sun);
cos_theta_zenith =sin(lat)*sin(dec)          /* see Iqbal */
		 + cos(lat)*cos(dec)*cos(w); 
if(cos_theta_zenith<=0.0)
  {
   cos_theta_zenith=0.0;
   return(0);
  }
return(cos_theta_zenith);
}

/* RM May 2008: Added function for the calculation of the solar azimuth angle */
/* no modifications needed for MFG */
double cos_azimuth(double lat,double lon, double dec,double sza,double lon_sun)
		/* cosine of the azimuth angle  */
{
double w;                                   /* hour angle*/
w=lon-lon_sun;
double(cos_solar_azimuth);
double x_azm,y_azm;
x_azm=sin(w)*cos(dec);
y_azm=(cos(lat)*sin(dec)-cos(w)*cos(dec)*sin(lat));
cos_solar_azimuth=cos(atan(x_azm/y_azm));

/* cos_solar_azimuth = (sin(sza)*sin(lat)-sin(dec))/(sin(sza)*cos(lat)); */
         /* approximation see wikipedia  cos(elevation)=sin(pi/2-elevation)= sin(sza)*/
/* cos_solar_azimuth = ((sin(dec)-cos(sza)*sin(lat))/(sin(sza)*cos(lat))) ;*/

/* if(w < 0.0)
  {
      cos_solar_azimuth=-1.0 * cos_solar_azimuth;
   return(0);
   } */
return(cos_solar_azimuth);
}


double azimuth(double lat,double lon, double dec,double sza,double lon_sun)
		/* cosine of the azimuth angle  */
{
double w;                                   /* hour angle*/
w=lon-lon_sun;
double(solar_azimuth);
double x_azm,y_azm;
x_azm=sin(w)*cos(dec);
y_azm=(cos(lat)*sin(dec)-cos(w)*cos(dec)*sin(lat));
solar_azimuth=atan(x_azm/y_azm);

/* cos_solar_azimuth = (sin(sza)*sin(lat)-sin(dec))/(sin(sza)*cos(lat)); */
         /* approximation see wikipedia  cos(elevation)=sin(pi/2-elevation)= sin(sza)*/
/*if(cos_solar_azimuth<=0.0)
  {
   cos_solar_azimuth=0.0;
   return(0);
   } */
return(solar_azimuth);
}

/* RM May 2006: Added function for the calculation of the scattering angle */

double scatt_angle(double sza, double azi, double zen, double sataz)
{
    double (scangle);
    double cosscangle;
    double relazi;
    relazi=azi+sataz;
    /* = cos(SAZ)*cos(sat-zenith-angle)+sin(sza)*sin(sat-zenith-angle)*cos(relative azimut) */
    cosscangle=cos(sza)*cos(zen)+sin(sza)*sin(zen)*cos(relazi);
    scangle=acos(cosscangle);
    return (scangle);
}



void xpifinfo(unsigned char *headline,struct bildinfo *temp)
/* liest aus dem XPIF-Header alle noetigen Informationen 
   und gibt sie als Struktur zurueck zur weiteren Verarbeitung */
{  
    (*temp).nrbytes=*((unsigned char *)(&headline[0])); /* Number of bytes per pixel  */
    (*temp).nrrecords=*((unsigned short *)(&headline[1]));
    (*temp).nrlines=*((unsigned long *)(&headline[3])); 
    (*temp).nrcolumns=*((unsigned long *)(&headline[7]));
    (*temp).nrcolbig=(*temp).nrrecords*RECLENG/(*temp).nrbytes; 
    /* image size with black frame, in case that nrcolums not equal
       nrrecords*RECLENG */

    (*temp).add_header=*((unsigned char *)(&headline[11]));
    (*temp).pif_comp=*((unsigned char *)(&headline[12]));
    (*temp).pif_ahx=*((unsigned short *)(&headline[13]));
    (*temp).pif_nb=*((unsigned char *)(&headline[15]));

    (*temp).nav_func=*((unsigned char *)(&headline[16]));
    (*temp).line_off=*((signed long *)(&headline[17]));
    (*temp).col_off=*((signed long *)(&headline[21]));
    (*temp).nav_lres=*((signed long *)(&headline[25]));  /* Line  Resolution  */
    (*temp).nav_cres=*((signed long *)(&headline[29]));  /* Column Resolution  */
    /* RM May 2008: the following line is commented out in the eumetsat version */
    /* strncat((*temp).nav_data,&headline[33],95); */

    (*temp).dat_f=*((unsigned char *)(&headline[128]));

    (*temp).dat_toff=*((signed short *)(&headline[190]));

    (*temp).src_year=*((signed short *)(&headline[208]));

    (*temp).doy=*((unsigned short *)(&headline[210]));   /* day of year  */
    (*temp).aq_time=*((unsigned short *)(&headline[212])); /* aquisition time  */
   
}

int fregrid(float lonrev[], float latrev[], int xrdim, int yrdim, int xdim, int ydim, float **var1, float **var2)
{
  float  value=-100.0;
  float delta[xrdim][yrdim];
  float delta1[xrdim][yrdim];
  float **value1,**value2;
  int i,j,kk,ll,ii,jj;
  double lat,lon;
  FILE *cdlout;
  char outf[100]="testcdl.out";

  value1=matrix(0,xrdim-1,0,yrdim-1);
  value2=matrix(0,xrdim-1,0,yrdim-1); 

  for (ii=0;ii<xrdim;ii++){
    for (jj=0;jj<yrdim;jj++){
      delta1[ii][jj]=1000.0;
      delta[ii][jj]=1000.0;   
    }
  }
  /* usee memset above instead of loop !! */

  /*  printf ("start to search for the near neighbour \n"); */
    for (i=0;i<xdim;i++){
     for (j=0;j<ydim;j++){
      geoloc(j,i,ydim,&lat,&lon);
      if (lat >= (latrev[0]) && lat <= (latrev[yrdim-1])  && 
	  lon >= lonrev[0] && lon <= lonrev[xrdim-1]) 
	{ /*printf("regrid inside if \n");*/
	  jj=borders(latrev,yrdim,lat);
	  ii=borders(lonrev,xrdim,lon);
       	  if (ii >= 2 && ii < (xrdim-2) && jj >= 2 && jj < (yrdim-2))
	    { 
	      /* printf("regrid after borders %f %f %f %f %d %d %f %f %d %d\n ",lon, lat, lonrev[ii],latrev[jj],ii,jj,var1[i][j],var2[i][j],i,j);*/
	      for (kk=ii-2;kk<=ii+2;kk++){    
		for (ll=jj-2;ll<=jj+2;ll++){
		  delta[kk][ll]=sqrt(pow((lon-lonrev[kk]),2)+pow((lat-latrev[ll]),2));
		  if (delta[kk][ll] < delta1[kk][ll]){
		    delta1[kk][ll]=delta[kk][ll];
		    value1[kk][ll]=var1[i][j];
		    value2[kk][ll]=var2[i][j]; 
		  }
		}
	      }
	    }
	  else{ 
	    delta[ii][jj]=sqrt(pow((lon-lonrev[ii]),2)+pow((lat-latrev[jj]),2));
	    /*    !write(*,*) "outif", ii,jj,delta1(ii,jj),delta(ii,jj), sis(i,j)*/
	    if (delta[ii][jj] < delta1[ii][jj]){
	      delta1[ii][jj]=delta[ii][jj];
	      value1[ii][jj]=var1[i][j];
	      value2[ii][jj]=var2[i][j];   
	      /* write(8,*) "edge", ii,jj,delta1(ii,jj),delta(ii,jj), var1(i,j), value(k,ii,jj)*/
	    }
	  }
	  /* printf("regrid the vals %f %f \n",value1[ii][jj],value2[ii][jj]);*/
	}
     }

 } 

  printf("start to write cdl file \n");  
  cdlout=fopen("testcdl.out","w") ;
    fprintf(cdlout,"netcdf tmp \{ \n"); 
    fprintf(cdlout,"dimensions:\n");
    fprintf(cdlout,"lat = %d, lon= %d, time=unlimited; \n", yrdim,xrdim);
    fprintf(cdlout,"variables: \n");
    fprintf(cdlout,"float lat(lat), lon(lon);\n");
    fprintf(cdlout,"float Z(lat,lon);\n");
    fprintf(cdlout,"// variable attributes\n");
    fprintf(cdlout,"lat:long_name = \"latitude\";\n");
    fprintf(cdlout,"lat:units = \"degree\";\n");
    fprintf(cdlout,  "lon:long_name = \"longitude\";\n");
    fprintf(cdlout, "lon:units = \"degrees\";\n");
    fprintf(cdlout, "Z:units = \"Watt\";\n");
    fprintf(cdlout, "Z:valid_range = 0., 1400.;\n");
    fprintf(cdlout,"data:\n");
    fprintf(cdlout,"lat=\n"); 
    for(j=0;j<yrdim-1;j++)
      {fprintf(cdlout,"%f ,\n",latrev[j]);}
    fprintf(cdlout,"%f ;\n",latrev[yrdim-1]);
    fprintf(cdlout,"lon=\n"); 
    for(i=0;i<xrdim-1;i++)
      {fprintf(cdlout,"%f , \n",lonrev[i]);}
    fprintf(cdlout,"%f ; \n",lonrev[xrdim-1]);
    fprintf(cdlout,"Z=\n");
    for(j=0;j<yrdim;j++){
      for(i=0;i<xrdim;i++){
	if (j==yrdim-1 && i==xrdim-1) {
	fprintf(cdlout,"%f ;\n",value1[i][j]);
        }
        else
        fprintf(cdlout,"%f ,\n",value1[i][j]);	    
      }
    }
    fprintf(cdlout,"\}"); 
    fclose(cdlout);
  

    free_matrix(value1,0,xrdim-1,0,yrdim-1);
    free_matrix(value2,0,xrdim-1,0,yrdim-1);
    
}
 
int fregrid2(float lonrev[], float latrev[], int xrdim, int yrdim, int xdim, int ydim, float **var1, float **var2, char cdlfile[])
{
  float  value=-100.0;
  float **delta;
  float **delta1;
  float **value1,**value2;
  int i,j,kk,ll,ii,jj;
  double lat,lon;
  FILE *cdloutf;
  /* define and open the file for the cdl/netcdf output for region specified in config file */
  cdloutf = fopen(cdlfile,"w");

  value1=matrix(0,yrdim-1,0,xrdim-1);
  value2=matrix(0,yrdim-1,0,xrdim-1); 
  delta1=matrix(0,yrdim-1,0,xrdim-1);
  delta=matrix(0,yrdim-1,0,xrdim-1);



  /* rm feb 2011 initialise the regridded field   */

  for(j=0;j<yrdim;j++){
    for(i=0;i<xrdim;i++){
      value1[j][i]=-1;
      value2[j][i]=-1;    
    }
  }
  

  /* initialise the delta array used to defines the nearest neighbour*/ 
  for (ii=0;ii<yrdim;ii++){
    for (jj=0;jj<xrdim;jj++){
      delta1[ii][jj]=1000.0;
      delta[ii][jj]=1000.0;   
    }
  }
  
  /*  printf ("start to search for the near neighbour \n"); */
    for (i=0;i<xdim;i++){
     for (j=0;j<ydim;j++){
      geoloc(j,i,ydim,&lat,&lon);
      if (lat >= (latrev[0]) && lat <= (latrev[yrdim-1])  && 
	  lon >= lonrev[0] && lon <= lonrev[xrdim-1] && lat != 999 && lon != 999) 
	{ /*printf("regrid inside if \n");*/
	  ii=borders(latrev,yrdim,lat);
	  jj=borders(lonrev,xrdim,lon);

       	  if (jj >= 2 && jj < (xrdim-2) && ii >= 2 && ii < (yrdim-2))
	    { 
	      /* printf("regrid after borders %f %f %f %f %d %d %f %f %d %d\n ",lon, lat, lonrev[ii],latrev[jj],ii,jj,var1[i][j],var2[i][j],i,j);*/
	      for (kk=ii-2;kk<=ii+2;kk++){    
		for (ll=jj-2;ll<=jj+2;ll++){
		  delta[kk][ll]=sqrt(pow((lon-lonrev[ll]),2)+pow((lat-latrev[kk]),2));
		  if (delta[kk][ll] < delta1[kk][ll]){
		    delta1[kk][ll]=delta[kk][ll];
		    value1[kk][ll]=var1[j][i];
		    value2[kk][ll]=var2[j][i]; 
		  }
		}
	      }
	    }
	  else{ 
	    delta[ii][jj]=sqrt(pow((lon-lonrev[jj]),2)+pow((lat-latrev[ii]),2));
	    /*    !write(*,*) "outif", ii,jj,delta1(ii,jj),delta(ii,jj), sis(i,j)*/
	    if (delta[ii][jj] < delta1[ii][jj]){
	      delta1[ii][jj]=delta[ii][jj];
	      value1[ii][jj]=var1[j][i];
	      value2[ii][jj]=var2[j][i];   
	      /* write(8,*) "edge", ii,jj,delta1(ii,jj),delta(ii,jj), var1(i,j), value(k,ii,jj)*/
	    }
	  }
	  /* printf("regrid the vals %f %f \n",value1[ii][jj],value2[ii][jj]);*/
	}
     }

 } 


  printf("start to write cdl file \n");  
  /* cdlout=fopen("testcdl.out","w") ; */
    fprintf(cdloutf,"netcdf tmp \{ \n"); 
    fprintf(cdloutf,"dimensions:\n");
    fprintf(cdloutf,"lat = %d, lon= %d, time=unlimited; \n", yrdim,xrdim);
    fprintf(cdloutf,"variables: \n");
    fprintf(cdloutf,"float lat(lat), lon(lon);\n");
    fprintf(cdloutf,"float G(lat,lon), B(lat,lon);\n");
    fprintf(cdloutf,"// variable attributes\n");
    fprintf(cdloutf,"lat:long_name = \"latitude\";\n");
    fprintf(cdloutf,"lat:units = \"degree\";\n");
    fprintf(cdloutf,  "lon:long_name = \"longitude\";\n");
    fprintf(cdloutf, "lon:units = \"degrees\";\n");
    fprintf(cdloutf, "G:units = \"Watt per square meter\";\n");
    fprintf(cdloutf, "G:valid_range = 0., 1400.;\n");
    fprintf(cdloutf, "B:units = \"Watt per square meter\";\n");
    fprintf(cdloutf, "B:valid_range = 0., 1200.;\n");
    fprintf(cdloutf,"data:\n");
    fprintf(cdloutf,"lat=\n"); 
    for(j=0;j<yrdim-1;j++)
      {fprintf(cdloutf,"%f ,\n",latrev[j]);}
    fprintf(cdloutf,"%f ;\n",latrev[yrdim-1]);
    fprintf(cdloutf,"lon=\n"); 
    for(i=0;i<xrdim-1;i++)
      {fprintf(cdloutf,"%f , \n",lonrev[i]);}
    fprintf(cdloutf,"%f ; \n",lonrev[xrdim-1]);
    fprintf(cdloutf,"G=\n");
 
    for(j=0;j<yrdim;j++){
      for(i=0;i<xrdim;i++){
	if (j==yrdim-1 && i==xrdim-1) {
	fprintf(cdloutf,"%f ;\n",value1[j][i]);
        }
        else
        fprintf(cdloutf,"%f ,\n",value1[j][i]);	    
      }
    }
    
    fprintf(cdloutf,"B=\n");
    for(j=0;j<yrdim;j++){
      for(i=0;i<xrdim;i++){
	if (j==yrdim-1 && i==xrdim-1) {
	  fprintf(cdloutf,"%f ;\n",value2[j][i]);
        }
        else
	  fprintf(cdloutf,"%f ,\n",value2[j][i]);	    
      }
    }
    fprintf(cdloutf,"\}"); 
    fclose(cdloutf);
  

    free_matrix(value1,0,yrdim-1,0,xrdim-1);
    free_matrix(value2,0,yrdim-1,0,xrdim-1);
    
}


int fcdloutCS(float lonrev[], float latrev[], int xrdim, int yrdim,float **var1, float **var2, char cdlfile[])
{

  int i,j;
  FILE *cdloutf;
  cdloutf = fopen(cdlfile,"w");
  printf("start to write cdl file \n");  
  /* cdlout=fopen("testcdl.out","w") ; */
  fprintf(cdloutf,"netcdf tmp \{ \n"); 
  fprintf(cdloutf,"dimensions:\n");
  fprintf(cdloutf,"lat = %d, lon= %d, time=unlimited; \n", yrdim,xrdim);
  fprintf(cdloutf,"variables: \n");
  fprintf(cdloutf,"float lat(lat), lon(lon);\n");
  fprintf(cdloutf,"float G(lat,lon), B(lat,lon);\n");
  fprintf(cdloutf,"// variable attributes\n");
  fprintf(cdloutf,"lat:long_name = \"latitude\";\n");
  fprintf(cdloutf,"lat:units = \"degree\";\n");
  fprintf(cdloutf,  "lon:long_name = \"longitude\";\n");
  fprintf(cdloutf, "lon:units = \"degrees\";\n");
  fprintf(cdloutf, "G:units = \"Watt per square meter\";\n");
  fprintf(cdloutf, "G:valid_range = 0., 1400.;\n");
  fprintf(cdloutf, "G:lon_name = \" clear sky solar irradiance \";\n");
  fprintf(cdloutf, "B:units = \"Watt per square meter\";\n");
  fprintf(cdloutf, "B:valid_range = 0., 1200.;\n");
  fprintf(cdloutf, "B:lon_name = \" clear sky beam irradiance \";\n");
  fprintf(cdloutf,"data:\n");
  fprintf(cdloutf,"lat=\n"); 
  for(j=0;j<yrdim-1;j++)
    {fprintf(cdloutf,"%f ,\n",latrev[j]);}
  fprintf(cdloutf,"%f ;\n",latrev[yrdim-1]);
  fprintf(cdloutf,"lon=\n"); 
  for(i=0;i<xrdim-1;i++)
    {fprintf(cdloutf,"%f , \n",lonrev[i]);}
  fprintf(cdloutf,"%f ; \n",lonrev[xrdim-1]);
  fprintf(cdloutf,"G=\n");
  
  for(j=0;j<yrdim;j++){
    for(i=0;i<xrdim;i++){
      if (j==yrdim-1 && i==xrdim-1) {
	fprintf(cdloutf,"%f ;\n",var1[j][i]);
      }
      else
        fprintf(cdloutf,"%f ,\n",var1[j][i]);	    
    }
  }
  
    fprintf(cdloutf,"B=\n");
    for(j=0;j<yrdim;j++){
      for(i=0;i<xrdim;i++){
	if (j==yrdim-1 && i==xrdim-1) {
	  fprintf(cdloutf,"%f ;\n",var2[j][i]);
        }
        else
	  fprintf(cdloutf,"%f ,\n",var2[j][i]);	    
      }
    }
    fprintf(cdloutf,"\}"); 
    fclose(cdloutf); 

}

/*------------------------------------------------------------------------------------------------*/
