/* ------------------------------------------------------------------------------------------------------ 

Funktionen fuer Heliosat3: Cloud-Index
zur Berechnung der Geometrie aus Sonne, Satellit und Beobachtungspunkt
und zur Korrektur der VIS-Kanaele


Dr. Annette Hammer July 2005
Modifications and Add Ons Dr. R. Mueller May 2008
---------------------------------------------------------------------------------------------------------- */

/* RM May 2005: added pointer for satellite azimuth and zenith angle in the function call */
int geomsg(int NAV_CRES,int NAV_LRES,int NAV_COFF,int NAV_LOFF,int lin,int col,double *lat_ptr,double *lon_ptr)
{
/* Laenge und Breite aus HeaderInformationen nach VCS-Konvention PIF bzw XPIF
   Achtung: gezaehlt werden Line und Column ab 1 (Fortran-Konvention) 
 RM May 1008: added r_sat and r_e needed for calculation of satellite azimuth and zenith angle 
 r_sat=altitude of SAT from core of the earth 
 r_e=radius of the earth */

/* define h 42164.0
#define r_e 6378.169
#define r_p 6356.5838
#define RE2H2 -1737121000.0
#define REP2 1.006803  diese Angaben von VCS */
#define r_e 6378.169
#define r_sat 42168

double x,y;
double alpha,beta,r_n;
double cosalpha,cosbeta,sinalpha,sinbeta,hcosacosb,klammer;     /* Abkuerzungen, Zwischenergebnisse */
double r1,r2,r3,r_xy;
double discr;
int hrv;

if(NAV_CRES==222){NAV_CRES=222.623596;NAV_LRES=222.623596;hrv=3;}
 else if(NAV_CRES==667){NAV_CRES=667.2044067;NAV_LRES=667.2044067;hrv=1;}
 else{
    fprintf(stderr,"satgeomsgnew: Line and Column Resolution not valid !!\n");
    exit(1);
   }


/* Position im Fulldisk-Format (Hohe (HRV) und niedrige Aufloesung (z.B. vis006)) */
x = 1856 + (NAV_COFF - col)/hrv ;
y = 1856 + (NAV_LOFF - lin)/hrv ;

/* Weltall? bei Abstand von der Kreismitte > 1801 (Pseudoradius)*/
if(( (x-1856.0)*(x-1856.0) + (y-1856.0)*(y-1856.0) ) > 3243601.0)
  return(1);


alpha = (NAV_COFF-col)*17.832*PI/(hrv*3712*180);
beta = (NAV_LOFF-lin)*17.832*PI/(hrv*3712*180);

cosalpha=cos(alpha);
sinalpha=sin(alpha);
sinbeta=sin(beta);
cosbeta=cos(beta);
hcosacosb=h*cosalpha*cosbeta;
klammer=cosbeta*cosbeta + REP2*sinbeta*sinbeta;

discr = hcosacosb*hcosacosb + RE2H2*klammer;
if(discr < 0.0)
 {
   return(2);
 }

r_n = (hcosacosb - sqrt(discr))/klammer;

r1 = h - r_n*cosalpha*cosbeta;
r2 = -r_n*sinalpha*cosbeta;
r3 = r_n*sinbeta;

r_xy = sqrt(r1*r1 + r2*r2);


*lat_ptr = atan(REP2*r3/r_xy);
*lon_ptr = atan(r2/r1);

return(0);
}
/*------------------------------------------------------------------------------------------------*/
double equation_of_time(int doy, double *gamma, double *f, double *dec)
{
  /* Equation of Time in Minuten: 
     Abweichung der Sonne um 12:00 Mittags von der Suedposition
     siehe: Iqbal, An Introduction to Solar Radiation  */

  double E_t;
  double cosgamma,cos2gamma,singamma,sin2gamma;
  
  *gamma = 2*PI* (doy - 1.0)/365.0;                    /* Tageswinkel in Rad */
  cosgamma=cos( *gamma);
  cos2gamma=cos(2* *gamma);
  singamma=sin( *gamma);
  sin2gamma=sin(2* *gamma);
  
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

double cos_zenit(double lat,double lon,double dec,double lon_sun)
		/* Cosinus des Sonnenzenitwinkels */
{
double w;                                   /* Stundenwinkel */
double(cos_theta_zenith);

w = (lon -lon_sun);
cos_theta_zenith =sin(lat)*sin(dec)          /* siehe Iqbal oder  */
		 + cos(lat)*cos(dec)*cos(w); /* Juergens Dissertation  */
if(cos_theta_zenith<=0.0)
  {
   cos_theta_zenith=0.0;
   return(0);
  }
return(cos_theta_zenith);
}

/* RM May 2008: Added function for the calculation of the solar azimuth angle */

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




/*------------------------------------------------------------------------------------------------*/

double timeoffset(int navlres, int navloff, int lin, int nom_rep_cycle)
{
/* berechnet den Abscannzeitpunkt einer Zeile in einem Satellitenbild
   bezogen auf das Bildende (==Bildname) 
   Abscannzeitpunkt=Bildende+timeoffset in Stunden (dezimale Zeit)
   Achtung: Diese Version funktioniert nur fuer full disk und Bildausschnitte daraus,
            aber nicht fuer Rapid Scan !! */

#define RPM (100.0) /* rotations per minute of satellite */

double zeit;
double td_Aeq_end=17.5;  /* time distance: Minutes between Aequator-Scan and Image-Name */
                         /* for nominal repeat cycle = 30 min */

if((nom_rep_cycle!=30)&&(nom_rep_cycle!=15)) /* Meteosat 4-7 und MSG (Meteosat 8) */
     {
      fprintf(stderr,"Funktion time-offset: No valid Nominal Repeat Cycle given\n");
      exit(1);
     }

zeit=((navloff-lin)*navlres/1000.0)/RPM/60.0-td_Aeq_end/60;  /* for nominal repeat cycle = 30 min */
/*   (navloff-lin)navlres/1000: Abstand vom Aequator pos. oder neg. */

zeit=zeit*(nom_rep_cycle/30.0);  /* for nominal repeat cycle = 15 min or 30 min */
if(nom_rep_cycle==15) zeit=zeit+0.25; /* 0.25 difference between image name and end of image */
   
return(zeit);
}

/*------------------------------------------------------------------------------------------------*/

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

/*------------------------------------------------------------------------------------------------*/



