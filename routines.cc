//===============================================================================
//
//                      routines.cc
//  
// program to calculate:
//
// 1. LIGHT AT SEA SURFACE
// clear sky above-surface incident
// irradiance (all wavelengths, W m-2) and clear sky
// below-surface incident irradiance (PAR, uEin m-2 s-1)
// from sun zenith angle (See Lumb64, Katsaros90 and Kraus94) 
// as in Toby light.for program.
//
// 2. DAILY TEMPERATURE ROUTINE
// calculation of daily temperature variations
//
// 3. AVERAGED IRRADIANCE THROUGH DEPTH
// calculation of the irradiance averaged with depth 
// using the Three Layer Model described in: Anderson 1993 (ANDE93).
// Once the light profile is calculated than photosinthesys
// is determined either with a Michaelis-Menten's function or
// with the Steele's function. The latter includes saturation 
// and inhibition, see Totterdell 1993 (pag 330), or Kirk (pag 
// 274) for a reference.
//
// 4. LIGHT INTENSITY AT DEPTH
// calculation of light at a certain depth
//
// 5. MIN / MAX ROUTINES
// calculation of the minimum between two numbers
//
// 6. CARBONATE SYSTEM ROUTINES
// calculation of air-sea flux exchanges and other parameters in seawater
// water, all relavant equations and procedure are as in:
// - Peng et al. 1987 (PENG87, in Appendix, pag. 454)   
// - Wanninkhof 1992 (WNNI92, ia Appendix, pag. 7379)
// - Millero 1995 (MILL95)
// good explanations and applications of the procedures 
// can aslo be found in:
// - Bissett et al. 1999 (BISS99, pag. 259)
// - Tyrrell and Taylor 1995 (TYRR95, pag. 603)
// - Walsh and Dieterle 1994 (WALS94, pag. 358) - applications on the Bering Sea!
// - Taylor et al. 1991 (TAYL91, pag. 155)
// 
//
//                        Agostino Merico
// 
//             (Routines last modified: 20 Sep 2002)
//
//
//===============================================================================



#include <iostream.h>
#include <fstream.h>
#include <math.h>
 
#include "param.h"

#define SOLARC 1373.0 // solar constant in W m-2, see Kirk at pag 27
#define WTOE 4.17     // Watts (W m-2) to Einstein (uEin m-2 s-1) conversion 
                      // factor, see Lalli and Parsons at pag 261

#define N 60          // number of intervals into which the mixed layer is split

#define C 0.75        // cloudiness
#define T 0.75        // transmittance




//=============================== LIGHT AT SEA SURFACE ==================================


//double get_light_at_surface(double lat, double lon, double date, double gmt)

double get_light_at_surface(double n)  // this works when 'n' is any hour of the year
{


  double gmt;
  double date;

  double ltime, tau, psi, delta, theta, Ed1, Ed2, beta, rlat;

  // Bering Sea
  double lat=58.0;   // 58  N
  double lon=-165.0; // 165 W

  // North Atlantic - to test PENG87
  //double lat=64.0;   // 64 N
  //double lon=-28.0;  // 28 W

  // ===== GIVEN ANY HOUR (n) OF THE YEAR, CALCULATE DAY (date) AND TIME (gmt) =====

  modf(n/24, &date); // date = integer part of n/24;
  gmt=n-(24*date);   // time of the day = date 

  // ===============================================================================


  // conversion of GMT to local time
  ltime=gmt+(12.0*lon/180.0);
  if(ltime > 24.0){
    ltime=ltime-24.0;
    date=date+1.0;
  }
  else if(ltime < 0.0){
    ltime=ltime+24.0;
    date=date-1.0;
  }

  // conversion of latitude, date and time of day to radians
  rlat=(PI/180.0)*lat;
  psi=2.0*PI*(date/365.0);
  tau=2.0*PI*(ltime/24.0);

  // calculation of delta, solar declination (see Kirk, 2.2 pag 35)
  delta=0.39637-(22.9133*cos(psi))+(4.02543*sin(psi))-(0.3872*cos(2.0*psi))+(0.052*sin(2.0*psi));

  // conversion of delta to radians
  delta=(PI/180.0)*delta;
  
  // calculation of theta, the instantaneous solar zenith angle 
  // beta is the elevation angle
  beta=asin(sin(rlat)*sin(delta)-cos(rlat)*cos(delta)*cos(tau));
  theta=(PI/2.0)-beta;
  
  // calculation of above-surface downward irradiance (W m-2, at all wavelengths).
  // 24% loss is assumed due to vertical transmission through atmosphere, 
  // proportionately more as pathlength through atmosphere increases.
  // As in Lumb 64 model.

  // solar irradiance (for ignoring atmosphere multiply SOLARC by cos(theta) only) 
  Ed1=SOLARC*cos(theta);            // with atmosphere: (1.0-(0.24*cos(theta)));         
  if(Ed1<0.0) Ed1=0.0; //take only positive values

  // calculation of the same quantity for PAR (400-700nm) only, in uEin m-2 s-1 and just 
  // below the sea surface
  //  - Ed1 is the irradiance at the top of the atmosphere
  //  - 0.48 is the % of total irradiance at PAR wavelengths (at the bottom of
  //    the atmosphere, i.e. at the sea surface)
  //  - C is the cloudiness as a constant fraction
  //  - T is the transmittance through the sea surface 
  
  Ed2=Ed1*(1-0.7*C)*T*0.48;  // in W m-2 (multiply by WTOE to get units of uEin m-2 s-1) 
  // NOTE:
  // in Sverdrup '42 or '47 you'll find (1-0.071*C)
  // but in that case C goes from 1 to 10, here it goes from 0 to 1


  // *********** TESTING ************

  //Ed2=cos(theta);
  //if(Ed2<0.0) Ed2=0.0;

  // ********************************

  return Ed2;   // in W m-2

  // cout<<" calculation done...\n\n";
  // cout<<" GMT is "<<gmt<<"\n";
  // cout<<" local time is "<<ltime<<"\n";
  // cout<<" delta is "<<(delta*180.0/PI)<<"\n";
  // cout<<" beta is "<<(beta*180.0/PI)<<"\n";
  // cout<<" sun zenith angle is "<<(theta*180.0/PI)<<"\n";
  // cout<<" Clear sky above-surface all wavelengths is "<<Ed1<<" W m-2\n";
  // cout<<" Clear sky above-surface PAR (400-700nm) is "<<Ed2<<" uEin m-2 s-1\n\n";
	
}


//============================ DAILY TEMPERATURE ROUTINE ================================


double get_temperature(double n, double temp)
{
  double varT;
  double power;
  power=(temp-TMAX)*0.030103;   // 0.030103 = log10(2)/10
  varT=pow(10,power);

  return varT;
}


//======================== AVARAGED IRRADIANCE THROUGH MLD ============================

double get_light(double irr_surf, double chloro, double d)   // d is MLD
{

  int i=0;
  int j=0;

  int z=0;         // depth

  double Iz=0.0;   // light at a certain deoth
  double psi=0.0;  // averaged light in the MLD NOT in a M-M limiting form

  double res=0.0;
  double res1=0.0;
  double res2=0.0;

  int ires=0;
  int ires1=0;
  int ires2=0;

  double c;         // square root of pigment concentration
  
  double k[3]={0.0, 0.0, 0.0};  // attenuation coefficients
  
  double b[3][6] = {
    {0.13096, 0.030969, 0.042644, -0.013738, 0.0024617, -0.00018059},
    {0.041025, 0.036211, 0.062297, -0.030098, 0.0062597, -0.00051944},
    {0.021517, 0.050150, 0.058900, -0.040539, 0.0087586, -0.00049476},
  };

 // --- Three layer model (Anderson's) ---

  // irr_surf is in W m-2

  c=sqrt(chloro); // chloro is = G in Tom's model. sqrt is to ensure a bias 
                  // toward smaller values (more frequent in nature)
 
  // calculate att. coeff. for the three layers (k1, k2, k3) 
  for(i=0;i<3;i++){
    for(j=0;j<6;j++) k[i]+=b[i][j]*pow(c,j);      
  }

  //case 1: M<=5
  if(d<=5){

    // calculate averaged light in layer 1 (0-5 m)
    for(z=1;z<=30;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*d/30.0));
      psi+=Iz;
    } 
  }
  //case 2: 5<M<=23
  if((d>5) && (d<=23)){
    res = (30.0/d)*5.0;  // number of intervals in the first layer
    ires = (int) res;    // take the integer part of res

    // calculate averaged light in layer 1 (0-5 m)
    for(z=1;z<=ires;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5.0/30.0));
      psi+=Iz; 
    } 
    // calculate averaged light in layer 2 (5-23 m)
    for(z=ires+1;z<=30;z++){
      Iz=irr_surf*exp(-5.0*k[0])*exp(-k[1]*((z-0.5)*d/30.0));
      psi+=Iz;
    }
  }
  //case 3: M>23
  if(d>23){
    res1 = (30.0/d)*5.0;    // number of intervals in the first layer
    res2 = (30.0/d)*18.0;   // number of intervals in the second layer
    ires1 = (int) res1;     // take the integer part of res1
    ires2 = (int) res2;     // take the integer part of res2

    // calculate averaged light in layer 1 (0-5 m)
    for(z=1;z<=ires1;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5.0/30.0));
      psi+=Iz;
    } 
    // calculate averaged light in layer 2 (5-23 m)
    for(z=ires1+1;z<=ires2;z++){
      Iz=irr_surf*exp(-5.0*k[0])*exp(-k[1]*((z-0.5)*23.0/30.0));
      psi+=Iz;
    }
    // calculate averaged light in layer 3 (23-60 m)
    for(z=ires2+1;z<=30;z++){
      Iz=irr_surf*exp(-5.0*k[0])*exp(-18.0*k[1])*exp(-k[2]*((z-0.5)*d/30.0));
      psi+=Iz;
    }
  }
  psi=psi/30.0;


  // ==========

  return psi;

  // ==========
}

//======================== LIMITING IRRADIANCE THROUGH DEPTH ============================


double get_averaged_light(double irr_surf, double chloro, double d)   // d is MLD
{

  int i=0;
  int j=0;

  int z=0;            // depth
  int nz=0;           // new depth

  double Iz2=0.0;     // irradiance calculated at depth z with 2w
  double Iz1=0.0;     // irradiance calculated at depth z with 1w
  double psi2w=0.0;   // M-M light limitation term (two wave-band model) 
  double psi1w=0.0;   // M-M light limitation term (single wave-band model) 
  double psi3l=0.0;   // M-M light limitation term (three-layer model) 
  double psi3lsat=0.0;// Sat light limitation term (three-layer model)
  double psi2wsat=0.0;// Sat light limitation term (two-waveband model)

  double psia=0.0;    // averaged light in the MLD

  double res=0.0;
  double res1=0.0;
  double res2=0.0;

  int ires=0;
  int ires1=0;
  int ires2=0;

  double Iz=0.0;    // irradiance calculated at depth z with 3l
  double integ=0.0; // irradiance integrated through depth (averaged light)
  double c;         // square root of pigment concentration
  
  double k[3]={0.0, 0.0, 0.0};  // attenuation coefficients
  
  double b[3][6] = {
    {0.13096, 0.030969, 0.042644, -0.013738, 0.0024617, -0.00018059},
    {0.041025, 0.036211, 0.062297, -0.030098, 0.0062597, -0.00051944},
    {0.021517, 0.050150, 0.058900, -0.040539, 0.0087586, -0.00049476},
  };

 // --- Three layer model (Anderson's) ---

  // irr_surf is in W m-2

  c=sqrt(chloro); // chloro is = G in Tom's model. sqrt is to ensure a bias 
                  // toward smaller values (more frequent in nature)
 
  // calculate att. coeff. for the three layers (k1, k2, k3) 
  for(i=0;i<3;i++){
    for(j=0;j<6;j++) k[i]+=b[i][j]*pow(c,j);      
  }
  //case 1: M<=5
  if(d<=5){
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=30;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*d/30.0));
      psi3l+=Iz/(Iz+IHD);                     // Michaelis-Menten's function (see TOTT93a pag 330 case ii)
      psi3lsat+=(Iz/ISAT)*exp(1.0-Iz/ISAT);   // Steele's function (see TOTT93a pag 330 case iii)
      psia+=Iz;
    } 
  }
  //case 2: 5<M<=23
  if((d>5) && (d<=23)){
    res = (30.0/d)*5.0;  // number of intervals in the first layer
    ires = (int) res;    // take the integer part of res
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=ires;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5.0/30.0));
      psi3l+=Iz/(Iz+IHD);                     // Michaelis-Menten's function (see TOTT93a pag 330 case ii)
      psi3lsat+=(Iz/ISAT)*exp(1.0-Iz/ISAT);   // Steele's function (see TOTT93a pag 330 case iii)
    } 
    // calculate light attenuation in layer 2 (5-23 m)
    for(z=ires+1;z<=30;z++){
      Iz=irr_surf*exp(-5.0*k[0])*exp(-k[1]*((z-0.5)*d/30.0));
      psi3l+=Iz/(Iz+IHD);                     // Michaelis-Menten's function (see TOTT93a pag 330 case ii)
      psi3lsat+=(Iz/ISAT)*exp(1.0-Iz/ISAT);   // Steele's function (see TOTT93a pag 330 case iii)
    }
  }
  //case 3: M>23
  if(d>23){
    res1 = (30.0/d)*5.0;    // number of intervals in the first layer
    res2 = (30.0/d)*18.0;   // number of intervals in the second layer
    ires1 = (int) res1;     // take the integer part of res1
    ires2 = (int) res2;     // take the integer part of res2
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=ires1;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5.0/30.0));
      psi3l+=Iz/(Iz+IHD);                     // Michaelis-Menten's function (see TOTT93a pag 330 case ii)
      psi3lsat+=(Iz/ISAT)*exp(1.0-Iz/ISAT);   // Steele's function (see TOTT93a pag 330 case iii)
    } 
    // calculate light attenuation in layer 2 (5-23 m)
    for(z=ires1+1;z<=ires2;z++){
      Iz=irr_surf*exp(-5.0*k[0])*exp(-k[1]*((z-0.5)*23.0/30.0));
      psi3l+=Iz/(Iz+IHD);                     // Michaelis-Menten's function (see TOTT93a pag 330 case ii)
      psi3lsat+=(Iz/ISAT)*exp(1.0-Iz/ISAT);   // Steele's function (see TOTT93a pag 330 case iii)
    }
    // calculate light attenuation in layer 3 (23-60 m)
    for(z=ires2+1;z<=30;z++){
      Iz=irr_surf*exp(-5.0*k[0])*exp(-18.0*k[1])*exp(-k[2]*((z-0.5)*d/30.0));
      psi3l+=Iz/(Iz+IHD);                     // Michaelis-Menten's function (see TOTT93a pag 330 case ii)
      psi3lsat+=(Iz/ISAT)*exp(1.0-Iz/ISAT);   // Steele's function (see TOTT93a pag 330 case iii)
    }
  }
  psi3l=psi3l/30.0;
  psi3lsat=psi3lsat/30.0;


  // --- Two waveband approximation ---

  for(z=1;z<=30;z++){     //1.9875 makes Chl in mmol/m3
    Iz2=0.5*irr_surf*(exp(-(KGR + KSS*chloro/1.9875)*((z-0.5)*(d/30.0))) + 
                      exp(-(KRE + KSS*chloro/1.9875)*((z-0.5)*(d/30.0)))); 
    psi2w+=Iz2/(Iz2+IHD);
    psi2wsat+=(Iz2/ISAT)*exp(1-(Iz2/ISAT));

  }
  psi2w=psi2w/30.0;
  psi2wsat=psi2wsat/30.0;

  // --- Single waveband approximation ---

  for(z=1;z<=30;z++){
    Iz1=irr_surf*exp(-(KW + KSS*chloro/1.9875)*((z-0.5)*(d/30.0)));  //1.9875 makes Chl in mmol/m3 
    psi1w+=Iz1/(Iz1+IHD);
  }
  psi1w=psi1w/30.0;


  // ===== Return - to all but Ehux =====

  //return psi3l;       // Michaelis-Menten's fashion
  return psi3lsat;    // Steele's fashion

  //return psi2w;
  //return psi2wsat;

  // ====================================

}


double get_averaged_light_eh(double irr_surf, double chloro, double d) //with IHEH required by ehux
{

  int i=0;
  int j=0;

  int z=0;          // depth
  int nz=0;         // new depth

  double Iz2=0.0;   // irradiance calculated at depth z with 2w
  double Iz1=0.0;   // irradiance calculated at depth z with 1w
  double psi2w=0.0;
  double psi1w=0.0;
  double psi3l=0.0;

  double psi3lsat=0.0;
  double psi2wsat=0.0;

  double res=0.0;
  double res1=0.0;
  double res2=0.0;

  int ires=0;
  int ires1=0;
  int ires2=0;

  double Iz=0.0;    // irradiance calculated at depth z with 3l
  double integ=0.0; // irradiance integrated through depth (averaged light)
  double c;         // square root of pigment concentration
  
  double k[3]={0.0, 0.0, 0.0};  // attenuation coefficients
  
  double b[3][6] = {
    {0.13096, 0.030969, 0.042644, -0.013738, 0.0024617, -0.00018059},
    {0.041025, 0.036211, 0.062297, -0.030098, 0.0062597, -0.00051944},
    {0.021517, 0.050150, 0.058900, -0.040539, 0.0087586, -0.00049476},
  };

 // --- Three layer model (Anderson's) ---

  c=sqrt(chloro); // chloro is = G in Tom's model. sqrt is to ensure a bias 
                  // toward smaller values (more frequent in nature)
 
  // calculate att. coeff. for the three layers (k1, k2, k3) 
  for(i=0;i<3;i++){
    for(j=0;j<6;j++) k[i]+=b[i][j]*pow(c,j);      
  }
  //case 1: M<=5
  if(d<=5){
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=30;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*d/30));
      psi3l+=Iz/(Iz+IHEH);
      psi3lsat+=(Iz/ISATEH)*exp(1.0-Iz/ISATEH);   // Steele's function (see TOTT93a pag 330 case iii)
    } 
  }
  //case 2: 5<M<=23
  if((d>5) && (d<=23)){
    res = (30/d)*5;  // number of intervals in the first layer
    ires = (int) res; // take the integer part of res
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=ires;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5/30));
      psi3l+=Iz/(Iz+IHEH);
      psi3lsat+=(Iz/ISATEH)*exp(1.0-Iz/ISATEH);   // Steele's function (see TOTT93a pag 330 case iii)
    } 
    // calculate light attenuation in layer 2 (5-23 m)
    for(z=ires+1;z<=30;z++){
      Iz=irr_surf*exp(-5*k[0])*exp(-k[1]*((z-0.5)*d/30));
      psi3l+=Iz/(Iz+IHEH);
      psi3lsat+=(Iz/ISATEH)*exp(1.0-Iz/ISATEH);   // Steele's function (see TOTT93a pag 330 case iii)
    }
  }
  //case 3: M>23
  if(d>23){
    res1 = (30/d)*5;    // number of intervals in the first layer
    res2 = (30/d)*18;   // number of intervals in the second layer
    ires1 = (int) res1; // take the integer part of res1
    ires2 = (int) res2; // take the integer part of res2
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=ires1;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5/30));
      psi3l+=Iz/(Iz+IHEH);
      psi3lsat+=(Iz/ISATEH)*exp(1.0-Iz/ISATEH);   // Steele's function (see TOTT93a pag 330 case iii)
    } 
    // calculate light attenuation in layer 2 (5-23 m)
    for(z=ires1+1;z<=ires2;z++){
      Iz=irr_surf*exp(-5*k[0])*exp(-k[1]*((z-0.5)*23/30));
      psi3l+=Iz/(Iz+IHEH);     
      psi3lsat+=(Iz/ISATEH)*exp(1.0-Iz/ISATEH);   // Steele's function (see TOTT93a pag 330 case iii)
    }
    // calculate light attenuation in layer 3 (23-60 m)
    for(z=ires2+1;z<=30;z++){
      Iz=irr_surf*exp(-5*k[0])*exp(-18*k[1])*exp(-k[2]*((z-0.5)*d/30));
      psi3l+=Iz/(Iz+IHEH);      
      psi3lsat+=(Iz/ISATEH)*exp(1.0-Iz/ISATEH);   // Steele's function (see TOTT93a pag 330 case iii)
    }
  }
  psi3l=psi3l/30.0;
  psi3lsat=psi3lsat/30.0;


  // --- Two waveband approximation ---
  
  for(z=1;z<=30;z++){
    Iz2=0.5*irr_surf*(exp(-(KGR + KSS*chloro/1.9875)*((z-0.5)*(d/30.0))) +  //1.9875 makes Chl in mmol/m3
		      exp(-(KRE + KSS*chloro/1.9875)*((z-0.5)*(d/30.0)))); 
    psi2w+=Iz2/(Iz2+IHEH);
    psi2wsat+=(Iz2/ISATEH)*exp(1-(Iz2/ISATEH));
  }
  psi2w=psi2w/30.0;
  psi2wsat=psi2wsat/30.0;
  
  
  // --- Single waveband approximation ---
  
  for(z=1;z<=30;z++){
    Iz1=irr_surf*exp(-(KW + KSS*chloro/1.9875)*((z-0.5)*(d/30.0)));  //1.9875 makes Chl in mmol/m3 
    psi1w+=Iz1/(Iz1+IHEH);
  }
  psi1w=psi1w/30.0;
  
 
  // ==== Return - only to Ehux ====

  //return psi3l;     // Michaelis-Menten's fashion
  return psi3lsat;  // Steele's fashion

  //return psi2w;
  //return psi2wsat;

  // ==============================
}


double get_averaged_light_cal(double irr_surf, double chloro, double d) //with IHEH required by ehux
{

  int i=0;
  int j=0;

  int z=0;          // depth
  int nz=0;         // new depth

  double Iz2=0.0;   // irradiance calculated at depth z with 2w
  double Iz1=0.0;   // irradiance calculated at depth z with 1w
  double psi2w=0.0;
  double psi1w=0.0;
  double psi3l=0.0;

  double psi2wsat=0.0;

  double res=0.0;
  double res1=0.0;
  double res2=0.0;

  int ires=0;
  int ires1=0;
  int ires2=0;

  double Iz=0.0;    // irradiance calculated at depth z with 3l
  double integ=0.0; // irradiance integrated through depth (averaged light)
  double c;         // square root of pigment concentration
  
  double k[3]={0.0, 0.0, 0.0};  // attenuation coefficients
  
  double b[3][6] = {
    {0.13096, 0.030969, 0.042644, -0.013738, 0.0024617, -0.00018059},
    {0.041025, 0.036211, 0.062297, -0.030098, 0.0062597, -0.00051944},
    {0.021517, 0.050150, 0.058900, -0.040539, 0.0087586, -0.00049476},
  };

 // --- Three layer model (Anderson's) ---

  c=sqrt(chloro); // chloro is = G in Tom's model. sqrt is to ensure a bias 
                  // toward smaller values (more frequent in nature)
 
  // calculate att. coeff. for the three layers (k1, k2, k3) 
  for(i=0;i<3;i++){
    for(j=0;j<6;j++) k[i]+=b[i][j]*pow(c,j);      
  }
  //case 1: M<=5
  if(d<=5){
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=30;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*d/30));
      psi3l+=Iz/(Iz+IHCA);
    } 
  }
  //case 2: 5<M<=23
  if((d>5) && (d<=23)){
    res = (30/d)*5;  // number of intervals in the first layer
    ires = (int) res; // take the integer part of res
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=ires;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5/30));
      psi3l+=Iz/(Iz+IHCA);
    } 
    // calculate light attenuation in layer 2 (5-23 m)
    for(z=ires+1;z<=30;z++){
      Iz=irr_surf*exp(-5*k[0])*exp(-k[1]*((z-0.5)*d/30));
      psi3l+=Iz/(Iz+IHCA);
    }
  }
  //case 3: M>23
  if(d>23){
    res1 = (30/d)*5;    // number of intervals in the first layer
    res2 = (30/d)*18;   // number of intervals in the second layer
    ires1 = (int) res1; // take the integer part of res1
    ires2 = (int) res2; // take the integer part of res2
    // calculate light attenuation in layer 1 (0-5 m)
    for(z=1;z<=ires1;z++){
      Iz=irr_surf*exp(-k[0]*((z-0.5)*5/30));
      psi3l+=Iz/(Iz+IHCA);
    } 
    // calculate light attenuation in layer 2 (5-23 m)
    for(z=ires1+1;z<=ires2;z++){
      Iz=irr_surf*exp(-5*k[0])*exp(-k[1]*((z-0.5)*23/30));
      psi3l+=Iz/(Iz+IHCA);     
    }
    // calculate light attenuation in layer 3 (23-60 m)
    for(z=ires2+1;z<=30;z++){
      Iz=irr_surf*exp(-5*k[0])*exp(-18*k[1])*exp(-k[2]*((z-0.5)*d/30));
      psi3l+=Iz/(Iz+IHCA);      
    }
  }
  psi3l=psi3l/30.0;

  // --- Two waveband approximation ---
  
  for(z=1;z<=30;z++){
    Iz2=0.5*irr_surf*(exp(-(KGR + KSS*chloro/1.9875)*((z-0.5)*(d/30.0))) +  //1.9875 makes Chl in mmol/m3
		      exp(-(KRE + KSS*chloro/1.9875)*((z-0.5)*(d/30.0)))); 
    psi2w+=Iz2/(Iz2+IHCA);
  }
  psi2w=psi2w/30.0;
  
  
  // --- Single waveband approximation ---
  
  for(z=1;z<=30;z++){
    Iz1=irr_surf*exp(-(KW + KSS*chloro/1.9875)*((z-0.5)*(d/30.0)));  //1.9875 makes Chl in mmol/m3 
    psi1w+=Iz1/(Iz1+IHCA);
  }
  psi1w=psi1w/30.0;
  
 
  // ==== Return - only to calcification ====

  return psi3l;     // Michaelis-Menten's fashion



  // ==============================
}



//=================== CALCULATE LIGHT INTENSITY AT A GIVEN DEPTH =======================

double get_light_intensity(double s_irrad, double chl){

  int i=0;
  int j=0;
  
  double f=0.0;         // irradiance calculate at depth z

  double c=0.0;         // square root of pigment concentration
  
  double k[3]={0.0, 0.0, 0.0};  // attenuation coefficients
  
  double b[3][6] = {
    {0.13096, 0.030969, 0.042644, -0.013738, 0.0024617, -0.00018059},
    {0.041025, 0.036211, 0.062297, -0.030098, 0.0062597, -0.00051944},
    {0.021517, 0.050150, 0.058900, -0.040539, 0.0087586, -0.00049476},
  };

  c=sqrt(chl); // chloro is = to G in Tom's model. sqrt is to ensure a bias 
               // toward smaller values (more frequent in nature)
 
  // calculate att. coeff. for the three layers (k1, k2, k3) 
  for(i=0;i<3;i++){
    for(j=0;j<6;j++) k[i]+=b[i][j]*pow(c,j);      
  }

  // calculate light attenuation in layer 1 (0-5 m) at 5 m
  f=s_irrad*exp(-k[0]*5);


  return f;
}
  
  

//================================= MINIMUM ROUTINE ====================================


double min(double x, double y)

{

  double minimum=0.0;
  
  if(x<=y) minimum=x;
  else minimum=y;
 
  return minimum;
  
}

//================================= MAXIMUM ROUTINE ====================================


double max(double x, double y)

{

  double maximum=0.0;
  
  if(x>=y) maximum=x;
  else maximum=y;
 
  return maximum;
  
}



//============================ CARBONATE SYSTEM ROUTINES ===============================


// Computation of gas transfer velocity using:
// Equation 8 and table A1 in Wanninkhof 1992 (WANNI92)

// ws    must be in: m s-1
// gastv will be in: cm h-1
double get_gas_transfer_velocity(double ws, double te)       
{
  double gastv=0.0; // gas transfer velocity in: cm hr-1
  double sc=0.0;    // Schmidt number

  sc = 2073.1 - 125.62*te + 3.6276*te*te - 0.043219*te*te*te;

  //gastv=0.31*ws*ws*pow((sc/660.0),-0.5);
  gastv = (2.5*(0.5246 + 1.6256*pow(10.0,-2)*te + 4.9946*pow(10.0,-4)*te*te) + 0.3*ws*ws)*pow((sc/660.0),-0.5);

  gastv=gastv/100.0; // transform gstv in: m hr-1

  return gastv; 
}


// Computation of CO2 solubility in Seawater using:
// Equation and values in table A2 in Wanninkhof 1992 (WANNI92)
// Method originally described by Weiss 1974
//
//
double get_co2_solubility(double sa, double te)          
{
  double csw=0.0; // CO2 solubility in seawater

  double tek=te+273.15;   // temperature in K

  // calculate CO2 solubility in mol kg-1 atm-1
  csw=exp(-60.2409 + 9345.17/tek + 23.3585*log(tek/100.0) + 
           sa*(0.023517 - 0.023656*(tek/100.0) + 0.0047036*(tek/100.0)*(tek/100.0)));

  // convert units into umol kg-1 atm-1
  csw=csw*1.0e6;

  return csw;
}


// Computation of a certain parameter in Seawater using methods described by:
//  - Peng et al. 1987 (PENG87) in Appendix at pag. 454
//  - Millero 1995 (MILL95)
//  - Archer at: "http://geosci.uchicago.edu/~archer/cgimodels/"
//
// Or look at the routine: '/users/phd/am1/model/carbo/carbconst.r'
// Or look at the routine: '/users/phd/am1/model/14c/carbof/pco2.f'
//
//
// Output flags:
//    r=1 -> get pCO2 in seawater
//    r=2 -> get [CO3]
//    r=3 -> get omega calcite
//    r=4 -> get omega aragonite
//    r=5 -> get pH
//    r=6 -> get [HCO3-]
//    r=7 -> get [Co2(aq)]
//
//                     salinity  temperature alkalinity    TCO2    silicate  output flag
double get_param_water(double sa, double te, double al, double co, double si, int r)
{

  double tek=0.0;

  double pco2=0.0;  // output of this function: pCO2

  double kc1, kc2, kb, kp1, kp2, kp3, kw; 
  kc1=kc2=kb=kp1=kp2=kp3=kw=0.0;
  
  double bo=0.0;  

  double fh=0.0;

  double ah1=0.0;

  double ac,ab,as,ap,aw;
  ab=as=ap=aw=0.0;
  
  double talk=0.0;

  double cp=0.0;
  double prat=0.0;
  double z=3.0;   // depth (to take into account effects of pressure on solubility)

  // get the proper units for calculations 
  al=al*1.0e-6;   // Eq m-3
  si=si*1.0e-6;   // mol m-3
  double po=0.0;  // mol m-3
  co=co*1.0e-6;   // mol m-3

  tek=te+273.15;        // temperature in K
  double pres = z/10.0; // pressure in bar

  // computing - dissociation constants for carbonic acid (Mehrbach et al., 1973)
  cp = (pres)/83.143/tek;    // P/RT, where R = 83.143 cm3 bar mol-1 K-1, ideal gas constant

  kc1 = 13.7201 - 0.031334*tek - 3235.76/tek - 1.3e-5*sa*tek + 0.1032*pow(sa,0.5);
  kc1 = pow(10,kc1);
  prat = (24.2 - 0.085*te)*cp;
  prat = exp(prat);
  kc1 = kc1*prat;

  kc2 = -5371.9645 - 1.671221*tek + 128375.28/tek + 2194.3055*log(tek)/2.30259 - 0.22913*sa - 
         18.3802*log(sa)/2.30259 + 8.0944e-4*sa*tek + 5617.11*log(sa)/tek/2.30259 - 2.136*sa/tek;
  kc2 = pow(10,kc2);
  prat = (16.4 - 0.04*te)*cp;
  prat = exp(prat);
  kc2 = kc2*prat;

  // computing - dissociation constants for boric acid (Lyman)
  kb = 2291.9/tek + 0.01756*tek - 3.385 - 0.32051*pow((sa/1.80655),(1.0/3.0));
  kb = pow(10,-kb);
  prat = (27.5 - 0.095*te)*cp;
  prat = exp(prat);
  kb = kb*prat;

  // computing - dissociation constants for silicic acid (PENG87)
  static double ks;
  ks = 1.0e-10;

  // computing - dissociation constants for phosphoric acid, kp2 and kp3 only
  kp2 = exp(-9.039 - 1450.0/tek);
  kp3 = exp(4.466 - 7276.0/tek);

  // computing - dissociation constants for water in seawater (PENG87, source: Millero 1979)
  kw  = exp(148.9802 - 13847.26/tek - 23.6521*log(tek) - 0.019813*sa + 
        pow(sa,0.5)*(-79.2447 + 3298.72/tek + 12.0408*log(tek)));

  // computing - total activity coefficient for hydrogen ion (PENG87)
  fh = 1.29 - 0.00204*tek + 4.6*1.0e-4*sa*sa - 1.48*1.0e-6*sa*sa*tek;

  // computing - total borate concentration (PENG87)
  bo = 4.106e-4*(sa/35);

  double c1 = kc1/2.0;
  double c2 = 1.0-4.0*kc2/kc1;
  double c4 = bo*kb;

  double aht = 1.0e-8;

  int icnt=1;

  // computing the several alkalinity species iteratively 
  // starting with an initial 'aht' trial value
  while(0.5e-4 < fabs(1.0 - (aht/ah1)) || icnt <= 100){

    ab = bo*kb/(aht + kb);                          // borate alkalinity
    as = si*4*1.0e-10/(aht + 4*1.0e-10);            // silicate alkalinity
    ap = po*(1/(1 + kp2/aht + kp2*kp3/(aht*aht)) +  // phosphate alkalinity 
	     2/(1 + aht/kp2 + kp3/aht) + 
	     3/(1 + aht/kp3 + aht*aht/kp2*kp3));
    aw = (kw*fh/aht) - (aht/fh);                    // water alkalinity

    ac = al - ab - as - ap - aw;                    // carbonate alkalinity  
    
    double X=ac/co;
    ah1= c1/X*(1.0 - X + sqrt(1.0 + c2*X*(-2.0 + X)));
    aht=ah1;

    icnt++;

  }

  double co3 = (ac - co)/(1.0 - (ah1*ah1)/(kc1*kc2));
  double hco3 = co/(1.0 + ah1/kc1 + kc2/ah1);
  double co2 = co/(1.0 + kc1/ah1 + kc1*kc2/(ah1*ah1));

  double Is;
  double khco2; // co2 solubility in seawater

  // saturation of CO2, from Edmond and Geiskes (Stumm and Morgan p 204)
  Is = 0.00147 + 0.03592*sa/2.0 + .000068*(sa/2.0)*(sa/2.0);
  khco2 = -2385.73/tek + 14.0184 - 0.0152642*tek + Is*(0.28596 - 6.167e-4*tek);
  
  khco2 = pow(10,-khco2);    // Kh = { HCO3* } / pCO2

  // calculate solubility of CO2 in seawater:  kco2  using the formulation of
  // Weiss (1974, Marine Chem., 2, 203-215) in mol kg-1 atm-1
  khco2 = exp(-60.2409 + 9345.17/tek + 23.3585*log(tek/100.0) + 
              sa*(0.023517 - 2.3656e-4*tek + 4.7036e-7*tek*tek));

  
  double atm = PCO2A/1.0e6;         // pCO2(air) in atm
  double satco2 = khco2*atm*1000.0; // in mol m-3

  // calculate Steady State pCO2
  pco2 = co2/khco2;                 // in atm

  //calculate h2co3
  double h2co3 = khco2*pco2;

  // calcite solubility, from Sayles 1980
  double kprime = 4.75e-7;  // mol2 kg-2
  double delv = -44.0;      // cm3 mol-1
  double rr = 83.143;       // cm3 bar K-1 mol-1 (ideal gas const.: R = 8.314 J mol-1 K-1)
  double dk = -0.0133;      // cm3 bar-1 mol-1
  double kpres = log(10.0);
  kpres = log(kprime) - delv/(rr*tek)*(pres) + 0.5*dk/(rr*tek)*pres*pres;
  kpres = exp(kpres);
  double csat = kpres/0.01;  

  double pH;
  double hplus = kc2*hco3/co3;
  pH = -log(hplus)/2.303;

  // calculate [ca++] from Millero pag. 270, eq. 127
  double calcium = 0.01028*(sa/35.0); 

  // calculate omega-calcite and omega-aragonite (from Mucci 1983, see also Millero pag. 249)
  double lnksp0cal = -395.8293 + 6537.773/tek + 71.595*log(tek) - 0.17959*tek;
  double lnksp0arag = -395.9180 + 6685.079/tek + 71.595*log(tek) - 0.17959*tek;
 
  double kcal = exp(lnksp0cal + (-1.78938 + 410.64/tek + 0.0065453*tek)*sqrt(sa) -
		           0.17755*sa + 0.0094979*pow(sa,1.5));

  double karag = exp(lnksp0arag + (-0.157481 + 202.938/tek + 0.0039780*tek)*sqrt(sa) -
	                     0.23067*sa + 0.0136808*pow(sa,1.5));

  double saltrat = sqrt(sa)/sqrt(35.0);

  // calcite volume correction
  double dvc = -45.464 + 0.3529*te - 4.985*te*te*pow(10,-3.0)*saltrat;        
  // calcite compressibility correction
  double dkc = (-13.70 + 0.1245*te + 0.0*te*te*pow(10,-3.0))*pow(10,-3.0)*saltrat;
  // aragonite volume correction
  double dva = -42.680 + 0.3529*te - 4.985*te*te*pow(10,-3.0)*saltrat;
  // aragonite compressibility correction
  double dka = (-13.70 + 0.1245*te + 0.0*te*te*pow(10,-3.0))*pow(10,-3.0)*saltrat;
  
  // Millero's book, pag 249
  //double dvc = -48.76 + 0.5302*te;
  //double dkc = (-11.76 + 0.3692*te)/1000;
  //double dva = -46.0 + 0.5304*te;
  //double dka = (-11.76 + 0.3692*te)/1000;

  kcal = kcal*exp(-dvc*cp + 0.5*dkc*(pres*pres)/rr/tek);
  karag = karag*exp(-dva*cp + 0.5*dka*(pres*pres)/rr/tek);

  double omega_cal = (calcium*co3)/kcal;
  double omega_arag = (calcium*co3)/karag;

  if(r==1) return pco2;         // return pCO2 in atm
  if(r==2) return co3*1.0e6;    // return [CO3=] in umol kg-1
  if(r==3) return omega_cal;    // return omega calcite
  if(r==4) return omega_arag;   // return omega aragonite
  if(r==5) return pH;           // return pH
  if(r==6) return hco3*1.0e6;   // return [HCO3-] in umol kg-1
  if(r==7) return co2*1.0e6;    // return [CO2(aq)] in umol kg-1
}

