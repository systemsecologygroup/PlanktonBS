//
//
//                           succession4new.cc
//
//                 3P-2Z-N-S-D-LA-LF-DIC-Alk version
//
// Program to study four phytoplankton species succession 
// in the Bering Sea with two zooplankton (micro, meso), three nutrients (nitrate, 
// ammonium, silicate), detritus, attached and free coccoliths, DIC and total alkalinity. 
// Carbonate concentration, Omega-calcite and Omega-aragonite cycles are also calculated. 
//
// NOTE: In the pre-1995 run, Ehux and coccoliths are set to zero so to avoid
//       the carbon system to be affected by the Ehux presence before 1995. 
//
//
//                         Agostino Merico
// 
//                     (Last modified Sep 2003) 
//  
//
//  (Previous modification Mar 2003 to make a multi-year run 
//  before 1995 followed by a transient run from 1995 to 2002)
//
//
//
//  to be compiled using g++ and together with:
//  
//                routines.cc  
//                nrutil.cc 
//
//  EXAMPLE: 
//  to compile type:  g++ succession4new.cc routines.cc nrutil.cc -o a.out -Wno-deprecated
//  to run type:      ./a.out 
//
//
//  INPUT FILES:  mldXX.in (or mldnew.dat for 'Fasham-modified' MLD) 
//                sstXX.in 
//                winXX.in
//                sal.in
// 
//  HEADER FILES: param.h
//                nrutil.h
//
//
//  CRUCIAL PARAMETERS: Y    (number of years to run the model)
//                      IGNY (number of initial years to ignore for steady-state)
//                      HOFY (hour of the year to consider for poincare' sections)  
//
//                      trans (set TRUE to look at transient results)
//                            (set FALSE to look at steady-state results)
//
//
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "param.h"     // parameters and prototype functions
#include "nrutil.h"    // required by function rk4

#define TRUE 1         // first year run
#define FALSE 0        // after first year run

#define NEQ 14         // number of ODEs
#define STEP 366       // number of steps in integration

#define Y 9            // number of years for which run the model (0 is one year cycle)

#define IGNY 0         // number of initial years to ignore (for steady-state)

#define HOFY 4320      // hour of the year to consider for poincare' sections

#define DSTEP 365      // number of steps [number of days in one year]
#define HSTEP 8760     // number of steps [number of hours in one year]
#define HHSTEP 17520   // number of steps [number of half hours in one year]

#define HYSTEP 61320   // number of steps [number of hours in seven years 1995-2001]

#define TIME 8758      // normalisation factor for time unit (24 gives X-axis in julian day)
                       //                                    (8760 gives X-axis in years)
                  
#define TI 1           // initial time [day, huor]
#define TD 366         // final time [day]
#define TH 8761        // final time [hour]
#define THH 17521      // final time [1/2 hour]


// ======== FUNCTIONS ======== 

void rkdriver(double vstart[], int nvar, double t1, double t2, int nstep, 
	      void (*derivs)(double, double [], double []));

void rk4(double y[], double dydt[], int n, double t, double h, double yout[], 
	 void (*derivs)(double, double [], double []));

void derivs(double t, double y[], double dydt[]);


// ===== GLOBAL VARIABLES =====
           
static int trans = TRUE; // set to TRUE  to look at transient results
                         // set to FALSE to look at steady-state

int yy;

//double MEH=0.0;

double ingEH=0.0;      // variable ingestion rate for Ehuxleyi (depends on omega)
double ingDI=0.0;      // variable ingesiton rate for diatoms (depends on silicate)

double pon=0.0;          // particulate organic nitrogen  
double photoeh=0.0;      // photosinthetic rate of Ehux
double calcieh=0.0;      // calcification rate of Ehuc
double newphypro=0.0;    // new primary production

double regphypro=0.0;    // regenerated primary production
double regdiapro=0.0;    // regenerated diatoms primary production
double regdinpro=0.0;    // regenerated dinofla primary production
double regflapro=0.0;    // regenerated flagell primary production
double regehupro=0.0;    // regenerated ehuxley primary production
double regtest=0.0;

double totphypro=0.0;    // primary production as given by the phytoplankton growth term
double totzoopro=0.0;    // total zooplankton biomass
double totphyloss=0.0;   // total phytoplankton loss
double totzooloss=0.0;   // total zooplankton loss
double totphymix=0.0;    // total phytoplankton mixing

double dianutgro=0.0;
double dinnutgro=0.0;
double flanutgro=0.0;
double ehunutgro=0.0;

double dialightgro=0.0;
double dinlightgro=0.0;
double flalightgro=0.0;
double ehulightgro=0.0;

double diagra=0.0;
double dingra=0.0;
double flagra=0.0;
double ehugra=0.0;
double micgra=0.0;

double callightgro=0.0;
double caltemgro=0.0;

double grazd=0.0;        // microzoo grazing on diatoms
double graze=0.0;        // microzoo grazing on Ehu

double esurf=0.0; 

double diff=0.0;   // cross-termocline mixing rate

double nbo=0.0;
double sbo=0.0;

double **y, *tt;   // required by rkdriver, for communicating back with main

double varTeh=0.0;

double varT=0.0;
double varM=0.0;
double varH=0.0;   // h+(t) = maximum(h(t), 0) as in FASH93, h(t) is d(mixed)/dt

double mixed=0.0;  // actual mixed layer depth as obtained from Levitus data

double mld[HSTEP]; // mixed layer depth
double mldo[HSTEP];

double *mldp95, *mldp95o, *sstp95, *win94;   // forcing prior 1995

double *mld95o, *mld96o, *mld97o, *mld98o, *mld99o, *mld00o, *mld01o;

double *mld95, *mld96, *mld97, *mld98, *mld99, *mld00, *mld01;
double *sst95, *sst96, *sst97, *sst98, *sst99, *sst00, *sst01;
double *par95, *par96, *par97, *par98, *par99, *par00, *par01;
double *win95, *win96, *win97, *win98, *win99, *win00, *win01;

double chltoc=0.0; // adaptive Chl:C ratio

double tem[HSTEP]; // temperature
double sal[HSTEP]; // salinity
double sir[HSTEP]; // irradiance at surface
double wsp[HSTEP]; // wind speed

double gtv=0.0;    // gas transfer velocity
double co2sol=0.0; // CO2 solubility in seawater
double pco2w=0.0;  // pCO2 in seawater
double co32=0.0;   // [CO3=] 
double o_cal=0.0;  // Omega calcite
double o_ara=0.0;  // Omega aragonite
double ph=0.0;     // pH
double bica=0.0;   // [HCO3-]
double co2aq=0.0;  // [CO2(aq)]

double psi=0.0;    // light limiting all phytoplankton
double psieh=0.0;  // light limiting Emiliania huxleyi
double psica=0.0;  // light limiting calcification

double li,gli,nc,gnc,glidf,gncdf,sc,gsc; 

double ber=0.0;
double los=0.0;

double yi1,yi2,yi3,yi4,yi5,yi6,yi7,yi8,yi9,yi10,yi11,yi12,yi13,yi14;  // initial reservois' values 


// open files for results
ofstream outinf("./results/info.dat");

// multi-year solution
ofstream out1("./results/diato.dat");
ofstream out2("./results/flage.dat");
ofstream out3("./results/nitra.dat");
ofstream out4("./results/silic.dat");
ofstream out5("./results/mesoz.dat");
ofstream out6("./results/detri.dat");
ofstream out7("./results/micro.dat");
ofstream out8("./results/dinof.dat");
ofstream out9("./results/ehuxl.dat");
ofstream out10("./results/ammon.dat");
ofstream out11("./results/acocc.dat");
ofstream out12("./results/fcocc.dat");
ofstream out13("./results/tdic.dat");
ofstream out14("./results/talk.dat");
ofstream out15("./results/pco2.dat");
ofstream out16("./results/co32.dat");
ofstream out17("./results/ocal.dat");
ofstream out18("./results/oara.dat");
ofstream out19("./results/tzoop.dat");
ofstream out20("./results/npratio.dat");

ofstream outa("./results/tempd.dat");
ofstream outb("./results/tempdf.dat");
ofstream outc("./results/airr.dat"); 
ofstream outd("./results/sirr.dat");

// one-year (last) solution
ofstream outl("./results/dia_d.dat");
ofstream outm("./results/fla_d.dat");
ofstream outn("./results/nit_d.dat");
ofstream outo("./results/mic_d.dat");
ofstream outp("./results/din_d.dat");
ofstream outq("./results/sil_d.dat");
ofstream outr("./results/mes_d.dat");
ofstream outs("./results/det_d.dat");
ofstream outw("./results/ehu_d.dat");
ofstream outx("./results/amm_d.dat");
ofstream outf("./results/aco_d.dat");
ofstream outg("./results/fco_d.dat");

ofstream outcp("./results/diagnoST-PROVA.dat");
ofstream outres("./results/resST-PROVA.dat");

ofstream outt("./results/phy_d.dat");
ofstream outu("./results/phyto.dat");
ofstream outy("./results/tzo_d.dat");

ofstream outv("./results/zp.dat");
ofstream outz("./results/ehc_d.dat");

ofstream outctochl("./results/cch_d.dat");
ofstream outluce("./results/ali_d.dat");

ofstream outdic("./results/dic_d.dat");
ofstream outalk("./results/alk_d.dat");
ofstream outpco("./results/pco_d.dat");
ofstream outco3("./results/co3_d.dat");
ofstream outoca("./results/oca_d.dat");
ofstream outora("./results/oar_d.dat");
ofstream outoph("./results/ph_d.dat");
ofstream outobi("./results/bic_d.dat");


ofstream outbe("./results/birth.dat");
ofstream outlo("./results/loss.dat");

ofstream outmi("./results/mixed.dat");

//=================================== MAIN ======================================


main()
{

  static int first_time = TRUE;  // for setting first year initial conditions


  // === LOAD INPUT FILES (MLD, TEMP, SAL, AND WIND SPEED VALUES) === 

  int h=0;  
  int hi=0;
  
  int t=0;

  mldp95=dvector(1,HSTEP);
  mldp95o=dvector(1,HSTEP);
  sstp95=dvector(1,HSTEP);
  win94=dvector(1,HSTEP);

  mld95=dvector(1,HSTEP);
  sst95=dvector(1,HSTEP);
  par95=dvector(1,HSTEP);
  win95=dvector(1,HSTEP);

  mld96=dvector(1,HSTEP);
  sst96=dvector(1,HSTEP);
  par96=dvector(1,HSTEP);
  win96=dvector(1,HSTEP);

  mld97=dvector(1,HSTEP);
  sst97=dvector(1,HSTEP);
  par97=dvector(1,HSTEP);
  win97=dvector(1,HSTEP);

  mld98=dvector(1,HSTEP);
  sst98=dvector(1,HSTEP);
  par98=dvector(1,HSTEP);
  win98=dvector(1,HSTEP);

  mld99=dvector(1,HSTEP);
  sst99=dvector(1,HSTEP);
  par99=dvector(1,HSTEP);
  win99=dvector(1,HSTEP);

  mld00=dvector(1,HSTEP);
  sst00=dvector(1,HSTEP);
  par00=dvector(1,HSTEP);
  win00=dvector(1,HSTEP);

  mld01=dvector(1,HSTEP);
  sst01=dvector(1,HSTEP);
  par01=dvector(1,HSTEP);
  win01=dvector(1,HSTEP);

  mld95o=dvector(1,HSTEP);
  mld96o=dvector(1,HSTEP);
  mld97o=dvector(1,HSTEP);
  mld98o=dvector(1,HSTEP);
  mld99o=dvector(1,HSTEP);
  mld00o=dvector(1,HSTEP);
  mld01o=dvector(1,HSTEP);

  //char mldp95f[20]=".input/mldnew2i.in";
  //char sstp95f[20]=".input/tem.in";

  char mldp95f[20]="./input/mldnew2i.in";
  char sstp95f[20]="./input/tem.in";
  char win94f[20]="./input/win94.in";


  char mld95f[20]="./input/mld95i2.in";
  char sst95f[20]="./input/sst95i.in";
  char par95f[20]="./input/par95n.in";
  char win95f[20]="./input/win95.in";

  char mld96f[20]="./input/mld96i.in";
  char sst96f[20]="./input/sst96i.in";
  char par96f[20]="./input/par96n.in";
  char win96f[20]="./input/win96.in";

  char mld97f[20]="./input/mld97i2.in";
  char sst97f[20]="./input/sst97i.in";
  char par97f[20]="./input/par97n.in";
  char win97f[20]="./input/win97.in";

  char mld98f[20]="./input/mld98i.in";
  char sst98f[20]="./input/sst98i.in";
  char par98f[20]="./input/par98n.in";
  char win98f[20]="./input/win98.in";

  char mld99f[20]="./input/mld99i.in";
  char sst99f[20]="./input/sst99i.in";
  char par99f[20]="./input/par99n.in";
  char win99f[20]="./input/win99.in";

  char mld00f[20]="./input/mld00i.in";
  char sst00f[20]="./input/sst00i.in";
  char par00f[20]="./input/par00n.in";
  char win00f[20]="./input/win00.in";

  char mld01f[20]="./input/mld01i2.in";
  char sst01f[20]="./input/sst01bi2.in";
  char par01f[20]="./input/par01n.in";
  char win01f[20]="./input/win01.in";


  //============== Forcing for post-1995 years =============
  
  if(trans){ 

    cout<<endl;
    cout<<"looking at transient"<<endl;
    cout<<endl;

    // === open 1995 ===
    ifstream inl1(mld95f);  
    if(!inl1){
      cout<<" Impossible to open 1995 MLD file\n";
      return 1;
    }
    ifstream inl2(sst95f); 
    if(!inl2){
      cout<<" Impossible to open 1995 TEM file\n";
      return 1;
    }
    ifstream inl3(par95f);  
    if(!inl3){
      cout<<" Impossible to open 1995 PAR file\n";
      return 1;
    }
    ifstream inlw1(win95f);  
    if(!inlw1){
      cout<<" Impossible to open 1995 WIN file\n";
      return 1;
    }

    // === open 1996 ===
    ifstream inl4(mld96f);  
    if(!inl4){
      cout<<" Impossible to open 1996 MLD file\n";
      return 1;
    }
    ifstream inl5(sst96f); 
    if(!inl5){
      cout<<" Impossible to open 1996 TEM file\n";
      return 1;
    }
    ifstream inl6(par96f);  
    if(!inl6){
      cout<<" Impossible to open 1996 PAR file\n";
      return 1;
    }
    ifstream inlw2(win96f);  
    if(!inlw2){
      cout<<" Impossible to open 1997 WIN file\n";
      return 1;
    }

    // === open 1997 ===
    ifstream inl7(mld97f);  
    if(!inl7){
      cout<<" Impossible to open 1997 MLD file\n";
      return 1;
    }
    ifstream inl8(sst97f); 
    if(!inl8){
      cout<<" Impossible to open 1997 TEM file\n";
      return 1;
    }
    ifstream inl9(par97f);  
    if(!inl9){
      cout<<" Impossible to open 1997 PAR file\n";
      return 1;
    }
    ifstream inlw3(win97f);  
    if(!inlw3){
      cout<<" Impossible to open 1997 WIN file\n";
      return 1;
    }

    // === open 1998 ===
    ifstream inl10(mld98f);  
    if(!inl10){
      cout<<" Impossible to open 1998 MLD file\n";
      return 1;
    }
    ifstream inl11(sst98f); 
    if(!inl11){
      cout<<" Impossible to open 1998 TEM file\n";
      return 1;
    }
    ifstream inl12(par98f);  
    if(!inl12){
      cout<<" Impossible to open 1998 PAR file\n";
      return 1;
    }
    ifstream inlw4(win95f);  
    if(!inlw4){
      cout<<" Impossible to open 1998 WIN file\n";
      return 1;
    }

    // === open 1999 ===
    ifstream inl13(mld99f);  
    if(!inl13){
      cout<<" Impossible to open 1999 MLD file\n";
      return 1;
    }
    ifstream inl14(sst99f); 
    if(!inl14){
      cout<<" Impossible to open 1999 TEM file\n";
      return 1;
    }
    ifstream inl15(par99f);  
    if(!inl15){
      cout<<" Impossible to open 1999 PAR file\n";
      return 1;
    }
    ifstream inlw5(win99f);  
    if(!inlw5){
      cout<<" Impossible to open 1999 WIN file\n";
      return 1;
    }

    // === open 2000 ===
    ifstream inl16(mld00f);  
    if(!inl16){
      cout<<" Impossible to open 2000 MLD file\n";
      return 1;
    }
    ifstream inl17(sst00f); 
    if(!inl17){
      cout<<" Impossible to open 2000 TEM file\n";
      return 1;
    }
    ifstream inl18(par00f);  
    if(!inl18){
      cout<<" Impossible to open 2000 PAR file\n";
      return 1;
    }
    ifstream inlw6(win95f);  
    if(!inlw6){
      cout<<" Impossible to open 2000 WIN file\n";
      return 1;
    }

    // === open 2001 ===
    ifstream inl19(mld01f);  
    if(!inl19){
      cout<<" Impossible to open 2001 MLD file\n";
      return 1;
    }
    ifstream inl20(sst01f); 
    if(!inl20){
      cout<<" Impossible to open 2001 TEM file\n";
      return 1;
    }
    ifstream inl21(par01f);  
    if(!inl21){
      cout<<" Impossible to open 2001 PAR file\n";
      return 1;
    }
    ifstream inlw7(win95f);  
    if(!inlw7){
      cout<<" Impossible to open 2001 WIN file\n";
      return 1;
    }

    // ==== load 1995 ====
    t=0;
    while(inl1){
      inl1>>mld95[t];
      //test
      mld95[t]=mld95[t];// + 3.0;
      mld95o[t]=mld95[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                 // calculating dM/dt, as in FASH93 at pag.493
      mld95[h]=(mld95[h+1]-mld95[h])/1.0;  // 1.0 is the time interval (1 hour)
    }
    t=0;
    while(inl2){
      inl2>>sst95[t];
      t++;
    }
    t=0;
    while(inl3){
      inl3>>par95[t];
      t++;
    }
    t=0;
    while(inlw1){
      inlw1>>win95[t];
      t++;
    }

    // ==== load 1996 ====
    t=0;
    while(inl4){
      inl4>>mld96[t];
      //test
      mld96[t]=mld96[t];// + 3.0;
      mld96o[t]=mld96[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                 // calculating dM/dt, as in FASH93 at pag.493
      mld96[h]=(mld96[h+1]-mld96[h])/1.0;  // 1.0 is the time interval (1 hour)
    }
    t=0;
    while(inl5){
      inl5>>sst96[t];
      t++;
    }
    t=0;
    while(inl6){
      inl6>>par96[t];
      t++;
    }
    t=0;
    while(inlw2){
      inlw2>>win96[t];
      t++;
    }

    // ==== load 1997 ====
    t=0;
    while(inl7){
      inl7>>mld97[t];
      //test
      mld97[t]=mld97[t];// - 5.0;
      mld97o[t]=mld97[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                 // calculating dM/dt, as in FASH93 at pag.493
      mld97[h]=(mld97[h+1]-mld97[h])/1.0;  // 1.0 is the time interval (1 hour)
    }
    t=0;
    while(inl8){
      inl8>>sst97[t];
      t++;
    }
    t=0;
    while(inl9){
      inl9>>par97[t];
      t++;
    }
    t=0;
    while(inlw3){
      inlw3>>win97[t];
      t++;
    }

    // ==== load 1998 ====
    t=0;
    while(inl10){
      inl10>>mld98[t];
      //test
      mld98[t]=mld98[t];// + 3.0;
      mld98o[t]=mld98[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                 // calculating dM/dt, as in FASH93 at pag.493
      mld98[h]=(mld98[h+1]-mld98[h])/1.0;  // 1.0 is the time interval (1 hour)
    }
    t=0;
    while(inl11){
      inl11>>sst98[t];
      t++;
    }
    t=0;
    while(inl12){
      inl12>>par98[t];
      t++;
    }
    t=0;
    while(inlw4){
      inlw4>>win98[t];
      t++;
    }

    // ==== load 1999 ====
    t=0;
    while(inl13){
      inl13>>mld99[t];
      //test
      mld99[t]=mld99[t];// + 3.0;
      mld99o[t]=mld99[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                 // calculating dM/dt, as in FASH93 at pag.493
      mld99[h]=(mld99[h+1]-mld99[h])/1.0;  // 1.0 is the time interval (1 hour)
    }
    t=0;
    while(inl14){
      inl14>>sst99[t];
      t++;
    }
    t=0;
    while(inl15){
      inl15>>par99[t];
      t++;
    }
    t=0;
    while(inlw5){
      inlw5>>win99[t];
      t++;
    }

    // ==== load 2000 ====
    t=0;
    while(inl16){
      inl16>>mld00[t];
      //test
      mld00[t]=mld00[t];// + 3.0;
      mld00o[t]=mld00[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                 // calculating dM/dt, as in FASH93 at pag.493
      mld00[h]=(mld00[h+1]-mld00[h])/1.0;  // 1.0 is the time interval (1 hour)
    }
    t=0;
    while(inl17){
      inl17>>sst00[t];
      t++;
    }
    t=0;
    while(inl18){
      inl18>>par00[t];
      t++;
    }
    t=0;
    while(inlw6){
      inlw6>>win00[t];
      t++;
    }

    // ==== load 2001 ====
    t=0;
    while(inl19){
      inl19>>mld01[t];
      //test
      mld01[t]=mld01[t];// + 3.0;
      mld01o[t]=mld01[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                 // calculating dM/dt, as in FASH93 at pag.493
      mld01[h]=(mld01[h+1]-mld01[h])/1.0;  // 1.0 is the time interval (1 hour)
    }
    t=0;
    while(inl20){
      inl20>>sst01[t];
      t++;
    }
    t=0;
    while(inl21){
      inl21>>par01[t];
      t++;
    }
    t=0;
    while(inlw7){
      inlw7>>win01[t];
      t++;
    }

    // === close input files ===

    inl1.close();
    inl2.close();
    inl3.close();
    inl4.close();
    inl5.close();
    inl6.close();
    inl7.close();
    inl8.close();
    inl9.close();
    inl10.close();
    inl11.close();
    inl12.close();
    inl13.close();
    inl14.close();
    inl15.close();
    inl16.close();
    inl17.close();
    inl18.close();
    inl19.close();
    inl20.close();
    inl21.close();
    inlw1.close();
    inlw2.close();
    inlw3.close();
    inlw4.close();
    inlw5.close();
    inlw6.close();
    inlw7.close();

  }
  else{

    cout<<endl;
    cout<<"looking at steady-state"<<endl;
    cout<<endl;

    ifstream in1(mld96f);  //in1("mldnew.dat");
    if(!in1){
      cout<<" Impossible to open MLD file\n";
      return 1;
    }
    ifstream in2(sst96f);  //in2("tem.dat");  
    if(!in2){
      cout<<" Impossible to open TEM file\n";
      return 1;
    }   
    ifstream in5(par96f);
    if(!in5){
      cout<<" Impossible to open PAR file\n";
      return 1;
    }
    ifstream inw(win96f);
    if(!inw){
      cout<<" Impossible to open WIN file\n";
      return 1;
    }

    t=0;
    while(in1){
      in1>>mld96[t];
      mld96o[t]=mld96[t];
      t++;
    }
    t=0;
    while(in2){
      in2>>sst96[t];
      //tem[t]+=5;
      t++;
    }
    t=0;
    while(in5){
      in5>>par96[t];
      t++;
    }
   t=0;
    while(inw){
      inw>>win96[t];
      t++;
    }
    for(h=0;h<=HSTEP;h++){                // calculating dM/dt, as in FASH93 at pag.493
      mld96[h]=(mld96[h+1]-mld96[h])/1.0; // 1.0 is the time interval (1 hour)
    }
  
    in1.close();
    in2.close();
    in5.close();
    inw.close();
  }

  //==================================================

  // ==== open pre-1995 ====
  ifstream in6(mldp95f);
  if(!in6){
    cout<<" Impossible to open pre-1995 MLD file\n";
    return 1;
  }
  ifstream in7(sstp95f);
  if(!in7){
    cout<<" Impossible to open pre-1995 SST file\n";
    return 1;
  }
  ifstream in8(par95f);
  if(!in8){
    cout<<" Impossible to open pre-1995 PAR file\n";
    return 1;
  }
  ifstream in9(win94f);
  if(!in9){
    cout<<" Impossible to open pre-1995 WIN file\n";
    return 1;
  }

  // ==== load pre-1995 ====
  t=0;
  while(in6){
    in6>>mldp95[t];
    mldp95o[t]=mldp95[t];
    t++;
  }
  t=0;
  while(in7){
    in7>>sstp95[t];
    t++;
  }
  t=0;
  while(in8){
    in8>>par95[t];
    t++;
  }
  t=0;
  while(in9){
    in9>>win94[t];
    t++;
  }


  for(h=0;h<=HSTEP;h++){                     // calculating dM/dt, as in FASH93 at pag.493
      mldp95[h]=(mldp95[h+1]-mldp95[h])/1.0; // 1.0 is the time interval (1 hour)
  }


  // == THESE SAL AND WSP SHOULD BE DELETED NOW ==

  // == open sal and WSP ==
  ifstream in3("./input/sal.in");
  if(!in3){
    cout<<" Impossible to open SAL file\n";
    return 1;
  }
  ifstream in4("./input/wsp.in");
  if(!in4){
    cout<<" Impossible to open WSP file\n";
    return 1;
  }
  // == load sal and WSP ==
  t=0;
  while(in3){
    in3>>sal[t];
    t++;
  }
  t=0;
  while(in4){
    in4>>wsp[t];
    t++;
  }
  
  in3.close();
  in4.close();

  // =============================================

  in6.close();
  in7.close();
  in8.close();
  in9.close();

  // =================================================
  
  
  int i=0;
  double *vstart;
  
  tt=dvector(1,HSTEP);
  y=dmatrix(1,NEQ,1,HSTEP);
  vstart=dvector(1,NEQ);
  
  
  // ======== set initial conditions ===========
  
  
  double i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14;
  
  if(first_time){  // setting first year initial conditions

    first_time = FALSE;
    
    i1=0.01;   // dia 0.01
    i2=0.01;   // fla 0.01
    i3=20.0;   // 15 nit  
    i4=35.0;   // 30 sil 
    i5=0.01;   // mes
    i6=0.05;   // det
    i7=0.01;   // mic
    i8=0.01;   // din 0.01
    i9=0.01;   //0.01;   // ehu 0.01
    i10=0.0001;// amm  (WHIT99, in Dynamics of Bering Sea) 
    i11=0.3;   //0.3;   // aco 0.3 (about 30 times the concentration of Ehux)
    i12=0.0001;//0.0001; // fco
    i13=2100.0;// dic  
    i14=2250.0;// alk

    outres<<"#jday "<<"month "<<"temp "<<"MLD "<<"sal "<<"Irr "<<"diato "<<"flage "<<"dino "
          <<"ehux "<<"microz "<<"mesoz "<<"totphy "<<"nit "<<"ammo "<<"sil "<<"DIC "<<"Alk "
	  <<"pCO2 "<<"CO3 "<<"omegacal "<<"omegaara "<<"acocco "<<"fcocco "<<"CO2(aq) "<<"HCO3 "
          <<"totzoo "<<endl;

    outcp<<"#jday  "<<"month  "<<"Pho:Cal ratio  "<<"f-ratio  "<<"Tot phy biomass  "
         <<"Tot phy prod (phyto growth terms)  "<<"PON  "<<"Tot zoo biomass  "
         <<"C:Chl ratio  "<<"TAlk  "<<"Sal"<<endl;
  }
  
  vstart[1]=i1;   // [1] - diatoms
  vstart[2]=i2;   // [2] - flagellates
  vstart[3]=i3;   // [3] - nitrogen
  vstart[4]=i4;   // [4] - silicate
  vstart[5]=i5;   // [5] - mesozooplankton
  vstart[6]=i6;   // [6] - detritus
  vstart[7]=i7;   // [7] - microzooplankton
  vstart[8]=i8;   // [8] - dinoflagellates
  vstart[9]=i9;   // [9] - ehux
  vstart[10]=i10; // [10] - ammonia
  vstart[11]=i11; // [11] - attached coccoliths
  vstart[12]=i12; // [12] - free coccoliths
  vstart[13]=i13; // [13] - total dissolved inorganic carbon
  vstart[14]=i14; // [14] - total alkalinity

  //vstart[10]=0.0;  // [10] - silicate mass balance management
  //vstart[11]=0.0;  // [11] - nitrogen mass balance management
  
  yi1=vstart[1];
  yi2=vstart[2];
  yi3=vstart[3];
  yi4=vstart[4];
  yi5=vstart[5];
  yi6=vstart[6];
  yi7=vstart[7];
  yi8=vstart[8];
  yi9=vstart[9];
  yi10=vstart[10];
  yi11=vstart[11];
  yi12=vstart[12];
  yi13=vstart[13];
  yi14=vstart[14];


  // ===========================================
  
  
  rkdriver(vstart,NEQ,TI,TH,HSTEP,derivs);   
  

  // ==== free all vectors ====
  
  free_dmatrix(y,1,NEQ,1,HSTEP);
  free_dvector(tt,1,HSTEP);
  free_dvector(vstart,1,NEQ);
  
  free_dvector(mldp95,1,HSTEP);
  free_dvector(mldp95o,1,HSTEP);
  free_dvector(sstp95,1,HSTEP);
  free_dvector(win94,1,HSTEP);

  free_dvector(mld95,1,HSTEP);
  free_dvector(sst95,1,HSTEP);
  free_dvector(par95,1,HSTEP);
  free_dvector(win95,1,HSTEP);

  free_dvector(mld96,1,HSTEP);
  free_dvector(sst96,1,HSTEP);
  free_dvector(par96,1,HSTEP);
  free_dvector(win96,1,HSTEP);

  free_dvector(mld97,1,HSTEP);
  free_dvector(sst97,1,HSTEP);
  free_dvector(par97,1,HSTEP);
  free_dvector(win97,1,HSTEP);

  free_dvector(mld98,1,HSTEP);
  free_dvector(sst98,1,HSTEP);
  free_dvector(par98,1,HSTEP);
  free_dvector(win98,1,HSTEP);

  free_dvector(mld99,1,HSTEP);
  free_dvector(sst99,1,HSTEP);
  free_dvector(par99,1,HSTEP);
  free_dvector(win99,1,HSTEP);

  free_dvector(mld00,1,HSTEP);
  free_dvector(sst00,1,HSTEP);
  free_dvector(par00,1,HSTEP);
  free_dvector(win00,1,HSTEP);

  free_dvector(mld01,1,HSTEP);
  free_dvector(sst01,1,HSTEP);
  free_dvector(par01,1,HSTEP);
  free_dvector(win01,1,HSTEP);

  free_dvector(mld95o,1,HSTEP);
  free_dvector(mld96o,1,HSTEP);
  free_dvector(mld97o,1,HSTEP);
  free_dvector(mld98o,1,HSTEP);
  free_dvector(mld99o,1,HSTEP);
  free_dvector(mld00o,1,HSTEP);
  free_dvector(mld01o,1,HSTEP);


  // ===== close all files =====
   
  outinf.close();
  
  out1.close();
  out2.close();
  out3.close();
  out4.close();
  out5.close();
  out6.close();
  out7.close();
  out8.close();
  out9.close();
  out10.close();
  out11.close();
  out12.close();
  out13.close();
  out14.close();
  out15.close();
  out16.close();
  out17.close();
  out18.close();
  out19.close();
  out20.close();

  outa.close();
  outb.close();
  outc.close();
  outd.close();
  
  outl.close();
  outm.close();
  outn.close();
  outo.close();
  outp.close();
  outq.close();
  outr.close();
  outs.close();
  outt.close();
  outu.close();
  outw.close();
  outv.close();
  outz.close();
  outx.close();
  outf.close();
  outg.close();
  outy.close();
  outcp.close();

  outres.close();

  outctochl.close();
  outluce.close();

  outdic.close();
  outalk.close();
  outpco.close();
  outco3.close();
  outoca.close();
  outora.close();
  outoph.close();
  outobi.close();


  outbe.close();
  outlo.close();

  outmi.close();

  return 0;  
  
}


//============================= ODE ROUTINES ================================


void rkdriver(double vstart[], int nvar, double t1, double t2, int nstep, 
	      void (*derivs)(double, double [], double []))
{

  int i,k;

  // temporary state variables
  double dia=0.0;
  double fla=0.0;
  double din=0.0;
  double ehu=0.0;
  double mic=0.0;
  double mes=0.0;

  double luce=0.0;

  double chlo=0.0;   // total chlorophyll - to feed in into the light routine (in mg Chl/m3)  
  
  // adaptive Chl:C ratios
  double chlcd=0.0;  // dia
  double chlcf=0.0;  // fla
  double chlcdf=0.0; // din
  double chlceh=0.0; // ehu
  
  // nutrient-dependent growth rates  
  double mud=0.0;    // dia
  double muf=0.0;    // fla
  double mudf=0.0;   // din
  double mueh=0.0;   // ehu

  double nit=0.0;    // nitrate  - to feed into the Chl:C routine
  double amm=0.0;    // ammonium - to feed into the Chl:C routine

  double sil=0.0;    // silicate - to feed into the carbonate routines   
  double tco2=0.0;   // TCO2 - to feed into the carbonate routines
  double alk=0.0;    // Alkalinity - to feed into the carbonate routines
  double temp=0.0;   // temperature - to feed into the carbonate routines
  double salin=0.0;  // salinity - to feed into the carbonate routines
  double wspeed=0.0; // wind speed - to feed into the carbonate routines
  
  double t,h;
  double *v,*vout,*dv;

  //int yy;                // actual year

  for(yy=0;yy<=Y;yy++){  // number of years

    diff=mm;

    if(yy==3) diff=mm95;
    if(yy==4) diff=mm96;
    if(yy==5) diff=mm97;
    if(yy==6) diff=mm98;
    if(yy==7) diff=mm99;
    if(yy==8) diff=mm00;
    if(yy==9) diff=mm01;

    nbo=N0;
    sbo=S0;

    if(yy==2){
      nbo=N094;
      sbo=S094;
    }

    if(yy==4){
      nbo=N095;
      sbo=S095;
    }

    if(yy==4){
      nbo=N096;
      sbo=S096;
    }

    if(yy==5){
      nbo=N097;
      sbo=S097;
    }

    if(yy==6){
      nbo=N098;
      sbo=S098;
    }

    if(yy==6){
      nbo=N099;
      sbo=S099;
    }

    if(yy==7){
      nbo=N000;
      sbo=S000;
    }

    if(yy==8){
      nbo=N001;
      sbo=S001;
    }


    // The model is run for a number of years to stabilize it
    // and in order to look at steady-state results. But
    // when a transient result is needed (last year forcing
    // different than all previous years) then the variable 
    // 'trans' is set TRUE and different MLD and TEM forcing 
    // functions are used for last year run. 

    // === transient ===
    if(yy<Y-6 && trans){
      cout<<" year before 1995"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mldp95[i];
	mldo[i]=mldp95o[i];
	tem[i]=sstp95[i];
	if(i>2880 && i<6720) sir[i]=par95[i];// + 0.0;
	else sir[i]=par95[i];
	sir[i]=par95[i];
	wsp[i]=win94[i];
      }
    }
    if(yy==Y-6 && trans){
      cout<<" year 1995"<<endl;    
      for(i=0;i<HSTEP;i++){
    	mld[i]=mld95[i];
    	mldo[i]=mld95o[i];
    	tem[i]=sst95[i];
	if(i>2880 && i<6720) sir[i]=par95[i];// - 8.0;
    	else sir[i]=par95[i];
	wsp[i]=win95[i];
      }
    }
    if(yy==Y-5 && trans){
      cout<<" year 1996"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mld96[i];
	mldo[i]=mld96o[i];
	tem[i]=sst96[i];
	if(i>2880 && i<6720) sir[i]=par96[i];// + 3.0;
	else sir[i]=par96[i];
	wsp[i]=win96[i];
      }
    }
    if(yy==Y-4 && trans){
      cout<<" year 1997"<<endl;    
      for(i=0;i<HSTEP;i++){
    	mld[i]=mld97[i];
    	mldo[i]=mld97o[i];
    	tem[i]=sst97[i];
    	if(i>2880 && i<6720) sir[i]=par97[i];// + 10.0;
	else sir[i]=par97[i];
	wsp[i]=win97[i];
      }
    }
    if(yy==Y-3 && trans){
      cout<<" year 1998"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mld98[i];
	mldo[i]=mld98o[i];
	tem[i]=sst98[i];
	if(i>2880) sir[i]=par98[i];// + 8.0;
	else sir[i]=par98[i];
	wsp[i]=win98[i];
      }
    }
    if(yy==Y-2 && trans){
      cout<<" year 1999"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mld99[i];
	mldo[i]=mld99o[i];
	tem[i]=sst99[i];
	if(i>2880 && i<6720) sir[i]=par99[i];// + 8.0;
	else  sir[i]=par99[i];
	wsp[i]=win99[i];
      }
    }
    if(yy==Y-1 && trans){
      cout<<" year 2000"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mld00[i];
	mldo[i]=mld00o[i];
	tem[i]=sst00[i];
	if(i>2880 && i<6720) sir[i]=par00[i];// + 10.0;
	else sir[i]=par00[i];
	wsp[i]=win00[i];
      }
    }
    if(yy==Y && trans){
      cout<<" year 2001"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mld01[i];
	mldo[i]=mld01o[i];
	tem[i]=sst01[i];
	if(i>2880 && i<6720) sir[i]=par01[i];// - 5.0;
	else sir[i]=par01[i];
	wsp[i]=win01[i];
      }
    }

    // === steady-state ===
    if(yy<Y && !trans){
      cout<<" year before 1995"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mldp95[i];
	mldo[i]=mldp95o[i];	
	tem[i]=sstp95[i];
	sir[i]=par95[i];
	wsp[i]=win94[i];
      }
    }
    if(yy==Y && !trans){
      cout<<" year 1996"<<endl;
      for(i=0;i<HSTEP;i++){
	mld[i]=mld96[i];
	mldo[i]=mld96o[i];
	tem[i]=sst96[i];
	sir[i]=par96[i];
	wsp[i]=win96[i];
      }
    }

    v=dvector(1,nvar);
    vout=dvector(1,nvar);
    dv=dvector(1,nvar);
    
    // note: nvar is the number of ODEs (i.e. NEQ)
    for(i=1;i<=nvar;i++){   // loading starting values
      v[i]=vstart[i];
      y[i][1]=v[i];
    }

    //cout<<"\n";
    //cout<<"diato  "<<vstart[1]<<"\n";
    //cout<<"dinof  "<<vstart[2]<<"\n";
    //cout<<"nitra  "<<vstart[3]<<"\n";
    //cout<<"silic  "<<vstart[4]<<"\n";
    //cout<<"zoopl  "<<vstart[5]<<"\n";
    //cout<<"detri  "<<vstart[6]<<"\n";
    //cout<<"\n";
    
    tt[1]=t1;
    t=t1;

    h=(t2-t1)/nstep;
    
    int conta=0;

    if(yy==0) chlcd=chlcdf=chlcf=chlceh=CHLTOC; //0.025;  // initial value for Chl:C ratio
      
    chlo=NTOC*(chlcd*v[1]+chlcdf*v[2]+chlcf*v[8]+chlceh*v[9]); // total chlorophyll in mg Chl/m3       

    nit=v[3];
    sil=v[4];
    amm=v[10];
    tco2=v[13];    
    alk=v[14];

    for(k=1;k<=nstep-3;k++){    // take nstep (for ex.: 8760, when in h-1) steps

      // ================= light system ==================

      varH=mld[k+1];                        // mixed layer depth variation, h(t)=dM/dt as in FASH93
      mixed=mldo[k+1];                      // mixed layer depth, M(t) in FASH93

      //esurf=get_light_at_surface(k+1);    // CALCULATED light at surface at time k of the year
      esurf=sir[k+1];                       // MEASURED   light at surface at time k of the year
    
      outd<<(k+1)<<"   "<<esurf<<endl;      // save light at surface (in W m-2)

      psi=get_averaged_light(esurf,chlo,mixed);      // light limitation for all phytopl either than Ehux
      psieh=get_averaged_light_eh(esurf,chlo,mixed); // light limitation for E. huxleyi
      psica=get_averaged_light_cal(esurf,chlo,mixed);// light limitation for Calcification

      li=get_light_intensity(esurf,chlo);           // light at a given depth (5 m)

      // ================ carbonate system ================

      temp=tem[k];
      salin=sal[k];      
      wspeed=wsp[k];

      gtv=get_gas_transfer_velocity(wspeed,temp);        // get gas transfer velocity
      co2sol=get_co2_solubility(salin,temp);             // get CO2 solubility 

      pco2w=get_param_water(salin,temp,alk,tco2,sil,1);  // get pCO2 in water
      co32=get_param_water(salin,temp,alk,tco2,sil,2);   // get [CO3=] 
      o_cal=get_param_water(salin,temp,alk,tco2,sil,3);  // get omega-calcite
      o_ara=get_param_water(salin,temp,alk,tco2,sil,4);  // get omega-aragonite
      ph=get_param_water(salin,temp,alk,tco2,sil,5);     // get pH
      bica=get_param_water(salin,temp,alk,tco2,sil,6);   // get [HCO3-]
      co2aq=get_param_water(salin,temp,alk,tco2,sil,7);  // get [CO2(aq)]

      //ingEH=90.0/exp(o_cal*o_cal);
      ingEH=10.0/(o_cal*o_cal*o_cal*o_cal); //16.45
      //MEH=1.2/(1.3*o_cal*o_cal*o_cal);//f(x)=1.2/(1.3*x*x*x)
      //ingEH=90000.0/(o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal*o_cal); 

      // =============== temperature system ===============

      varT=exp(0.063*temp);          // compute growth limitation with temperature (EPPL72)
      varTeh=exp(0.063*temp);

      // =================== Chl:C system =================  // Cloern et al. 1995 L&O:40(7) 1313-1321
      
      mud=min((nit/NHD+amm/AHD)/(1+nit/NHD+amm/AHD), sil/(SH+sil));
      muf=(nit/NHF+amm/AHDF)/(1+nit/NHF+amm/AHF);
      mudf=(nit/NHDF+amm/AHDF)/(1+nit/NHF+amm/AHDF);
      mueh=(nit/NHEH+amm/AHEH)/(1+nit/NHEH+amm/AHEH);

      // adapted Chl:C ratio (CLOE95)
      // luce=esurf/(KW*mixed)*(1-exp(-KW*mixed)); // 0.4 is to transform Wm-2 into mol quanta m-2 d-1

      //luce=get_light(esurf,chlo,mixed);                          // average light in the MLD  
      //luce=0.3*luce;     // 0.34 is to transfotm luce from W m-2 to mol quanta m-2 d-1
      //chlcd = 0.003+0.0154*exp(0.05*temp)*exp(-0.059*luce)*mud   // Chl:C in diatoms (24 for units of h)
      //chlcf = 0.003+0.0154*exp(0.05*temp)*exp(-0.059*luce)*muf;  // Chl:C in flagellates
      //chlcdf = 0.003+0.0154*exp(0.05*temp)*exp(-0.059*luce)*mudf;// Chl:C in dinoflagellates
      //chlceh = 0.003+0.0154*exp(0.05*temp)*exp(-0.059*luce)*mueh;// Chl:C in E. huxleyi

      // constant Chl:C ratio
      chlcd=CHLTOC;
      chlcf=CHLTOC;
      chlcdf=CHLTOC;
      chlceh=CHLTOC;

      // ==================================================


      (*derivs)(t,v,dv);      
      rk4(v,dv,nvar,t,h,vout,derivs);
      
      if((double)(t+h) == t) nrerror(" Step size too small in routine rkdriver ");
      t+=h;
      tt[k+1]=t;              // store intermediate steps

      for(i=1;i<=14;i++) vout[i]=fabs(vout[i]);

      for(i=1;i<=nvar;i++){ 
	v[i]=vout[i];
	y[i][k+1]=v[i];      
	
	if(fmod(k,24)==0 && yy==Y){
	  nc=v[3];
	  sc=v[4];
	}

	// (tt[k+1]+HSTEP*yy)/TIME
	if(i==1){ 
	  if(fmod(k,24)==0) out1<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save diatoms   	
	  dia=y[i][k+1];
	}
	if(i==2){
	  if(fmod(k,24)==0) out2<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save flage
	  fla=y[i][k+1];	
	}
	if(i==3){
	  if(fmod(k,24)==0) out3<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save nitrate	
	}
	if(i==4){
	  if(fmod(k,24)==0) out4<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save silicate 	
	}      
	if(i==5){
	  if(fmod(k,24)==0) out5<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save mesozoo
	  mes=y[i][k+1];
	}
	if(i==6){
	  if(fmod(k,24)==0) out6<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save detritus
	}
     	if(i==7){
	  if(fmod(k,24)==0) out7<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save microzoo
	  mic=y[i][k+1];
	}
	if(i==8){
	  if(fmod(k,24)==0) out8<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save dinofla
	  din=y[i][k+1];
	}
	if(i==9){
	  if(fmod(k,24)==0) out9<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save ehux
	  ehu=y[i][k+1];
	}
	if(i==10){
	  if(fmod(k,24)==0) out10<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save ammonia
	}
	if(i==11){
	  if(fmod(k,24)==0) out11<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save acocc
	}
	if(i==12){
	  if(fmod(k,24)==0) out12<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save fcocc
	}
	if(i==13){
	  if(fmod(k,24)==0) out13<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save tdic
	}
	if(i==14){
	  if(fmod(k,24)==0) out14<<tt[k+1]+HSTEP*yy<<"  "<<y[i][k+1]<<endl; // save talk
	}
      }


      //if(fmod(k,24)==0){ // start saving since firts year

      if(yy>2 && fmod(k,24)==0){  // start saving after third-year run 

      //if(yy==Y && fmod(k,24)==0){  // start saving after year before last

      //if(yy==Y){

	// NTOC = 12 * CTON  (CTON = 6.625)
	outl<<tt[k+1]/24.0<<"  "<<y[1][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // diatom (mmol N m-3)
	outm<<tt[k+1]/24.0<<"  "<<y[2][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // flagel (mmol N m-3)
	outn<<tt[k+1]/24.0<<"  "<<y[3][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // nitrat (mmol N m-3) 
	outq<<tt[k+1]/24.0<<"  "<<y[4][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // silica (mmol Si m-3)
	outr<<tt[k+1]/24.0<<"  "<<y[5][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // mesozo (mmol N m-3) 
	outs<<tt[k+1]/24.0<<"  "<<y[6][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // detrit (mmol N m-3)
	outo<<tt[k+1]/24.0<<"  "<<y[7][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // microz (mmol N m-3)
	outp<<tt[k+1]/24.0<<"  "<<y[8][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // dinofl (mmol N m-3)
	outw<<tt[k+1]/24.0<<"  "<<y[9][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;        // ehuxle (mmol N m-3)
	outx<<tt[k+1]/24.0<<"  "<<y[10][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;       // ammoni (mmol N m-3)
	outz<<tt[k+1]/24.0<<"  "<<12*CTON*y[9][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;// ehuxle (mg org-C m-3)
	outf<<tt[k+1]/24.0<<"  "<<12*y[11][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;    // at coc (mg cal-C m-3)
	outg<<tt[k+1]/24.0<<"  "<<12*y[12][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;    // fr coc (mg cal-C m-3)
	outdic<<tt[k+1]/24.0<<"  "<<y[13][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;     // Total CO2 (umol C kg-1)
	outalk<<tt[k+1]/24.0<<"  "<<y[14][k+1]<<"  "<<tt[k+1]*0.00126+1<<endl;     // Total Alk (uEq kg-1)
	outctochl<<tt[k+1]/24.0<<"  "<<1.0/chlcd<<"  "<<1.0/chlcdf<<"  "<<1.0/chlcf
                 <<"  "<<1.0/chlceh<<"  "<<tt[k+1]*0.00126+1<<endl;  // C:Chl seasonal ratio

	outluce<<tt[k+1]/24.0<<"  "<<luce<<"  "<<tt[k+1]*0.00126+1<<endl;


	// ============= DIAGNOSTIC OUTPUT =============
	//
	// to obtain PP in units of mmol C m-2 d-1 multiply by 24.0*mixed*CTON
	outcp<<tt[k+1]/24.0<<"  "<<tt[k+1]+HSTEP*yy<<"  "<<calcieh/(CTON*photoeh)
             <<"  "<<newphypro/(newphypro+regphypro)<<"  "<<totphypro<<"  "<<(newphypro+regphypro)
             <<"  "<<pon<<"  "<<totzoopro<<"  "<<1.0/(chlcd+chlcf+chlcdf+chlceh)*4
             <<"  "<<y[14][k+1]<<"  "<<salin<<"  "<<newphypro<<"  "<<regphypro
             <<"  "<<totphyloss<<"  "<<totzooloss<<"  "<<totphymix
             <<"  "<<dialightgro<<"  "<<dinlightgro<<"  "<<flalightgro<<"  "<<ehulightgro
             <<"  "<<dianutgro<<"  "<<dinnutgro<<"  "<<flanutgro<<"  "<<ehunutgro
             <<"  "<<callightgro<<"  "<<caltemgro<<"  "<<diagra<<"  "<<dingra<<"  "<<flagra
             <<"  "<<ehugra<<"  "<<micgra<<endl;

	outinf<<regdiapro<<"  "<<regdinpro<<"  "<<regflapro<<"  "<<regehupro<<"  "<<amm<<endl;

	// N:P ratio
	//       jday                  month                      N/P                     P                N
	out20<<tt[k+1]/24<<"  "<<tt[k+1]*0.00126+1<<"  "<<y[3][k+1]/y[10][k+1]<<"  "<<y[10][k+1]<<"  "<<y[3][k+1]<<endl;

	// Water pCO2
	outpco<<tt[k+1]/24.0<<"  "<<pco2w*1.0e6<<"  "<<tt[k+1]*0.00126+1<<endl; // water pCO2 (uatm) 
	outco3<<tt[k+1]/24.0<<"  "<<co32<<"  "<<tt[k+1]*0.00126+1<<endl;        // [CO3=] (umol C kg-1)
	outoca<<tt[k+1]/24.0<<"  "<<o_cal<<"  "<<tt[k+1]*0.00126+1<<endl;       // omega-calcite 
	outora<<tt[k+1]/24.0<<"  "<<o_ara<<"  "<<tt[k+1]*0.00126+1<<endl;       // omega-aragonite
        outoph<<tt[k+1]/24.0<<"  "<<ph<<"  "<<tt[k+1]*0.00126+1<<endl;          // pH 
        outobi<<tt[k+1]/24.0<<"  "<<bica<<"  "<<tt[k+1]*0.00126+1<<endl;        // [HCO3-] (umol C kg-1)


	// ================ MAIN OUTPUT ================
	//
	outres<<tt[k+1]/24.0<<"  "<<tt[k+1]+HSTEP*yy<<"  "<<temp<<"  "<<mixed<<"  "<<salin<<"  "<<esurf
               <<"  "<<wspeed<<"  "<<NTOC*chlcd*y[1][k+1]<<"  "<<NTOC*chlcf*y[2][k+1]<<"  "<<NTOC*chlcdf*y[8][k+1]
               <<"  "<<NTOC*chlceh*y[9][k+1]<<"  "<<NTOCZ*y[7][k+1]<<"  "<<NTOCZ*y[5][k+1]<<"  "
	       <<NTOC*(chlcd*y[1][k+1]+chlcdf*y[2][k+1]+chlcf*y[8][k+1]+chlceh*y[9][k+1])<<"  "
               <<"  "<<y[3][k+1]<<"  "<<y[10][k+1]<<"  "<<y[4][k+1]<<"  "<<y[13][k+1]<<"  "
               <<y[14][k+1]<<"  "<<pco2w*1.0e6<<"  "<<co32<<"  "<<o_cal<<"  "<<o_ara<<"  "
	       <<12*y[11][k+1]<<"   "<<12*y[12][k+1]<<"  "<<co2aq<<"  "<<bica<<"  "
	       <<NTOCZ*(y[5][k+1]+y[7][k+1])<<"  "<<grazd<<"  "<<graze<<"  "<<ph<<"  "
	      <<ingDI<<"  "<<ingEH<<"  "<<wspeed<<"  "<<gtv<<"  "<<(y[3][k+1]+y[10][k+1])<<endl;

	// total phytoplankton in: ug Chl L-1 (assumed = mg Chl m-3)
	// From:
	// C:N = (106/16)*12 = 79.5 (NTOC in param.h - 12 is to go from mol to weight)
	// Chl:C is calculated by adaptation to light, temperature and 
	// nutrient growth rate (Cloern et al., 1995)
	// Therefore to transform units from N to Chl use factor:
	// (C:N)*(Chl:C) =  79.5*Chl:C
        outt<<tt[k+1]/24<<"  "<<NTOC*(chlcd*y[1][k+1]+chlcdf*y[2][k+1]+chlcf*y[8][k+1]+chlceh*y[9][k+1])
            <<"  "<<tt[k+1]*0.00126+1<<endl; 

	// total zooplankton in: ug C/l (assumed = mg C/m3)
	outy<<tt[k+1]/24<<"  "<<NTOCZ*(y[5][k+1]+y[7][k+1])<<"  "<<tt[k+1]*0.00126+1<<endl;

	//outa<<tt[k+1]/24<<"   "<<MUD0*varT*24<<endl;   // save max growth vs. time for diatoms
	//outb<<tt[k+1]/24<<"   "<<MUDF0*varT*24<<endl;  // save max growth vs. time for flagellates
	outa<<temp<<"   "<<MUD0*varT*24<<endl;     // save max growth vs. temperature for diatoms
	outb<<temp<<"   "<<MUDF0*varT*24<<endl;    // save max growth vs. temperature for flagellates
	
	outc<<tt[k+1]/24<<"   "<<tt[k+1]*0.00126+1<<"  "<<psi<<endl;    // averaged light intensity vs. time 
	outmi<<tt[k+1]/24<<"  "<<tt[k+1]*0.00126+1<<"   "<<mixed<<endl; // mixed layer depth
	//outd<<tt[k+1]/24<<"  "<<tt[k+1]*0.00126+1<<"   "<<esurf/4.17<<endl;// light at surf (W m-2)  
      
      }                               
                            
      // save multi-year results daily
      if(fmod(k,24)==0){
	//outd<<tt[k+1]+HSTEP*yy)/TIME<<"   "<<esurf/4.17<<endl;// light at surf (W m-2)
	outu<<(tt[k+1]+HSTEP*yy)/TIME<<"  "<<dia+fla+din+ehu<<endl;  // save total phyto in mmol N m-3
	out19<<(tt[k+1]+HSTEP*yy)/TIME<<"  "<<mic+mes<<endl;         // save total zoopl in mmol N m-3
	out15<<(tt[k+1]+HSTEP*yy)/TIME<<"  "<<pco2w*1.0e6<<endl;     // save pCO2 in seawater 
	out16<<(tt[k+1]+HSTEP*yy)/TIME<<"  "<<co32<<endl;            // save [CO32-]  
	out17<<(tt[k+1]+HSTEP*yy)/TIME<<"  "<<o_cal<<endl;           // save omega-calcite
	out18<<(tt[k+1]+HSTEP*yy)/TIME<<"  "<<o_ara<<endl;           // save omega-aragonite
      }
      
      // save poincare' sections
      if(yy>IGNY && k==HOFY){
	outv<<y[1][k+1]+y[2][k+1]<<"  "<<y[5][k+1]<<endl;  // save Z-P
      }
      
      // new initial conditions
      if(k==nstep-3){
	vstart[1]=y[1][k+1];
	vstart[2]=y[2][k+1];
	vstart[3]=y[3][k+1];
	vstart[4]=y[4][k+1];
	vstart[5]=y[5][k+1];
	vstart[6]=y[6][k+1];
	vstart[7]=y[7][k+1];
	vstart[8]=y[8][k+1];
	vstart[9]=y[9][k+1];
	vstart[10]=y[10][k+1];
	vstart[11]=y[11][k+1];
	vstart[12]=y[12][k+1];
	vstart[13]=y[13][k+1];
	vstart[14]=y[14][k+1];

	yi1=vstart[1];
	yi2=vstart[2];
	yi3=vstart[3];
	yi4=vstart[4];
	yi5=vstart[5];
	yi6=vstart[6];
	yi7=vstart[7];
	yi8=vstart[8];
	yi9=vstart[9];
	yi10=vstart[10];
	yi11=vstart[11];
	yi12=vstart[12];
	yi13=vstart[13];
	yi14=vstart[14];
      }
      
      conta++;
      
      chlo=NTOC*(chlcd*y[1][k+1]+chlcdf*y[2][k+1]+chlcf*y[8][k+1]+chlceh*y[9][k+1]);

      sil=y[4][k+1];
      amm=y[10][k+1];
      tco2=y[13][k+1];   
      alk=y[14][k+1];

    }
    
    free_dvector(dv,1,nvar);
    free_dvector(vout,1,nvar);
    free_dvector(v,1,nvar);
  
  
  }
}

void rk4(double y[], double dydt[], int n, double t, double h, double yout[],
	 void (*derivs)(double, double [], double []))
{
  int i;
  double th, hh, h6;
  double *dym, *dyt, *yt;

  dym=dvector(1,n);
  dyt=dvector(1,n);
  yt=dvector(1,n);

  hh=h*0.5;
  h6=h/6.0;
  th=t+hh;

  for(i=1;i<=n;i++) yt[i]=y[i]+hh*dydt[i];  // first step
  (*derivs)(th,yt,dyt);                     // second step
  for(i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
  (*derivs)(th,yt,dym);                     // third step
  for(i=1;i<=n;i++){
    yt[i]=y[i]+h*dym[i];
    dym[i]+=dyt[i];
  }
  (*derivs)(t+h,yt,dyt);                    // fourth step 
  // accumulate increments with proper weights
  for(i=1;i<=n;i++) yout[i]=y[i]+h6*(dydt[i]+dyt[i]+2.0*dym[i]);

  free_dvector(yt,1,n); 
  free_dvector(dyt,1,n);
  free_dvector(dym,1,n);
}


//============================= DERIVS ROUTINE ================================


void derivs(double t, double y[], double dydt[])
{

  double phid, phif, phidf,phieh,phis;
  double phinr, phiar; 

  double sinkd=0.0;
  double sinko=0.0;

  // grazing terms
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  double g4=0.0;
  double g5=0.0;
  double g7=0.0;
  double g8=0.0;
  double g9=0.0;  

  double varHp;

  double calc=0.0;
  double cocpereh=0.0;
  double detach=0.0;

  double ad=0.0;
  double af=0.0;
  double adf=0.0;
  double aeh=0.0;

  varHp=max(varH,0.0);

  //cout<<varHp<<"  "<<varH<<endl;

  //========================= AMMONIA ==========================

  //phid=min(y[3]/(y[3]+NHD),y[10]/(y[10]+PHD));
  //phif=min(y[3]/(y[3]+NHF),y[10]/(y[10]+PHF));
  //phidf=min(y[3]/(y[3]+NHDF),y[10]/(y[10]+PHDF));
  //phieh=min(y[3]/(y[3]+NHEH),0.7+(0.3*y[10]/(y[10]+PHEH)));

  double qd1=0.0;
  double qd2=0.0;
  double qf1=0.0;
  double qf2=0.0;
  double qdf1=0.0;
  double qdf2=0.0;
  double qeh1=0.0;
  double qeh2=0.0;
  
  int i;
  for(i=1;i<=14;i++) y[i]=fabs(y[i]);

  qd1=(y[3]/NHD)/(1.0 + y[3]/NHD + y[10]/AHD);
  qd2=(y[10]/AHD)/(1.0 + y[3]/NHD + y[10]/AHD);

  qf1=(y[3]/NHF)/(1.0 + y[3]/NHF + y[10]/AHF);
  qf2=(y[10]/AHF)/(1.0 + y[3]/NHF + y[10]/AHF);
  
  qdf1=(y[3]/NHDF)/(1.0 + y[3]/NHDF + y[10]/AHDF);
  qdf2=(y[10]/AHDF)/(1.0 + y[3]/NHDF + y[10]/AHDF);
  
  //if(yy<Y-6){ 
  //  qeh1=0.0;
  //  qeh2=0.0;
  //}
  //else{
  qeh1=(y[3]/NHEH)/(1.0 + y[3]/NHEH + y[10]/AHEH);
  qeh2=(y[10]/AHEH)/(1.0 + y[3]/NHEH + y[10]/AHEH);
  //}

  phid = qd1 + qd2;    //(y[3]/NHD + y[10]/AHD)/(1 + y[3]/NHD + y[10]/AHD);
  phif = qf1 + qf2;    //(y[3]/NHF + y[10]/AHF)/(1 + y[3]/NHF + y[10]/AHF);
  phidf = qdf1 + qdf2; //(y[3]/NHDF + y[10]/AHDF)/(1 + y[3]/NHDF + y[10]/AHDF);
  phieh = qeh1 + qeh2; //(y[3]/NHEH + y[10]/AHEH)/(1 + y[3]/NHEH + y[10]/AHEH);

  //============================================================


  phis=y[4]/(y[4]+SH); 

  phid=min(phid,phis);

  // microzooplankton grazing 
  // g2: on flagellates
  // g5: on Ehuxleyi
  // g7: on diatoms
  // g8: on detritus
  //g2=ZMIF*P2*y[2]*y[2]*y[7]/(KMIG*(P2*y[2]+P5*y[9]+P7*y[1])+(P2*y[2]*y[2]+P5*y[9]*y[9]+P7*y[1]*y[1])); 
  //g5=ZMIE*P5*y[9]*y[9]*y[7]/(KMIG*(P2*y[2]+P5*y[9]+P7*y[1])+(P2*y[2]*y[2]+P5*y[9]*y[9]+P7*y[1]*y[1])); 
  //g7=ZMID*P7*y[1]*y[1]*y[7]/(KMIG*(P2*y[2]+P5*y[9]+P7*y[1])+(P2*y[2]*y[2]+P5*y[9]*y[9]+P7*y[1]*y[1]));
  //g7=0.0;

  ingDI=0.45/(tanh(y[4])+y[4]);

  if(yy<Y-6){ 
    g7=0.0;//ZMID*P7d*y[1]*y[1]*y[7]/(KMIG*(P2d*y[2]+P5d*y[9]+P7d*y[1])+(P2d*y[2]*y[2]+P5d*y[9]*y[9]+P7d*y[1]*y[1]));
    g5=ZMIE*P5d*y[9]*y[9]*y[7]/(KMIG*(P2d*y[2]+P5d*y[9]+P7d*y[1])+(P2d*y[2]*y[2]+P5d*y[9]*y[9]+P7d*y[1]*y[1]));
    g2=ZMIF*P2d*y[2]*y[2]*y[7]/(KMIG*(P2d*y[2]+P5d*y[9]+P7d*y[1])+(P2d*y[2]*y[2]+P5d*y[9]*y[9]+P7d*y[1]*y[1]));
  }  
  else{
    if(y[4]<3.0){ //If silicate is less than 3uM    
      g7=ZMID*P7*y[1]*y[1]*y[7]/(KMIG*(P2*y[2]+P5*y[9]+P7*y[1])+(P2*y[2]*y[2]+P5*y[9]*y[9]+P7*y[1]*y[1]));  
      g5=ZMIE*P5*y[9]*y[9]*y[7]/(KMIG*(P2*y[2]+P5*y[9]+P7*y[1])+(P2*y[2]*y[2]+P5*y[9]*y[9]+P7*y[1]*y[1])); 
      g2=ZMIF*P2*y[2]*y[2]*y[7]/(KMIG*(P2*y[2]+P5*y[9]+P7*y[1])+(P2*y[2]*y[2]+P5*y[9]*y[9]+P7*y[1]*y[1])); 
    }
    else{
      g7=0.0;//ZMID*P7d*y[1]*y[1]*y[7]/(KMIG*(P2d*y[2]+P5d*y[9]+P7d*y[1])+(P2d*y[2]*y[2]+P5d*y[9]*y[9]+P7d*y[1]*y[1]));
      g2=ZMIF*P2d*y[2]*y[2]*y[7]/(KMIG*(P2d*y[2]+P5d*y[9]+P7d*y[1])+(P2d*y[2]*y[2]+P5d*y[9]*y[9]+P7d*y[1]*y[1])); 
      g5=ZMIE*P5d*y[9]*y[9]*y[7]/(KMIG*(P2d*y[2]+P5d*y[9]+P7d*y[1])+(P2d*y[2]*y[2]+P5d*y[9]*y[9]+P7d*y[1]*y[1])); 
    }
  }

  // mesozooplankton grazing
  // g1: on diatoms
  // g3: on dinoflagellates
  // g4: on microzzoplankton
  // g9: on detritus
  g1=ZMED*P1*y[1]*y[1]*y[5]/(KMEG*(P1*y[1]+P3*y[8]+P4*y[7])+(P1*y[1]*y[1]+P3*y[8]*y[8]+P4*y[7]*y[7])); 
  g3=ZMEDF*P3*y[8]*y[8]*y[5]/(KMEG*(P1*y[1]+P3*y[8]+P4*y[7])+(P1*y[1]*y[1]+P3*y[8]*y[8]+P4*y[7]*y[7])); 
  g4=ZMEMI*P4*y[7]*y[7]*y[5]/(KMEG*(P1*y[1]+P3*y[8]+P4*y[7])+(P1*y[1]*y[1]+P3*y[8]*y[8]+P4*y[7]*y[7]));

  // diatoms sinking accelerates as silicate is depleted
  
  // Fasham style  (Rosa's upgrade)
  //sink=VD*(SH/(y[4] + SH));
       
  // Toby's style (TYRR96)
  // diatoms
  sinkd=VD;
  if(y[4]<2.0){
    sinkd=VD*(1.0+(7.0*(2.0-y[4])/2.0));
  }
  // others
  sinko=VDO;

  // Eslinger's style (ESLI01)
  //sink=5.0*(1.0-tanh(0.1375*y[4]));

  // Pondaven's style (POND99)
  //if(y[3]>NHD && y[4]>SH) sink=0.0/24.0;
  //if(y[3]<NHD || y[4]<SH) sink=5.0/24.0;
  if(yy<Y-6) calc=0.0;
  else calc=CALMAX*varT*psica; // CALMAX in: mmol cal-C (mmol org-C)-1 h-1 = mg cal-C (mg org-C)-1 h-1

  // transfer from attached (coccosphere) liths to free liths is
  // governed by the average number of liths usually found per cell
  // = ((y[11]/COCCAR) / (y[9]/EHOCAR)).
  // The ratio of liths in the coccosphere to E.hux cells should not
  // exceed this number.  It is assumed that there will always be
  // some (perhaps small) detachment of free liths
  // y[11]  must be in: mmol cal-C m-3 
  // y[9]   must be in: mmol org-C m-3
  // COCCAR must be in: mmol cal-C coccolith-1
  // EHOCAR must be in: mmol org-C cell-1
  
   if(yy<Y-6){
     cocpereh=0.0;
     detach=0.0;
   }
   else{
     //cocpereh=(y[11]/(CTON*y[9]))/(COCCAR/EHOCAR);
     //detach=((cocpereh-COCMAX)*((CTON*y[9])/EHOCAR)*COCCAR)+(MEH*y[11]);
     //if(detach<(DETMIN*y[11])) detach=DETMIN*y[11];
     detach=max(DET*(y[11]-(COCMAX*COCCAR*(CTON*y[9]/EHOCAR))), (DETMIN*y[11]));
   }
  // detachment of coccoliths (number of cocco. getting detached from the cell as in Eq. 10 in TYRR96)
  // in mmol calcite C m-3 day-1
  // using CTON* instead of ICOC*
  // double detach=max((y[11]-(COCMAX*COCCAR*(6.0/30.0)*(CTON*y[9]/EHOCAR))),(DETMIN*y[11]));

  // growth terms
  ad = MUD0*varT*psi*phid;      // diatoms
  af = MUF0*varT*psi*phif;      // flagellates
  adf = MUDF0*varT*psi*phidf;   // dinoflagellates

  //if(yy<Y-6) aeh = 0.0;
  //else 
  aeh = MUEH0*varTeh*psieh*phieh; // Ehuxleyi


  // -- [1] -- ODE FOR DIATOMS -- in: mmol N m-3
  
  dydt[1] = ad*y[1] - g1 - g7 - MD*y[1] - ((sinkd+diff+varHp)/mixed)*y[1]; 


  // -- [2] -- ODE FOR FLAGELLATES -- in: mmol N m-3
  
  dydt[2] = af*y[2] - g2 - MF*y[2] - ((sinko+diff+varHp)/mixed)*y[2];  

  
  // -- [3] -- ODE FOR NITRATE -- in: mmol N m-3
  
  dydt[3] = - MUD0*varT*psi*(qd1/(qd1+qd2))*phid*y[1] - MUF0*varT*psi*qf1*y[2] - MUDF0*varT*psi*qdf1*y[8] - 
              MUEH0*varTeh*psieh*qeh1*y[9] + NIT*y[10] + ((diff+varHp)/mixed)*(nbo-y[3]); 


  // -- [4] -- ODE FOR SILICATE -- in: mmol Si m-3

  dydt[4] = - ad*y[1] + ((diff+varHp)/mixed)*(sbo-y[4]);  

  
  // -- [5] -- ODE FOR MESOZOOPLANKTON -- in: mmol N m-3 (graze on: diatom, dinofla, microzoo, detritus)

  dydt[5] = B1*g1 + B3*g3 + B4*g4 + B9*g9 - EXME*y[5] - MZME*y[5]*y[5] - (varH/mixed)*y[5];    


  // -- [6] -- ODE FOR DETRITUS -- in: mmol N m-3

  dydt[6] = (1-B1)*g1 + (1-B2)*g2 + (1-B3)*g3 + (1-B4)*g4 + (1-B5)*g5 + (1-B7)*g7 + (1-B8)*g8 + (1-B9)*g9 +
            MD*y[1] + MF*y[2] + MDF*y[8] + MEH*y[9] - g8 - g9 - MDE*y[6] - ((diff+varHp+VDT)/mixed)*y[6];  
           

  // -- [7] -- ODE FOR MICROZOOPLANKTON -- in: mmol N m-3 (graze on: flage, Ehux, free cocco, detritus) 

  dydt[7] = B2*g2 + B5*g5 + B7*g7 + B8*g8 - EXMI*y[7] - MZMI*y[7]*y[7] - g4 - (varH/mixed)*y[7];


  // -- [8] -- ODE FOR DINOFLAGELLATES -- in: mmol N m-3
  
  dydt[8] = adf*y[8] - g3 - MDF*y[8] - ((sinko+diff+varHp)/mixed)*y[8];  


  // -- [9] -- ODE FOR EMILIANIA HUXLEYI -- in: mmol N m-3 here, in output file also in mmol C m-3

  //if(yy<Y-6) dydt[9] = 0.0;
  //else 
  dydt[9] = aeh*y[9] - g5 - MEH*y[9] - ((sinko+diff+varHp)/mixed)*y[9];  

  // -- [10] -- ODE FOR AMMONIUM -- in: mmol N m-3

  dydt[10] = - MUD0*varT*psi*(qd2/(qd1+qd2))*phid*y[1] - MUF0*varT*psi*qf2*y[2] - 
               MUDF0*varT*psi*qdf2*y[8] - MUEH0*varTeh*psieh*qeh2*y[9] +
               (EXME*y[5] + EXMI*y[7] + FZRME*MZME*y[5]*y[5] + FZRMI*MZMI*y[7]*y[7] + MDE*y[6]) - 
               NIT*y[10] - ((diff+varHp)/mixed)*y[10]; 


  // -- [11] -- ODE FOR ATTACHED COCCOLITHS -- in: mmol calcite-C m-3 here
  //
  // Attached coccoliths: calcification (i.e. newly produced coccoliths, attached) - grazing - 
  //                      cell mortality - detachment - mixing
  if(yy<Y-6) dydt[11] = 0.0;
  else dydt[11] = calc*CTON*y[9] - (g5/y[9])*y[11] - MEH*y[11] - detach - ((diff+varHp)/mixed)*y[11]; 


  // -- [12] -- ODE FOR FREE COCCOLITHS -- in: mmol calcite-C m-3 here
  //
  // Free coccoliths: detached cocc. + empty coccosphere (when cells die) + 
  //                  fraction of cocco not ingested during grazing -
  //                  grazing on free coccoliths - dissolution - mixing 

  if(yy<Y-6) dydt[12] = 0.0;
  else dydt[12] = detach + MEH*y[11] + 0.1*(g5/y[9])*y[11] - DISSOL*y[12] - ((diff+varHp)/mixed)*y[12];
  //0.5*(g5/y[9])*y[12]

  // -- [13] -- ODE FOR DISSOLVED INORGANIC CARBON -- in: umol C m-3
  //
  dydt[13] = - CTON*(ad*y[1] + af*y[2] + adf*y[8] + aeh*y[9] + calc*y[9]) + CTON*MDE*y[6] + 
               CTON*(EXME*y[5] + EXMI*y[7] + FZRMI*MZMI*y[7]*y[7] + FZRME*MZME*y[5]*y[5]) + 
               DISSOL*y[12] + gtv*co2sol*(PCO2A-pco2w)/mixed + ((diff+varHp)/mixed)*(DIC0-y[13]); 
  
  //                                NOTE:
  // 
  //                           co2sol in: umol kg-1 atm-1
  //                            PCO2A in: atm
  //                            pco2w in: atm
  //                              gtv in: m h-1
  //                            mixed in: m
  // GTV*co2sol*(PCO2A - pco2w)/mixed in: umol kg-1 h-1


  // -- [14] -- ODE FOR TOTAL ALKALINITY -- in: uEq m-3
  //
  // Total Alkalinity:  + 2 x dissolution - 2 x calcification                    (CO32- carries two negative charges)
  //                    + nitrate upake - ammonium uptake                        (NO3- is negative and NH4+ is positive)  
  //                    + excretions + mortality fractions + detritus breakdown  (all three terms increase NH4+ pool)
  //
  //
  // NOTE EMAIL FROM TOBY: dAlk/dt = -2*calcification + 2*dissolution + diffusion*(Alk_depth - Alk)
  //
  // ALL THE REST: nitrate upake, ammonium uptake, ammonification, etc. is negligible!
  //
  dydt[14] = - 2.0*calc*CTON*y[9] + 2.0*DISSOL*y[12] + ((diff+varHp)/mixed)*(ALK0-y[14]);
    //           MUD0*varT*psi*(qd1/(qd1+qd2))*phid*y[1] + MUF0*varT*psi*qf1*y[2] + 
    //           MUDF0*varT*psi*qdf1*y[8] + MUEH0*varTeh*psieh*qeh1*y[9] -
    //           (MUD0*varT*psi*(qd2/(qd1+qd2))*phid*y[1] + MUF0*varT*psi*qf2*y[2] + 
    //           MUDF0*varT*psi*qdf2*y[8] + MUEH0*varTeh*psieh*qeh2*y[9] + NIT*y[10]) +
    //           (EXMI*y[7] + FZRMI*MZMI*y[7]*y[7] + EXME*y[5] + FZRME*MZME*y[5]*y[5] + MDE*y[6]) + 
    //           ((diff+varHp)/mixed)*(ALK0-y[14]);


  // =========== DIAGNOSTIC VARIABLES ============

  grazd=g7/y[1];         // microzoo grazing on diatoms
  if(yy<Y-6) graze=0.0;  // microzoo grazing on Ehux before 1995
  else graze=g5/y[9];    // microzoo grazing on Ehux after 1995

  regphypro = MUD0*varT*psi*(qd2/(qd1+qd2))*phid*y[1] + MUF0*varT*psi*qf2*y[2] + 
              MUDF0*varT*psi*qdf2*y[8] + MUEH0*varT*psieh*qeh2*y[9]; 

  regdiapro = (y[10]/AHF)/(1.0 + y[3]/NHF + y[10]/AHF);//MUD0*varT*psi*(qd2/(qd1+qd2))*phid*y[1];
  regdinpro = y[10]/AHF;//MUDF0*varT*psi*qdf2*y[8];
  regflapro = 1.0 + y[3]/NHF + y[10]/AHF;//MUF0*varT*psi*qf2*y[2];
  regehupro = y[3];//MUEH0*varT*psieh*qeh2*y[9];
  regtest = y[10];

  newphypro = MUD0*varT*psi*(qd1/(qd1+qd2))*phid*y[1] + MUF0*varT*psi*qf1*y[2] + 
              MUDF0*varT*psi*qdf1*y[8] + MUEH0*varT*psieh*qeh1*y[9]; 

  totphypro = ad*y[1] + af*y[2] + adf*y[8] + aeh*y[9];

  totphyloss = g1+g7+MD*y[1]+((sinkd+diff+varHp)/mixed)*y[1] + g2+MF*y[2]+((sinko+diff+varHp)/mixed)*y[2] +
               g3+MDF*y[8]+((sinko+diff+varHp)/mixed)*y[8]+g5 + MEH*y[9]+((sinko+diff+varHp)/mixed)*y[9];

  totphymix = ((sinkd+diff+varHp)/mixed)*y[1]+((sinko+diff+varHp)/mixed)*y[2]+
              ((sinko+diff+varHp)/mixed)*y[8]+((sinko+diff+varHp)/mixed)*y[9];

  pon = y[1]+y[2]+y[8]+y[9]+y[7]+y[5]+y[6]; // phy + zoo + det 

  calcieh = calc*CTON*y[9];   // PIC in: mmol inorganic C m-3 h-1

  photoeh = aeh*CTON*y[9];    // POC in: mmol organic C m-3 h-1

  totzoopro = B1*g1 + B2*g2 + B3*g3 + B4*g4 + B5*g5 + B7*g7;

  totzooloss = EXME*y[5]+MZME*y[5]*y[5]+(varH/mixed)*y[5] + EXMI*y[7]+MZMI*y[7]*y[7]+g4+(varH/mixed)*y[7];

  dianutgro=MUD0*varT*phid;
  dinnutgro=MUDF0*varT*phidf;
  flanutgro=MUF0*varT*phif;
  ehunutgro=MUEH0*varT*phieh;

  dialightgro=MUD0*varT*psi;
  dinlightgro=MUDF0*varT*psi;
  flalightgro=MUF0*varT*psi;
  ehulightgro=MUEH0*varT*psieh;

  diagra=(g1+g7)/y[1];
  dingra=g3/y[8];
  flagra=g2/y[2];
  ehugra=g5/y[9];
  micgra=g4/y[7];

  callightgro=CTON*calc;
  caltemgro=(detach+MEH*y[11]+0.1*(g5/y[9])*y[11]);

  // ============ MASS BALANCE CHECK =============


  //double mbs=0.0;    
  //double mbn=0.0;
  //double diff=0.0;

  //diff=(mm+varHp)/mixed;

  // -- silicate cycle --
  //dydt[7] = 0.5*ZDG*y[5]*y[1] + 0.8*MDMAX*y[1]*y[1]/(KD+y[1]) + (VD/mixed)*y[1] + diff*y[1] - diff*(S0-y[4]);
  //mbs = yi1 + yi4 - (y[1] + y[4] + y[7]);  // silicate mass balance

  // -- nitrogen cycle --
  //dydt[8] = (VD/mixed)*y[1] + (VDF/mixed)*y[2] + (VDT/mixed)*y[6] + 
  //          diff*y[1] + diff*y[2] + (varH/mixed)*y[5] + diff*y[6] - diff*(N0-y[3]);

  //mbn = yi1 + yi2 + yi3 + yi5 + yi6 - (y[1] + y[2] + y[3] + y[5] + y[6] + y[8]); // nitrogen mass balance
  
  //if(mbs>1.0e-12) cout<<" - warning! - mass balance violation in the silicate cycle "<<"\n";
  //if(mbn>1.0e-12) cout<<" - warning! - mass balance violation in the nitrogen cycle "<<"\n";          


  // =============================================


}










