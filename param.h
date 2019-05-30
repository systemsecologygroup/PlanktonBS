// 
//                      param.h
//
//                    header file
//
//                  Agostino Merico
//
//
//  Definition of parameters and function declarations
//

#define PI 3.141592653589793 


// ===== MORTALITY, EXCRETION AND RESPIRATION =====

//#define MMAX 0.2    // max mortality in d-1, as in TAYL91
//#define MMIN 0.1    // min mortality in d-1, as in TAYL91

#define MD 0.04/24.0 //0.045 in pre-1995 run 0.08 diato max specific mortality rate (0.05 d-1, as in FASH93)
#define MF 0.04/24.0 //0.045 in pre-1995 run 0.08 flage max specific mortality rate (0.05 d-1, as in FASH93)
#define MDF 0.04/24.0//0.045 in pre-1995 run 0.08 dinof max specific mortality rate (0.05 d-1, as in FASH93)
#define MEH 0.04/24.0//0.045 in pre-1995 run 0.08 ehuxl max specific mortality rate (0.05 d-1, as in FASH95)

#define MZMI 0.05/24.0// (late eh 0.05) (sta. 0.05) micr max specific mort. rate (0.3 d-1 as in FASH95)
#define MZME 0.2/24.0 // (late eh 0.13) (sta. 0.2)  meso max specific mort. rate (0.3 d-1 as in FASH95)

#define EXMI 0.025/24.0// (late eh 0.028) (sta. 0.025) microzooplankton excretion rate (in d-1)
#define EXME 0.1/24.0  // (late eh 0.130) (sta. 0.1)   mesozooplankton  excretion rate (in d-1)

#define MDE 0.1/24.0  // 0.05 detrital breakdown rate in h-1 (0.05 d-1, 5%), as in FASH95a  

#define KZMI 0.2       // 0.2 microzooplankton mortality half saturation constant (in mmol/m3 as in FASH93)
#define KZME 0.2       // 0.2 mesozooplankton  mortality half saturation constant (in mmol/m3 as in FASH93)


// ==================== GRAZING ====================

// in [(mg C m-3)-1 d-1]
// 1 mmol N m-3 = 79.2 mg C m-3

// Fasham style (as in FASH93)

#define B1 0.75       // mesozooplankton assimilation efficiency for diatoms
#define B3 0.75       // mesozooplankton assimilation efficiency for dinofla
#define B4 0.75       // mesozooplankton assimilation efficiency for microzo
#define B9 0.0        // mesozooplankton assimilation efficiency for detritu

#define B2 0.75       // microzooplankton assimilation efficiency for flagell
#define B5 0.75       // microzooplankton assimilation efficiency for ehuxley
#define B7 0.75       // microzooplankton assimilation efficiency for diatoms
#define B8 0.0        // microzooplankton assimilation efficiency for detritu

#define P1 0.33       // mesozooplankton feeding preference for diatoms
#define P3 0.33       // mesozooplankton feeding preference for dinofla 
#define P4 0.33       // mesozooplankton feeding preference for microzo

#define P2 0.5        // microzooplankton feeding preference for flag - SIL < 3uM
#define P5 0.33       // microzooplankton feeding preference for ehux - SIL < 3uM
#define P7 0.5        // microzooplankton feeding preference for diat - SIL < 3uM

#define P2d 0.5       // microzooplankton feeding preference for flag - SIL > 3uM
#define P5d 0.33      // microzooplankton feeding preference for ehux - SIL > 3uM
#define P7d 0.5       // microzooplankton feeding preference for diat - SIL > 3uM

#define KMIG 1.0      // microzoopl feeding half-saturation constant (in mmol/m3)
#define KMEG 1.0      // mesozoopla feeding half-saturation constant (in mmol/m3)

#define ZMID 0.175/24.0   // SR 0.7 or 0.0 microzoo maximum ingestion rate on diat (in d-1)
#define ZMIF 0.7/24.0   // SR 0.7 or 0.7 microzoo maximum ingestion rate on flag (in d-1)
#define ZMIE 0.175/24.0 // SR 0.175 or 0.7 microzoo maximum ingestion rate on ehux (in d-1)

#define ZMIAC 0.3/24.0 // 0.3 or 0.7 microzoo maximum ingestion rate on attach. coccoliths 

#define ZMED 0.7/24.0  // 0.7 mesozoo maximum ingestion rate on diatom (in d-1)
#define ZMEDF 0.7/24.0 // 0.7 mesozoo maximum ingestion rate on dinofl (in d-1)
#define ZMEMI 0.7/24.0 // 0.7 mesozoo maximum ingestion rate on microz (in d-1)

#define FZRME 0.1      // Fraction of dead mesozooplankton recycled into ammonium within the mixed layer
                       // the remaining is assumend to sink and exported rapidly out of the MLD
#define FZRMI 0.1      // Fraction of dead microzooplankton recycled into ammonium within the mixed layer
                       // the remaining is assumend to sink and exported rapidly out of the MLD

#define FZRNOEH 0.2


#define ED 0.5
#define EF 0.5
#define EDE 0.5
#define GGD 1.0/24.0
#define GGF 1.0/24.0
#define GGDE 1.0/24.0


// =================== COCCOLITHS ==================

#define CALMAX 0.2/24.0 // max calcification rate in mg cal-C (mg org-C)-1 day-1 as in Fernandez et al. 1993
#define DISSOL 0.08/24.0// dissolution rate of calcite (CaCO3) in within the MLD in day-1 (0.05 d-1 in TYRR96)
#define COCMAX 15.0     // max num of coccoliths attached to a cell in: coccoliths cell-1 as THETA in TYRR96
#define COCCAR 0.25e-12 // inogranic carbon content of a coccolith in: g calcite-C cell-1 as C_lith in TYRR96
#define EHOCAR 10.0e-12 // organic carbon content of an Ehux cell in: g organic-C cell-1 as C_eh in TYRR96
#define DETMIN 0.1/24.0 // minimum rate of coccoliths detachment (in day-1 as in TYRR96)
#define DET 23.5/24.0   // detachment rate (in day-1) 
#define GAL 0.2         // linear lose of attached coccoliths due to grazing 

#define ICOC 0.75      // calcite-C/coccolith to organic-C/cell ratio (30 coccoliths/cell are assumed)
                       // (calcite-C * 30/cell)* (cell/organic-C) = 30 * calcite-C/organic-C

// ================== TEMPERATURE ==================

#define TMAX 14.0     // max temperature in degree C, as in WOD98 for the Bering Sea
#define TMIN 1.0      // min temperature in degree C, as in WOD98 for the Bering Sea 


// ===================== LIGHT =====================

#define ISAT 100.0   //120 or SR 100 in W m-2 Light intensity that gives max photosynth. rate for all phyto. 
                     //(500 in uEin m-2 s-1 as in KIRK pag 274)

#define ISATEH 280.0 //SR 280 in W m-2 Light intensity that gives max photosynth. rate for Ehux. 
                     //(1500 in uEin m-2 s-1 in NANN96 pag 199)

#define KRE 0.3     // red light extinction coeff. in m-1, (as in TAYL93a, table I)
#define KGR 0.058   // green light extinction coeff. in m-1, (as in TAYL93a, table I)
#define KSS 0.03    // (standard value 0.03) phytoplankton self-shading coefficient (in m2 mmol-1)
#define KW 0.04     // extincion coefficient due to water in m-1 (FASH90)

#define IHD 15.0         //(9.60 W m-2) All phy. light half-sat. const. (83.4 uEin m-2 s-1, not as in TYRR96)
#define IHEH 45.0 //50.0 //(24.0 W m-2) Ehuxleyi light half-sat. const. (83.4 uEin m-2 s-1, not as in TYRR96)

#define IHCA 40.0// was 15.0 for resSTANDARD.dat //(9.60 W m-2) Calcification light half-sat. const. (40.0 uEin m-2 s-1, as in TYRR96)

#define G 1.0       // pigment concentration (1-15 mg m-3 as in ANDE93)


// ============= MIXED LAYER DEPTH =================

#define H 20.0      // case with constant mixed layer depth (in m)


// ======== MAX GROWTH RATES AT 0 degree C =========

#define MUD0 1.2/24.0     // BL 1.35 or SR 1.20 diatoms maximum growth rate in d-1 (1.287 in AKSN94)
#define MUF0 0.65/24.0    // BL 0.70 or SR 0.65 flagell maximum growth rate in d-1 (0.852 in AKSN94)
#define MUDF0 0.6/24.0    // BL 0.60 or SR 0.60 dinofla maximum growth rate in d-1
#define MUEH0 1.15/24.0   // BL 0.70 or SR 1.15 Ehuxley maximum growth rate in d-1 (0.852 in AKSN94)

// ==== Case without Ehuxleyi ====

#define MUD0_NOEH 1.2/24.0  // 1.3 diato maximum growth rate in d-1 (1.287 in AKSN94)
#define MUF0_NOEH 0.5/24.0  // 0.7 flage maximum growth rate in d-1 (0.852 in AKSN94)
#define MUDF0_NOEH 0.6/24.0 // 0.4 dinof maximum growth rate in d-1


//=============================================================
//========================== UNITS ============================
//
//
// note: 1 uM = 1 umol/L = 1 umol/dm3 = 1 mmol/m3
//
// note:        1 mmol/L = 1 umol/kg   
// (with 1.025 = seawater density, here assumed equal to 1)
//
// note: 1 ug-atom/L     = 1 uM       = 1 mg-atoms/m3
// 
// 
// note:         1 mg/m3 = 1 ug/L 
//
//=============================================================
//=============================================================


// ========== REDFIELD AND OTHER RATIOS ==========

// Cellular elemental composition (mol-ratio) as in AKSN94

#define R 0.063      // P:N ratio (1:16) 
#define RD 0.063     // P:N ratio in diatoms (1:16) 
#define RF 0.063     // P:N ratio in flagellates (1:16)
#define RDF 0.063    // P:N ratio in dinoflagel. (1:16)
#define REH 0.031    // P:N ratio in Ehuxleyi (0.5:16)

#define NTOC 79.5      // phytoplankton transform. factor from (mmol N m-3) to (mg C m-3) using C:N=6.625
#define NTOCZ 67.5     // zooplankton transform.   factor from (mmol N m-3) to (mg C m-3) using C:N=5.625
#define CTOCHL 50.0    // C:Chl constant ratio
#define CHLTOC 1.0/50.0// Chl:C constant ratio

#define CTONZ 5.625     // C:N ratio in zooplankton
#define CTON 6.625      // C:N ratio (106:16)


// ============== CARBONATE SYSTEM ================

#define PCO2A 358.0*1.0e-6// Constant pCO2 in air in atm as from MURA02, (345 is a value between PENG87 and TAYL91)
#define DIC0 2000.0       // 2100 Total CO2  (DIC) below the MLD in mmol/m3 (= umol/kg) (as in WALS94)
#define ALK0 2100.0       // Total alkalinity below the MLD in ueq kg-1 (as in WALS94)

//#define RE 0.2       // recycling efficiency as in TYRR95
//#define NRE 0.1      // Nutrients recycling efficiency as in TYRR95
//#define FDC 0.05     // fraction of detritus going into the DIC pool
//#define FCC 0.15     // fraction of free coccoliths going into the DIC pool
//#define FCA 0.15     // fraction of free coccoliths going into the TAlk pool 

// C:N Redfield ratio - conversion from (mmol N m-3) of seawater to (umol C kg-1) 
#define CRI 6.625/1.025 // = C/N/(seawater density) reduction/increase of [C] in 
                        //   water when 1 (mmol N m-3) is produced (TYRR95) 

#define WD 1.0/1.025 // 1.025 is water density, conv. factor from (mmol N m-3) of seawater to (uEq kg-1)
#define GTV 4.2/24.0 // constant gas transfer velocity (or piston velocity) in m/d (as in TYRR95)


// ================= NUTRIENTS ====================

#define NIT 0.05/24.0 // nitrification rate (ammonium becoming nitrate) 5% day-1

#define A0 0.0      // ammonium conc. below mixed layer in mmol/m3 as in WALS94

#define N0 20.0     // 20 peng87=13.0     // nitrogen conc. below mixed layer in mmol/m3

#define N094 20.0
#define N095 20.0 
#define N096 20.0
#define N097 20.0
#define N098 20.0
#define N099 20.0 // was 20 (high N) not good for ehux in 2000
#define N000 20.0
#define N001 20.0

#define S0 35.0     // 35 peng87=7.0      // silicate conc. below mixed layer in mmol/m3

#define S094 35.0
#define S095 35.0  // 95 sil was 20 (low) but nit was 20 (high), which gave a big ehux in 96 (less silicate around helpd?) 
#define S096 35.0 
#define S097 35.0 
#define S097 35.0
#define S098 35.0 
#define S099 35.0 // was 40.0
#define S000 35.0
#define S001 35.0

#define D0 0.0       // detritus conc. below mixed layer in mmol/m3 as in eh60.dat

#define NHD 1.5      // SR 1.5 diatoms nitrate half sat. constant in mmol/m3 as in ESLI01
#define NHF 1.5      // SR 1.5 flagell nitrate half sat. constant in mmol/m3 as in ESLI01
#define NHDF 1.5     // SR 1.5 dinofla nitrate half sat. constant in mmol/m3 as in ESLI01
#define NHEH 1.5     // SR 1.5 ehuxley nitrate half sat. constant in mmol/m3 as in ESLI01

#define AHD 0.05     // SR 0.05 diatoms ammonium half sat. constant in mmol/m3 as in TYRR96
#define AHF 0.05     // SR 0.05 flagell ammonium half sat. constant in mmol/m3 as in TYRR96
#define AHDF 0.05    // SR 0.05 dinofla ammonium half sat. constant in mmol/m3 as in TYRR96
#define AHEH 0.05    // SR 0.05 ehuxley ammonium half sat. constant in mmol/m3 as in TYRR96

#define SH 3.5       // SR 3.5 silicate half saturation constant in mmol/m3 as in TYRR96

//#define PHD 0.25     // diatoms phosphate half sat. constant in mmol/m3 as in TYRR96
//#define PHF 0.25     // flagell phosphate half sat. constant in mmol/m3 as in TYRR96
//#define PHDF 0.25    // dinofla phosphate half sat. constant in mmol/m3 as in TYRR96
//#define PHEH 0.025   // ehuxley phosphate half sat. constant in mmol/m3 as in TYRR96



// ==================== SINKING ======================

#define VD 0.5/24.0    // sinking speed for diatoms (minimum value) in m/h
#define VDO 0.0001/24.0// sinking speed for phytoplankton other than diatoms in m/h

#define VDT 1.0/24.0   // 0.4 sinking speed for detritus (1.0 m/d in FASH90) THIS CONTROL THE EXPORT PRODUCTION

// ===================== MIXING ======================

//#define K 0.6        // vertical mixing coefficient in m/day as in TAYL91
//#define K 0.025      // vertical mixing coefficient in m/h

#define mm 0.01/24.0              // (0.01 is the SR) cross-thermocline mixing rate in m/h as in FASH93

#define mm95 0.01/24.0        
#define mm96 0.01/24.0       
#define mm97 0.01/24.0
#define mm98 0.01/24.0
#define mm99 0.01/24.0
#define mm00 0.01/24.0
#define mm01 0.01/24.0


#define TSTEP_FLAG 1 // flag for time step: 1 is day, 2 is hour, 3 is 1/2 hour




// ============= function prototypes =================


double get_light_at_surface(double);                // Calculate light at surface with astronomicae formulas
double get_light(double,double,double);             // Calculate averaged light in the MLD          
double get_averaged_light(double,double,double);    // Calculate limiting averaged light in the MLD
double get_averaged_light_eh(double,double,double); // As above but for E. huxleyi only
double get_averaged_light_cal(double,double,double); // As above but for Calcification only
double get_light_intensity(double,double);          // Calculate light at a given depth

double get_gas_transfer_velocity(double,double);       
double get_co2_solubility(double,double);           
double get_param_water(double,double,double,double,double,int); 

double min(double, double);
double max(double, double);


// ===========================================

