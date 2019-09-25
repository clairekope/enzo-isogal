/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GALAXY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, Feb, 2004
/  modified1:  Elizabeth Tasker, Oct, 2006 (tidied up)
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EnzoTiming.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define Mpc (3.0856e24)         //Mpc [cm] 
#define SolarMass (1.989e33)    //Solar Mass [g]
#define GravConst (6.67e-8)     //Gravitational Constant [cm3g-1s-2]
#define pi (3.14159)
#define mh (1.67e-24)           //Mass of Hydrogen [g]
#define kboltz (1.381e-16)      //Boltzmann's Constant [ergK-1]
#define kboltzKeV (8.617e-8)
#define mu (0.6)
#define CM_PER_KM (1.0e5)
#define CM_PER_KPC (3.0856e21)

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
                      float *TemperatureUnits, float *TimeUnits,
                      float *VelocityUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/* Internal routines */
void setup_chem(float density, float temperature, int equilibrate,
		float& DEdest, float&  HIdest, float& HIIdest,
		float& HeIdest, float& HeIIdest, float& HeIIIdest,
		float& HMdest, float& H2Idest, float& H2IIdest,
		float& DIdest, float& DIIdest, float& HDIdest);
float gasvel(FLOAT radius, float DiskDensity, FLOAT ExpansionFactor, 
             float GalaxyMass, FLOAT ScaleHeightR, FLOAT ScaleHeightz, 
             float DMConcentration, FLOAT Time);
float gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, 
                 FLOAT inv [3][3], float DiskDensity,
                 FLOAT ScaleHeightR, FLOAT ScaleHeightz,
                 FLOAT cellwidth);
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot,
                 FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3]);
double bilinear_interp(double x, double y, 
                       double x1, double x2, double y1, double y2,
                       double f_x1y1, double f_x1y2, 
                       double f_x2y1, double f_x2y2);

/* Internal Routines for CGM setup */
float HaloGasDensity(FLOAT R);
float HaloGasTemperature(FLOAT R);

/* Internal Routines for Disk Potential Setup */
float DiskPotentialCircularVelocity(FLOAT cellwidth,
                                    FLOAT z,FLOAT density,
                                    FLOAT &temperature);
double trapzd(double (func)(), double a, double b, int n);
double qromb(double (*func)(double), double a, double b);
void polint(double xa[],double ya[],int n,double x,double *y,double *dy);
static double drcyl;
static double r2;

static float DensityUnits, LengthUnits, TemperatureUnits = 1,
             TimeUnits, VelocityUnits, MassUnits;

double gScaleHeightR, gScaleHeightz, densicm, MgasScale, Picm,
       TruncRadius, SmoothRadius, SmoothLength,Ticm;

/* Global variables (within this file) for circumgalactic medium setup 
   (also used a bit for disk potential setup) */
int GalaxySimulationGasHalo, EquilibrateChem;
double GalaxySimulationGasHaloScaleRadius,
  GalaxySimulationGasHaloDensity, GalaxySimulationGasHaloDensity2,
  GalaxySimulationGasHaloTemperature, GalaxySimulationGasHaloAlpha,
  GalaxySimulationGasHaloZeta, GalaxySimulationGasHaloZeta2,
  GalaxySimulationGasHaloCoreEntropy, GalaxySimulationGasHaloGalaxyMass,
  GalaxySimulationGasHaloDMConcentration,
  GalaxySimulationGasHaloMetallicity,
  GalaxySimulationDiskMetallicityEnhancementFactor;

/* struct to carry around data required for circumgalactic media
   if we need to generate radial profiles of halo quantities via 
   numerical integration */
struct CGMdata {
  double *n_rad, *T_rad, *rad, R_outer, dr;
  int nbins;
};
struct CGMdata CGM_data;

/* declarations for a bunch of functions needed to generate radial profiles
   of halo quantities via numerical integration - see the actual functions for
   descriptions. */
double halo_S_of_r(double r, double n, grid* Grid);
double halo_dSdr(double r);
double halo_dn_dr(double r, double n, grid* Grid);
double halo_g_of_r(double r);
double halo_galmass_at_r(double r);
void halo_init(grid* Grid);
void halo_clean(void);

int grid::GalaxySimulationInitializeGrid(FLOAT DiskRadius,
           FLOAT GalaxyMass,
           FLOAT GasMass,
           FLOAT DiskPosition[MAX_DIMENSION], 
           FLOAT ScaleHeightz,
           FLOAT ScaleHeightR,
           FLOAT GalaxyTruncationRadius, 
           FLOAT DMConcentration,
           FLOAT DiskTemperature,
           FLOAT InitialTemperature,
           FLOAT UniformDensity,
           int   EquilChem,
           int   GasHalo,
           FLOAT GasHaloScaleRadius,
           FLOAT GasHaloDensity,
           FLOAT GasHaloDensity2,
           FLOAT GasHaloTemperature,
           FLOAT GasHaloAlpha,
           FLOAT GasHaloZeta,
           FLOAT GasHaloZeta2,
           FLOAT GasHaloCoreEntropy,
           FLOAT GasHaloMetallicity,
           int   UseHaloRotation,
           FLOAT RotationScaleVelocity,
           FLOAT RotationScaleRadius,
           FLOAT RotationPowerLawIndex,
           FLOAT DiskMetallicityEnhancementFactor,
           FLOAT AngularMomentum[MAX_DIMENSION],
           FLOAT UniformVelocity[MAX_DIMENSION], 
           int UseMetallicityField, 
           FLOAT GalaxySimulationInflowTime,
           FLOAT GalaxySimulationInflowDensity,
           int level,
           FLOAT GalaxySimulationInitialBfield[MAX_DIMENSION],
           int GalaxySimulationInitialBfieldTopology,
           FLOAT GalaxySimulationCR
          )
{
 /* declarations */

int dim, i, j, k, m, field, disk, size, MetalNum, MetalIaNum, vel;
int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
    H2IINum, DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
float DiskDensity, DiskVelocityMag;
int CRNum, DensNum;

/* global-scope variables for disk potential functions (would be better if not global) */

gScaleHeightR = ScaleHeightR;
gScaleHeightz = ScaleHeightz;
densicm = UniformDensity;  // gas density if no halo is used
MgasScale = GasMass;
Ticm = InitialTemperature;  // gas temperature if no halo is used
Picm = kboltz*UniformDensity*Ticm/(mu*mh);  // gas pressure if no halo is used
TruncRadius = GalaxyTruncationRadius;
SmoothRadius = TruncRadius*.02/.026;
SmoothLength = TruncRadius - SmoothRadius;

/* set all of the gas halo quantities to variables that are global
   within this file, so the calls to HaloGasDensity() and HaloGasTemperature()
   are unchanged. */
EquilibrateChem = EquilChem;
GalaxySimulationGasHalo = GasHalo;   // integer, >= 0
GalaxySimulationGasHaloScaleRadius = GasHaloScaleRadius;  // in mpc
GalaxySimulationGasHaloDensity = GasHaloDensity; // in grams/cm^3
GalaxySimulationGasHaloDensity2 = GasHaloDensity2;
GalaxySimulationGasHaloTemperature = GasHaloTemperature;  // in Kelvin
GalaxySimulationGasHaloAlpha = GasHaloAlpha;  // power-law index; unitless
GalaxySimulationGasHaloZeta = GasHaloZeta;
GalaxySimulationGasHaloZeta2 = GasHaloZeta2;
GalaxySimulationGasHaloCoreEntropy = GasHaloCoreEntropy;  // power-law index; unitless
GalaxySimulationGasHaloGalaxyMass = GalaxyMass;
GalaxySimulationGasHaloDMConcentration = DMConcentration;
GalaxySimulationGasHaloMetallicity = GasHaloMetallicity; // Zsun
GalaxySimulationDiskMetallicityEnhancementFactor = DiskMetallicityEnhancementFactor; // w.r.t to halo

/*  initializes halo radius, density, temperature profiles 
    for circumgalactic medium if needed (i.e., for CGM profiles that
    require integration to get quantities we care about. */
halo_init(this);
 
/* create fields */
NumberOfBaryonFields = 0;
DensNum = NumberOfBaryonFields;
FieldType[NumberOfBaryonFields++] = Density;
FieldType[NumberOfBaryonFields++] = TotalEnergy;
if (DualEnergyFormalism)
  FieldType[NumberOfBaryonFields++] = InternalEnergy;
vel = NumberOfBaryonFields;
FieldType[NumberOfBaryonFields++] = Velocity1;
if (GridRank > 1) 
  FieldType[NumberOfBaryonFields++] = Velocity2;
if (GridRank > 2)
  FieldType[NumberOfBaryonFields++] = Velocity3;
if (UseMHD) {
  FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
  FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
  FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
  }
if(HydroMethod == MHD_RK ){
  FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
  }
if (UseDivergenceCleaning) {
  FieldType[NumberOfBaryonFields++] = Phi_pField;
  }

/* If cosmic rays present, set up field */
CRNum = NumberOfBaryonFields;
if( CRModel )
  FieldType[NumberOfBaryonFields++] = CRDensity;

if (MultiSpecies) {
  FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
  FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
  FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
  FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
  FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
  FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
  if (MultiSpecies > 1) {
    FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
    FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
    FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
  }
  if (MultiSpecies > 2) {
    FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
    FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
    FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
  }
}

if (UseMetallicityField)
  FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity; /* fake it with metals */
if (StarMakerTypeIaSNe)
  FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;

/* Return if this doesn't concern us. */

if (ProcessorNumber != MyProcessorNumber) 
  return SUCCESS;

/* Set various units. */

float CriticalDensity = 1, BoxLength = 1;
FLOAT a, dadt, ExpansionFactor = 1;
if (ComovingCoordinates) {
  CosmologyComputeExpansionFactor(Time, &a, &dadt);
  ExpansionFactor = a/(1.0+InitialRedshift);
  CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                    &TimeUnits, &VelocityUnits, Time);
  CriticalDensity = 2.78e11*POW(HubbleConstantNow, 2); // in Msolar/Mpc^3
  BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
} else if( PointSourceGravity ){
  ENZO_FAIL("ERROR IN GALAXY SIM GRID INITIALIZE: non-cosmology units not supported for point source gravity");
} else {
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                 &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  } // end get units error if  
} // end units if/else
/*
//TEST
 printf("Testing cooling rate\n");
 int s = 6;

 float *rate = new float[s];
 float *dens = new float[s];
 float *temp = new float[s];
 float *velx = new float[s];
 float *vely = new float[s];
 float *velz = new float[s];
 float *HIdens = new float[s];
 float *HIIdens = new float[s];
 float *HeIdens = new float[s];
 float *HeIIdens = new float[s];
 float *HeIIIdens = new float[s];
 float *H2Idens = new float[s];
 float *H2IIdens = new float[s];
 float *HMdens = new float[s];
 float *edens = new float[s];
 float *metal = new float[s];

 temp[0] = 1e4/TemperatureUnits/((Gamma-1.0)*mu);
 temp[1] = 5e4/TemperatureUnits/((Gamma-1.0)*mu);
 temp[2] = 1e5/TemperatureUnits/((Gamma-1.0)*mu);
 temp[3] = 5e5/TemperatureUnits/((Gamma-1.0)*mu);
 temp[4] = 1e6/TemperatureUnits/((Gamma-1.0)*mu);
 temp[5] = 5e6/TemperatureUnits/((Gamma-1.0)*mu);

 for (i=0; i<s; ++i) {
   dens[i] = 1.007947 * 1.660538921e-24 / DensityUnits;
   velx[i] = vely[i] = velz[i] = 0.0;
   HIdens[i] = 1e-20 * dens[i];
   HIIdens[i] = 0.76 * dens[i];
   HeIdens[i] = (1.0 - 0.76) * dens[i];
   HeIIdens[i] = 1e-20 * dens[i];
   HeIIIdens[i] = 1e-20 * dens[i];
   H2Idens[i] = 1e-20 * dens[i];
   H2IIdens[i] = 1e-20 * dens[i];
   HMdens[i] = 1e-20 * dens[i];
   edens[i] = HIIdens[i] + HeIIdens[i]/4.0 + HeIIIdens[i]/2.0;
   metal[i] = 0.02041 * dens[i]; //mass fraction
 }

 this->GrackleCustomCoolRate(1, &s, rate, dens, temp,
			     velx, vely, velz, HIdens, HIIdens,
			     HeIdens, HeIIdens, HeIIIdens, edens,
			     HMdens, H2Idens, H2IIdens,
			     NULL, NULL, NULL,
			     metal);

 printf("Cooling Rates: (erg cm^3 s^-1)\n");
 for (i=0; i<s; ++i){
   rate[i] *= POW(mh,2) * POW(LengthUnits,2) / ( POW(TimeUnits,3) * DensityUnits);
   printf("%e ", rate[i]);
 }
 printf("\n");
 
 delete [] rate;
 delete [] dens;
 delete [] temp;
 delete [] velx;
 delete [] vely;
 delete [] velz;
 delete [] HIdens;
 delete [] HIIdens;
 delete [] HeIdens;
 delete [] HeIIdens;
 delete [] HeIIIdens;
 delete [] edens;
 delete [] HMdens;
 delete [] H2Idens;
 delete [] H2IIdens;
 delete [] metal;

 exit(0);
// END TEST
*/ 
/* correct background density if it's not given in code units */
if( UniformDensity < 1.0E-10 ){
  UniformDensity /= DensityUnits;
  if( debug && MyProcessorNumber == ROOT_PROCESSOR ) 
    fprintf(stdout,"Converting GalaxySimulationUniformDensity = %"GSYM" from CGS to code units\n",UniformDensity);
} // end uniform density if

/* Set up inflow */
if (GalaxySimulationInflowTime > 0.0){
  TimeActionType[0] = 2;
  TimeActionParameter[0] = GalaxySimulationInflowDensity*DensityUnits;
  TimeActionTime[0] = GalaxySimulationInflowTime*1e9/TimeUnits;
}

/* Scale gas halo rotation quantities to code units.
 * gas halo rotation variable are NOT global */
RotationScaleVelocity *= CM_PER_KM; // km/s to cm/s
RotationScaleVelocity /= LengthUnits/TimeUnits; // cm/s to code length/code time
RotationScaleRadius *= CM_PER_KPC;  // kpc to cm
RotationScaleRadius /= LengthUnits;  // cm to code length

/* compute size of fields */
size = 1;
for (dim = 0; dim < GridRank; dim++)
  size *= GridDimension[dim];

/* allocate fields */
this->AllocateGrids();

/* I'm commenting this out because the metal field should
   be set during grid initialization rather than just setting
   it as a constant color field. -- DWS */
// /* set metals to small value */
//  if (UseMetallicityField)
//    for (i = 0; i < size; i++)
//      BaryonField[MetalNum][i] = 1.0e-10;
 
/* Loop over the mesh. */
float density, disk_dens;
FLOAT halo_vmag, disk_vel[MAX_DIMENSION], Velocity[MAX_DIMENSION];
FLOAT temperature, disk_temp, init_temp, initial_metallicity;
FLOAT r_sph, x, y = 0, z = 0;
int n = 0, iter;

for (k = 0; k < GridDimension[2]; k++)
  for (j = 0; j < GridDimension[1]; j++)
    for (i = 0; i < GridDimension[0]; i++, n++) {

    if (UseMetallicityField) {
      /* Set a background metallicity value that will scale with density.
      If the cell is in the disk, this wifll be increased by a factor
      of 3.  This should really be a parameter that is read in -- DWS */ 
      initial_metallicity = GalaxySimulationGasHaloMetallicity;
    }

    /* Compute position */

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
    if (GridRank > 1)
      y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
    if (GridRank > 2)
      z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      Velocity[dim] = 0;
      disk_vel[dim] = 0;

    /* Find distance from center. */

    r_sph = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
                 POW(fabs(y-DiskPosition[1]), 2) +
                 POW(fabs(z-DiskPosition[2]), 2) );
    r_sph = max(r_sph, 0.1*CellWidth[0][0]);
    
    /*r_cyl = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
                 POW(fabs(y-DiskPosition[1]), 2) );
*/
    density = HaloGasDensity(r_sph)/DensityUnits;
    temperature = disk_temp = init_temp = HaloGasTemperature(r_sph);


    FLOAT xpos, ypos, zpos, zheight, drad; 
    float CellMass;
    FLOAT rp_hat[3];
    FLOAT yhat[3];

    /* Loop over dims if using Zeus (since vel's face-centered). */

    for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
         dim++) {

      /* Compute position. */

      xpos = x-DiskPosition[0]-(dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
      ypos = y-DiskPosition[1]-(dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
      zpos = z-DiskPosition[2]-(dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

      /* Compute z and r_perp (AngularMomentum is angular momentum 
         and must have unit length). */    

      /* magnitude of z = r.L in L direction */

      zheight = AngularMomentum[0]*xpos + 
                AngularMomentum[1]*ypos +
                AngularMomentum[2]*zpos;

      /* position in plane of disk */

      rp_hat[0] = xpos - zheight*AngularMomentum[0];
      rp_hat[1] = ypos - zheight*AngularMomentum[1];
      rp_hat[2] = zpos - zheight*AngularMomentum[2];
      drad = sqrt(rp_hat[0]*rp_hat[0] + rp_hat[1]*rp_hat[1] + rp_hat[2]*rp_hat[2]);
      drcyl = drad;

      /* Normalize the vector r_perp = unit vector pointing along plane of disk */

      rp_hat[0] = rp_hat[0]/drad;
      rp_hat[1] = rp_hat[1]/drad;
      rp_hat[2] = rp_hat[2]/drad;
      
      /* If requested, calculate velocity for CGM halo.
       * Will be replaced wtih disk velocity later if appropriate */
      if (UseHaloRotation){
          halo_vmag = RotationScaleVelocity 
                      * POW(r_sph/RotationScaleRadius, 
                            RotationPowerLawIndex);

        /* Cylindrical velocity */
        Velocity[0] = halo_vmag * (AngularMomentum[1]*rp_hat[2] -
                                   AngularMomentum[2]*rp_hat[1]);
        Velocity[1] = halo_vmag * (AngularMomentum[2]*rp_hat[0] -
                                   AngularMomentum[0]*rp_hat[2]);
        Velocity[2] = halo_vmag * (AngularMomentum[0]*rp_hat[1] -
                                   AngularMomentum[1]*rp_hat[0]);
      }
      
      if (r_sph < DiskRadius) {

        /* Find another vector perpendicular to r_perp and AngularMomentum */

        yhat[0] = AngularMomentum[1]*rp_hat[2] - AngularMomentum[2]*rp_hat[1];
        yhat[1] = AngularMomentum[2]*rp_hat[0] - AngularMomentum[0]*rp_hat[2];
        yhat[2] = AngularMomentum[0]*rp_hat[1] - AngularMomentum[1]*rp_hat[0];

        /* generate rotation matrix */
        FLOAT inv[3][3],temp;
        int i,j;

        // matrix of basis vectors in coordinate system defined by the galaxy
        inv[0][0] = rp_hat[0];
        inv[0][1] = yhat[0];
        inv[0][2] = AngularMomentum[0];
        
        inv[1][0] = rp_hat[1];
        inv[1][1] = yhat[1];
        inv[1][2] = AngularMomentum[1];
        
        inv[2][0] = rp_hat[2];
        inv[2][1] = yhat[2];
        inv[2][2] = AngularMomentum[2];

        // Matrix is orthogonal by construction so inverse = transpose

        for (i=0;i<3;i++)
          for (j=i+1;j<3;j++){
            temp = inv[i][j];
            inv[i][j] = inv[j][i];
            inv[j][i] = temp;
          }

        if( fabs(drcyl*LengthUnits/Mpc) > TruncRadius ){
          disk_dens = 0.0;
          break;
        }

        DiskDensity = (GasMass * SolarMass
                  / (8.0*pi*ScaleHeightz*Mpc*POW(ScaleHeightR*Mpc,2.0)))
                  / DensityUnits;   //Code units (rho_0) 

        if (PointSourceGravity > 0 )
          DiskVelocityMag = gasvel(drad, DiskDensity, ExpansionFactor,
                                  GalaxyMass, ScaleHeightR,
                                  ScaleHeightz, DMConcentration, Time);
        else if( DiskGravity > 0 ){
          CellMass = gauss_mass(drad*LengthUnits, zheight*LengthUnits,
                                xpos*LengthUnits, ypos*LengthUnits,
                                zpos*LengthUnits, inv, 
                                DiskDensity*DensityUnits,
                                ScaleHeightR*Mpc, ScaleHeightz*Mpc, 
                                CellWidth[0][0]*LengthUnits);

          disk_dens = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;

          DiskVelocityMag = DiskPotentialCircularVelocity(
                                                    CellWidth[0][0],
                                                    zheight*LengthUnits,
                                                    disk_dens, disk_temp);
        }
        if (PointSourceGravity*DiskGravity != FALSE ) 
          ENZO_FAIL("Cannot activate both PointSource and Disk gravity options for Isolated Galaxy");

        if (dim == 0) {
          CellMass = gauss_mass(drad*LengthUnits, zheight*LengthUnits,
                                xpos*LengthUnits, ypos*LengthUnits,
                                zpos*LengthUnits, inv, 
                                DiskDensity*DensityUnits,
                                ScaleHeightR*Mpc, ScaleHeightz*Mpc,
                                CellWidth[0][0]*LengthUnits);
          disk_dens = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;
        }

        /* If we're above the disk, then exit. */

        if (disk_dens < density)
          break;

        /* Compute velocity magnitude (divided by drad). 
           This assumes PointSourceGravityPosition and Disk center 
           are the same. */

        /* Compute velocty: L x r_perp. */
        if (dim == 0 || dim == 1)
          disk_vel[0] = DiskVelocityMag*(AngularMomentum[1]*rp_hat[2] -
                                         AngularMomentum[2]*rp_hat[1]);
        if (dim == 0 || dim == 2)
          disk_vel[1] = DiskVelocityMag*(AngularMomentum[2]*rp_hat[0] -
                                         AngularMomentum[0]*rp_hat[2]);
        if (dim == 0 || dim == 3)
          disk_vel[2] = DiskVelocityMag*(AngularMomentum[0]*rp_hat[1] -
                                         AngularMomentum[1]*rp_hat[0]);

      } // end: if (r_sph < DiskRadius)

      /* Replace CGM ("Halo") defaults with disk if dense enough; i.e.
       * replace 'density', 'temperature', 'initial_metallicity', and
       * 'Velocity' (which are currently set to CGM values) with their
       * appropriate disk values */
       
      if (disk_dens > density && fabs(drcyl*LengthUnits/Mpc) <= TruncRadius){
        density = disk_dens;
        
        /* temperature, disk_temp & init_temp start at the temp returned
         * by HaloGasTemperature(r_sph). If (r_sph < DiskRadius),
         * disk_temp may be modified by DiskPotentialCircularVelocity.
         * If it *hasn't* been modified, set it to DiskTemperature.
         * Then, replace 'temperature' with 'disk_temp' and impose
         * a ceiling */
        if (disk_temp == init_temp)
          disk_temp = DiskTemperature; 
        temperature = disk_temp;
        if( temperature > 1.0e7 )
          temperature = init_temp;
        
        /* Here we're setting the disk to be X times more enriched -- DWS */
        if( UseMetallicityField )
          initial_metallicity *= GalaxySimulationDiskMetallicityEnhancementFactor;
          
        /* Replace default/CGM velocity with disk velocity */
        Velocity[0] = disk_vel[0];
        Velocity[1] = disk_vel[1];
        Velocity[2] = disk_vel[2];
      }

    } // end: loop over dims 

    /* Set density. */

    BaryonField[0][n] = density;

    if (UseMetallicityField) {
      BaryonField[MetalNum][n] = initial_metallicity 
                                 * CoolData.SolarMetalFractionByMass 
                                 * density;
    }

    /* This should probably be scaled with density in some way to be
       a proper metallicity -- DWS (loop redundancy addressed by CEK) */
    if (StarMakerTypeIaSNe)
      BaryonField[MetalIaNum][n] = 1.0e-10;
   
    for (dim = 0; dim < GridRank; dim++)
      BaryonField[vel+dim][n] = Velocity[dim] + UniformVelocity[dim];

    /* Set energy (thermal and then total if necessary). */

    BaryonField[1][n] = temperature/TemperatureUnits / ((Gamma-1.0)*mu);

    if (DualEnergyFormalism)
      BaryonField[2][n] = BaryonField[1][n];

    if (HydroMethod != Zeus_Hydro)
      for (dim = 0; dim < GridRank; dim++)
        BaryonField[1][n] += 0.5*POW(BaryonField[vel+dim][n], 2);

    if (BaryonField[1][n] <= 0.0)
      printf("G_GSIC: negative or zero energy  n = %"ISYM"  temp = %"FSYM"   e = %"FSYM"\n",
             n, temperature, BaryonField[1][n]);

    if ( UseMHD ){
      switch ( GalaxySimulationInitialBfieldTopology ){
        case 0: //uniform
          for (dim = 0; dim < GridRank; dim++) {
            if( UseMHDCT ){
              MagneticField[dim][n] = GalaxySimulationInitialBfield[dim];
            }
            BaryonField[B1Num+dim][n] = GalaxySimulationInitialBfield[dim];
          }
          break;
          default:
          ENZO_FAIL("undefined value of GalaxySimulationInitialBfieldTopology");
      }
      BaryonField[1][n] += 0.5*(BaryonField[B1Num][n]*BaryonField[B1Num][n]
                               +BaryonField[B2Num][n]*BaryonField[B2Num][n]
                               +BaryonField[B3Num][n]*BaryonField[B3Num][n])/
                                BaryonField[0][n];
    }//UseMHD
    if( CRModel )
      BaryonField[CRNum][n] = BaryonField[DensNum][n] * GalaxySimulationCR;

      // Set multispecies fields!
      // this attempts to set them such that species conservation is maintained,
      // using the method in CosmologySimulationInitializeGrid.C
    if(MultiSpecies){
      if (MultiSpecies == 3)
	setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		   BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		   BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		   BaryonField[HMNum][n], BaryonField[H2INum][n], BaryonField[H2IINum][n],
		   BaryonField[DINum][n], BaryonField[DIINum][n], BaryonField[HDINum][n]);
      else if (MultiSpecies == 2) {
	float temp;
	setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		   BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		   BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		   BaryonField[HMNum][n], BaryonField[H2INum][n], BaryonField[H2IINum][n],
		   temp, temp, temp);
      }
      else {
	float temp;
	setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		   BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		   BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		   temp, temp, temp,
		   temp, temp, temp);
      }
    } // if(MultiSpecies)

  } // end loop over grids

  halo_clean(); // deletes halo-related arrays if needed (for circumgalactic medium)
 
  return SUCCESS;

} // end Grid::GalaxySimulationInitializeGrid

void setup_chem(float density, float temperature, int equilibrate,
		float& DEdest, float& HIdest, float& HIIdest,
		float& HeIdest, float& HeIIdest, float& HeIIIdest,
		float& HMdest, float& H2Idest, float& H2IIdest,
		float& DIdest, float& DIIdest, float& HDIdest)
{
  if (equilibrate) {
    /*  What temperature and density bins does the cell fall between? 
     *  'density' is in code units; 'temperature' is K
     *  Table should be in CGS
     */
        
    // Start by assuming values are larger than those in table;
    // set to dim_size-1 for highest available value
    bool interpolate = true;
    int dens_indx, temp_indx, iter;
    dens_indx = temp_indx = EquilibriumTable.dim_size-1;
    
    for (iter=0; iter < EquilibriumTable.dim_size; ++iter) {
      if (density < EquilibriumTable.density[iter]) {
	dens_indx = iter-1;
	break;
      }
    }
    
    for (iter=0; iter<EquilibriumTable.dim_size; ++iter) {
      if (temperature < EquilibriumTable.temperature[iter]) {
	temp_indx = iter-1;
	break;
      }
    }

    // Density or temperature lower than in table
    if (dens_indx == -1) {
      dens_indx = 0;
      interpolate = false;
    }
    if (temp_indx == -1){
      temp_indx = 0;
      interpolate = false;
    }

    // Density or temp higher than table; unchanged from inital value
    if (dens_indx == EquilibriumTable.dim_size-1 ||
	temp_indx == EquilibriumTable.dim_size-1)
      interpolate = false;

    if (interpolate) {
      HIdest = bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HIIdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIIdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIIIdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      
      DEdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.de[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.de[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      if (MultiSpecies > 1) {
	HMdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HM[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HM[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	H2Idest =  bilinear_interp(density, temperature,
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	H2IIdest =  bilinear_interp(density, temperature,
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      }
      if (MultiSpecies > 2) {
	DIdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.DI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.DI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	DIIdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.DII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.DII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	HDIdest =  bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      }
    } // end interpolate
    else { // don't interpolate; density and/or temp at edge of table
      HIdest =  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HIIdest = EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIdest = EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIIdest = EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIIIdest = EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      DEdest =  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx];

      if (MultiSpecies > 1) {
	HMdest =  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx];

	H2Idest = EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx];

	H2IIdest = EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx];
      }
      if (MultiSpecies > 2) {
	DIdest =  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx];

	DIIdest = EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx];

	HDIdest = EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx];
      }
    } // end no interpolation
  } // end if equilibrate
  else {
    HIdest = TestProblemData.HI_Fraction * density *
      TestProblemData.HydrogenFractionByMass;

    HeIdest = TestProblemData.HeI_Fraction * density *
      (1.0-TestProblemData.HydrogenFractionByMass);
       
    HeIIdest = TestProblemData.HeII_Fraction * density *
      (1.0-TestProblemData.HydrogenFractionByMass);
       
    HeIIIdest = (1.0 - TestProblemData.HydrogenFractionByMass) *
      density - HeIdest - HeIIdest;

    if(MultiSpecies > 1){
      HMdest = TestProblemData.HM_Fraction *
	TestProblemData.HydrogenFractionByMass * density;
   
      H2Idest = 2 * TestProblemData.H2I_Fraction *
	TestProblemData.HydrogenFractionByMass * density;
   
      H2IIdest = 2 * TestProblemData.H2II_Fraction 
	* TestProblemData.HydrogenFractionByMass * density;
    }

    // HII density is calculated by subtracting off the various ionized fractions
    // from the total
    HIIdest = TestProblemData.HydrogenFractionByMass * density - HIdest;
    if (MultiSpecies > 1)
      HIIdest -= (HMdest + H2IIdest + H2Idest);

    // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
    // density for convenience) is calculated by summing up all of the ionized species.
    // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
    // calculating mass density, not number density (because the BaryonField values are 4x as
    // heavy for helium for a single electron)
    DEdest = HIIdest + 0.25*HeIIdest + 0.5*HeIIIdest;
    
    if (MultiSpecies > 1)
      DEdest += 0.5*H2IIdest - HMdest;
       
    DEdest = max(DEdest, tiny_number);
       
    // Set deuterium species (assumed to be a negligible fraction of the total, so not
    // counted in the conservation)
    if(MultiSpecies > 2){
      DIdest = HIdest * TestProblemData.DeuteriumToHydrogenRatio;
      DIIdest = HIIdest *	TestProblemData.DeuteriumToHydrogenRatio;
      HDIdest = H2Idest *	0.75 * TestProblemData.DeuteriumToHydrogenRatio;
    }

  } // end not equilibrate
}

float gasvel(FLOAT radius, float DiskDensity, FLOAT ExpansionFactor, float GalaxyMass, FLOAT ScaleHeightR, FLOAT ScaleHeightz, float DMConcentration, FLOAT Time)
{

 double OMEGA=OmegaLambdaNow+OmegaMatterNow;                 //Flat Universe

 double r = radius*LengthUnits/100;    // Radius [m]

 double M_200 = GalaxyMass*SolarMass/1000.0;      // Virial Mass [kg]

 double H = sqrt(HubbleConstantNow*100*HubbleConstantNow*100*(OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor,-3)-(OMEGA-1.)*POW(ExpansionFactor,-2)));                                

 double r_200 = (1.63e-2*POW(GalaxyMass,1.0/3.0)*POW((OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor, -3)-(OMEGA-1.0)*POW(ExpansionFactor,-2)),-1.0/3.0)*ExpansionFactor*POW(H,-2.0/3.0)*POW(100,2.0/3.0))*Mpc/1.0e5;
 //virial radius [m]: M_200/M_Solar = GalaxyMass

 double M_gas, M_DM, M_Tot, Acc, V_Circ;
 double f_C = log(1.0+DMConcentration)-DMConcentration/(1.0+DMConcentration);
 double r_s = r_200/DMConcentration;  //[m]

 // Mass of gas disk and DM at given radius

     M_gas=8.0*M_PI*ScaleHeightz*Mpc/100*ScaleHeightR*Mpc/100*ScaleHeightR*Mpc/100*DiskDensity*DensityUnits*1000*PEXP(-r/(ScaleHeightR*Mpc/100))*(PEXP(r/(ScaleHeightR*Mpc/100))-r/(ScaleHeightR*Mpc/100)-1.0);

     M_DM=(M_200/f_C)*(log(1.0+r/r_s)-(r/r_s)/(1.0+r/r_s));

     if (SelfGravity==1){
  M_Tot=M_DM+M_gas;
     }
     else{
  M_Tot=M_DM;
     }

  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1;
  double MassUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
         &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  double MassUnitsDouble=1.0;

  if(ComovingCoordinates)
    MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);

  // Set the point source gravity parameters.  This is the DM mass (in g)
  //   within rs.  The core radius to rs in cm.
  //
  // BWO 10 July 2009: Both of these values are now converted to code units, because 
  // otherwise the values go over 32-bit precision.  This is used in
  // Grid::ComputeAccelerationFieldExternal, and converted back to CGS where needed.
  //

  PointSourceGravityConstant = (M_200/f_C)*(log(1.0+1.0)-1.0/(1.0+1.0))*1000.0 / MassUnitsDouble;
  PointSourceGravityCoreRadius = r_s*100.0 / LengthUnits;

  /*
  fprintf(stderr,"Grid::GalaxySimulationInitializeGrid:  %d  %e  %e\n",MyProcessorNumber,MassUnitsDouble, LengthUnits);
  fprintf(stderr,"  PointSourceGravityConstant = %e  %d\n",PointSourceGravityConstant,MyProcessorNumber);
  fprintf(stderr,"  PointSourceGravityCoreRadius = %e  %d\n",PointSourceGravityCoreRadius,MyProcessorNumber);
  */

 // Force per unit mass on disk (i.e. acceleration) [ms-2]

     Acc=((GravConst/1000.0)*M_Tot)/(r*r);

 // Magnitude of Circular Velocity of disk 

     V_Circ = sqrt(r*Acc)*100;       //cms-1

     /*      printf("r = %g  M_Tot = %g  Acc = %g  M_DM = %g  M_gas = %g  f_C = %g\n",
       r, M_Tot, Acc, M_DM, M_gas, f_C);
     printf("r_s = %g  DMConcentration = %g  r_200 = %g  r/r_s = %g\n",
       r_s, DMConcentration, r_200, r/r_s);
     printf("EF = %g  H = %g  OMEGA = %g\n", ExpansionFactor, H, OMEGA);
     printf("radius = %g  v_circ = %g\n", radius, V_Circ);  */

     return (V_Circ/VelocityUnits);  //code units
}


// Computes the total mass in a given cell by integrating the density profile using 5-point Gaussian quadrature
float gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, FLOAT cellwidth)
{
  
  FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
  FLOAT Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
  FLOAT xResult [5];
  FLOAT yResult [5];
  float Mass = 0;
  FLOAT xrot,yrot,zrot;
  int i,j,k;
  FLOAT rrot;
  
  for (i=0;i<5;i++) {

      xResult[i] = 0.0;
      for (j=0;j<5;j++) {

    yResult[j] = 0.0;
    for (k=0;k<5;k++) {

        rot_to_disk(xpos+EvaluationPoints[i]*cellwidth/2.0,ypos+EvaluationPoints[j]*cellwidth/2.0,zpos+EvaluationPoints[k]*cellwidth/2.0,xrot,yrot,zrot,inv);
        rrot = sqrt(POW(xrot,2)+POW(yrot,2));

        if( PointSourceGravity > 0 )
    yResult[j] += cellwidth/2.0*Weights[k]*PEXP(-rrot/ScaleHeightR)/POW(cosh(zrot/(2.0*ScaleHeightz)),2);
        else if( DiskGravity > 0 ){
    if( rrot/Mpc < SmoothRadius )
      yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz);
    else if( rrot/Mpc < TruncRadius )
      yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz)*0.5*(1.0+cos(pi*(rrot-SmoothRadius*Mpc)/(SmoothLength*Mpc)));
        } // end disk gravity if

    }
    xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
      }
      Mass += cellwidth/2.0*Weights[i]*xResult[i];
  }  
  Mass *= DiskDensity;
  return Mass;
}

//Finds coordinates in rotated coordinate system
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3])
{
  xrot = xpos*inv[0][0] + ypos*inv[0][1] + zpos*inv[0][2];
  yrot = xpos*inv[1][0] + ypos*inv[1][1] + zpos*inv[1][2];
  zrot = xpos*inv[2][0] + ypos*inv[2][1] + zpos*inv[2][2];
}


double DiskPotentialDarkMatterMass(FLOAT R){
/*
 *  computes dark matter mass enclosed within spherical radius R
 *  for potential in Mori & Burkert 2000, consistent with eq
 *
 *    rho = rho0 * r0**3 / ( (r + r0)*(r**2 + r0**2 ) )
 *      
 *  Parameters:
 *  -----------
 *    R - Spherical radius (code units)
 *
 *  Returns: Mass, in grams
 */
  FLOAT R0 = DiskGravityDarkMatterR*Mpc,x=R/R0*LengthUnits;
  double M0 = pi*DiskGravityDarkMatterDensity*R0*R0*R0;

  return M0*(-2.0*atan(x)+2.0*log(1+x)+log(1.0+x*x));
} // end DiskPotentialDarkMatterMass


/* -------------------- BEGINNING OF Routines used for initializing the circumgalactic medium -------------------- */
/* 
   Computes halo gas density values assuming a variety of user-specifiable models
   for the CGM, toggled by the variable GalaxySimulationGasHalo.  Depending on the
   specific model chosen, different global parameters are needed (as set near the beginning
   of Grid::GalaxySimluationInitializeGrid).  Halo types are:

   GalaxySimulationGasHalo = 0  -- "zero CGM" - sets to a very low density/temperature
   GalaxySimulationGasHalo = 1  -- assuming hydrostatic equilibrium of CGM given an NFW dark matter halo 
                                   and a temperature as a function of radius set by the virial theorem.
   GalaxySimulationGasHalo = 2  -- assumes density, temperature set according to T = Tvir and entropy
                                   as a power-law function of radius.
   GalaxySimulationGasHalo = 3  -- as #2, but the entropy distribution has a floor value, so S = S_f + S_0 (r/r_0)^alpha
   GalaxySimulationGasHalo = 4  -- assumes a hydrostatic equilibrium of CGM given an NFW dark matter halo
                                   and an entropy that is a power-law function of radius.
   GalaxySimulationGasHalo = 5  -- as #4, but the entropy distribution has a floor value, so S = S_f + S_0 (r/r_0)^alpha
   GalaxySimulationGasHalo = 5  -- as #4, but the entropy distribution follows that for a precipitation-regulated NFW halo
                                   in Voit 2019 (ApJ)
        

   Inputs:  R - spherical radius, code units

   Returns:  density, grams/cm^3

   Note: using global variables w/following units:

   GalaxySimulationGasHalo: integer, >= 0
   GalaxySimulationGasHaloScaleRadius, units of Mpc
   GalaxySimulationGasHaloDensity, units of grams/cm^3
   GalaxySimulationGasHaloTemperature, units of Kelvin
   GalaxySimulationGasHaloAlpha, power-law index; unitless
   GalaxySimulationGasHaloCoreEntropy, units of keV cm^2
   GalaxySimulationGasHaloMetallicity, units of Zsun
*/
float HaloGasDensity(FLOAT R){

  if(GalaxySimulationGasHalo < 1){
    /* "zero CGM" - sets a very low density */
   
    return densicm;

  } else if(GalaxySimulationGasHalo == 1){
    /* gets density assuming hydrostatic equilibrium using a temperature 
       as a function of radius given by virial theorem */
    
    double T0,haloDensity;
    T0 = HaloGasTemperature(GalaxySimulationGasHaloScaleRadius*Mpc/LengthUnits);
    haloDensity = GalaxySimulationGasHaloDensity*(T0/HaloGasTemperature(R));
    haloDensity /= POW((R*LengthUnits/GalaxySimulationGasHaloScaleRadius/Mpc),3);
    return min(haloDensity,GalaxySimulationGasHaloDensity);
    
  } else if(GalaxySimulationGasHalo == 2){
    /* assumes entropy is a power-law function of radius and T = Tvir, so
       n(r) = n_0 * (r/r_0)**(-alpha/(gamma-1))
       where n_0 is (user-supplied) number density at (user-supplied) radius r_0,
       alpha is (user-supplied) power-law exponent,
       gamma is adiabatic index.
    */
    double scale_radius_cgs, this_radius_cgs, power_law_exponent;
    
    scale_radius_cgs = GalaxySimulationGasHaloScaleRadius*Mpc;
    this_radius_cgs = R*LengthUnits;
    power_law_exponent = -1.0*GalaxySimulationGasHaloAlpha/(Gamma-1.0);
    
    return GalaxySimulationGasHaloDensity*POW(this_radius_cgs/scale_radius_cgs, power_law_exponent);
    
  } else if(GalaxySimulationGasHalo == 3){
    /* assumes entropy is a  power-law function of radius and T = Tvir that has a core (i.e., minimum 
       entropy value), so:

       n(r) = (Tvir / (Score + S_0*(r/r_0)^alpha))^(1/(gamma-1))

       where n_0 is (user-supplied) number density at (user-supplied) radius r_0,
       Tvir is user-supplied temperature,
       Score is a user-supplied entropy,
       alpha is (user-supplied) power-law exponent,
       gamma is adiabatic index.
    */

    double scale_radius_cgs, this_radius_cgs, power_law_exponent, S_0, T_kev, n_0, this_number_density;

    // get radii in common set of units
    scale_radius_cgs = GalaxySimulationGasHaloScaleRadius*Mpc;
    this_radius_cgs = R*LengthUnits;

    // now get number density using expression above

    T_kev = GalaxySimulationGasHaloTemperature*8.6174e-08;  // halo temperature in keV
    n_0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert n_0 to electron number density 
    S_0 = T_kev / POW(n_0,Gamma-1.0);   // S_0 in units of kev cm^2

    // get number density at this radius giving the requested info
    this_number_density = T_kev/(GalaxySimulationGasHaloCoreEntropy + S_0*POW(this_radius_cgs/scale_radius_cgs, power_law_exponent));
    this_number_density = POW(this_number_density, 1.0/(Gamma-1.0));

    return this_number_density*mu*mh;  // return physical density
    
  } else if(GalaxySimulationGasHalo == 4 || GalaxySimulationGasHalo == 5 || GalaxySimulationGasHalo == 6){
    /* assumes entropy is a power-law function of radius OR a cored power-law function
       of radius and gas is in hydrostatic equilibrium w/the NFW halo.  */

    double this_radius_cgs;
    int index;
    this_radius_cgs = R*LengthUnits;  // radius in CGS
    index = int(this_radius_cgs/CGM_data.dr+1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    return CGM_data.n_rad[index]*mu*mh;  // return physical density

  } else if(GalaxySimulationGasHalo == 7){
    /* Eqn 24 in the Appendix of Voit 2019; a fit to the theoretical density profile of a precipitation-regulated NFW halo.
       Equation gives the electron number density, but pull the same trick as methods 2 & 3 and assume n_e = n */
    double this_radius_kpc, this_number_density;
    this_radius_kpc = R*LengthUnits/CM_PER_KPC;
    
    this_number_density = POW( POW(this_radius_kpc,GalaxySimulationGasHaloZeta) / GalaxySimulationGasHaloDensity, 2);
    this_number_density += POW( POW(this_radius_kpc/100,GalaxySimulationGasHaloZeta2) / GalaxySimulationGasHaloDensity2, 2);
    this_number_density = POW(this_number_density, -0.5);

    return this_number_density*mu*mh;  // return physical density
    
  } else {
    ENZO_FAIL("Grid::GalaxySimulationInitializeGrid - invalid choice of GalaxySimulationGasHalo in HaloGasDensity().");
  }
  
} // end HaloGasDensity


/* 
   Computes halo gas temperature values assuming a variety of user-specifiable models
   for the CGM, toggled by the variable GalaxySimulationGasHalo.  The properties of each of the
   models are described immediately above this in the comments for the function HaloGasDensity(). 

   Inputs:  R - spherical radius, code units

   Returns:  Temperature, Kelvin
*/

float HaloGasTemperature(FLOAT R){

  if(GalaxySimulationGasHalo < 1){
    /* "zero CGM" - sets a very low temperature */
   
    return Ticm;

  } else if(GalaxySimulationGasHalo == 1){
    /* gets temperature as a function of radius given by virial theorem */

    return GravConst*DiskPotentialDarkMatterMass(R)*mu*mh/(3.0*kboltz*R*LengthUnits);
    
  } else if(GalaxySimulationGasHalo == 2){
    /* assumes entropy is a power-law function of radius and T = Tvir */

    return GalaxySimulationGasHaloTemperature;
    
  } else if(GalaxySimulationGasHalo == 3){

    /* assumes entropy is a cored power-law function of radius and T = Tvir */

    return GalaxySimulationGasHaloTemperature;

  } else if(GalaxySimulationGasHalo == 4 || GalaxySimulationGasHalo == 5 || GalaxySimulationGasHalo == 6){
    /* assumes entropy is a power-law function of radius and gas is in hydrostatic equilibrium */

    double this_radius_cgs;
    int index;
    this_radius_cgs = R*LengthUnits;  // radius in CGS
    index = int(this_radius_cgs/CGM_data.dr+1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    return CGM_data.T_rad[index];  // return temperature in Kelvin

  } else if(GalaxySimulationGasHalo == 7){
    /* Theoretical temperature profile of a precipitation-regulated NFW halo, using fits to n(r) and S(r) */
    double this_radius_kpc, this_number_density, this_entropy;
    this_radius_kpc = R*LengthUnits/CM_PER_KPC;
    
    this_number_density = POW( POW(this_radius_kpc,GalaxySimulationGasHaloZeta) / GalaxySimulationGasHaloDensity, 2);
    this_number_density += POW( POW(this_radius_kpc/100,GalaxySimulationGasHaloZeta2) / GalaxySimulationGasHaloDensity2, 2);
    this_number_density = POW(this_number_density, -0.5);

    this_entropy = GalaxySimulationGasHaloCoreEntropy * POW(this_radius_kpc, GalaxySimulationGasHaloAlpha);

    return this_entropy * POW(this_number_density, Gamma-1.0) / kboltzKeV; // units of K
    
  } else {
    ENZO_FAIL("Grid::GalaxySimulationInitializeGrid - invalid choice of GalaxySimulationGasHalo in HaloGasTemperature().");
  }
  
}

/* Initializes arrays of number density, temperature, and radius for 
   choices of circumgalactic medium that require numerical integration based
   on user-defined parameters.  These quantities are stored in a global struct
   for convenience (global within this file, at least). */
void halo_init(grid* Grid){

  if(GalaxySimulationGasHalo < 4 || GalaxySimulationGasHalo > 6) return;

  double k1, k2, k3, k4, this_n, this_radius, temperature;
  double M, Rvir, rho_crit = 1.8788e-29*0.49, Tvir, n0, r0, dr;

  int index;
  
  CGM_data.nbins = 8192;
  
  CGM_data.n_rad = new double[CGM_data.nbins];
  CGM_data.T_rad = new double[CGM_data.nbins];
  CGM_data.rad   = new double[CGM_data.nbins];

  for(int i=0;i<CGM_data.nbins;i++) CGM_data.n_rad[i]=CGM_data.T_rad[i]=CGM_data.rad[i]=-1.0;
  
  M = GalaxySimulationGasHaloGalaxyMass * SolarMass;  // halo total mass in CGS
  
  Rvir = pow(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS

  CGM_data.R_outer = Rvir;  // integrate out to the virial radius of halo

  CGM_data.dr = CGM_data.R_outer / double(CGM_data.nbins);  // stepsize for RK4 integration and radial bins
  
  // set some quantities based on user inputs; this defines our integration
  Tvir = GalaxySimulationGasHaloTemperature;
  n0 = GalaxySimulationGasHaloDensity / (mu*mh);
  r0 = GalaxySimulationGasHaloScaleRadius*Mpc;

  // used for our numerical integration
  dr = CGM_data.dr;
  this_n = n0;
  this_radius = r0;

  // set the bin that we start at (otherwise it doesn't get set!)
  index = int(this_radius/dr+1.0e-3);
  CGM_data.n_rad[index] = this_n;
  CGM_data.T_rad[index] = Tvir;
  CGM_data.rad[index] = this_radius;

  /* starting at the point where the user has defined the radius, density, and 
     temperature, use RK4 to integrate the number density outward to R_outer using the expression 
     for dn_dr in another function.  Calculate the temperature using the entropy at this radius. */
  while(this_radius <= CGM_data.R_outer){
    
    // calculate RK4 coefficients.
    k1 = halo_dn_dr(this_radius,         this_n, Grid);
    k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1, Grid);
    k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2, Grid);
    k4 = halo_dn_dr(this_radius + dr,     this_n + dr*k3, Grid);

    // update density and radius
    this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4);
    this_radius += dr;  // new radius

    // calculate temperature at this radius using known entropy 
    temperature = halo_S_of_r(this_radius,this_n,Grid) * POW(this_n,Gamma-1.0);

    // store everything in the struct
    index = int(this_radius/dr+1.0e-3);    
    CGM_data.n_rad[index] = this_n;
    CGM_data.T_rad[index] = temperature;
    CGM_data.rad[index] = this_radius;
  }


  /* now we do the same thing as above, but integrating intward to zero radius. */
  this_n = n0;
  this_radius = r0;
  dr *= -1.0;

  while(this_radius > 0.0){
    
    // calculate RK4 coefficients.
    k1 = halo_dn_dr(this_radius,         this_n, Grid);
    k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1, Grid);
    k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2, Grid);
    k4 = halo_dn_dr(this_radius + dr,     this_n + dr*k3, Grid);

    // update density and radius
    this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4);
    this_radius += dr;  // new radius

    // calculate temperature at this radius using known entropy 
    temperature = halo_S_of_r(this_radius, this_n, Grid) * POW(this_n,Gamma-1.0);

    // store everything in the struct
    index = int(this_radius/(-1.0*dr)+1.0e-3);

    if(index >= 0){
      CGM_data.n_rad[index] = this_n;
      CGM_data.T_rad[index] = temperature;
      CGM_data.rad[index] = this_radius;
    }
  }

  // this integration acts a little squirrelly around r=0 because the mass values are garbage.  Cheap fix.
  CGM_data.rad[0]=CGM_data.rad[1];
  CGM_data.n_rad[0]=CGM_data.n_rad[1];
  CGM_data.T_rad[0]=CGM_data.T_rad[1];
  
  return;
}

/* If we declared these arrays, clean them up at the end of the problem initialization. 
   Software campsite principle. */
void halo_clean(void){

  if(GalaxySimulationGasHalo < 4 || GalaxySimulationGasHalo > 6) return;

  delete [] CGM_data.n_rad;
  delete [] CGM_data.T_rad;
  delete [] CGM_data.rad;

  return;
}

/* Halo entropy as a function of radius for the user-specified CGM types that require numerical
   integration. 

   Input is radius in CGS units.  output is entropy in CGS units (Kelvin cm^2) 
*/
double halo_S_of_r(double r, double n, grid* Grid){

  double Tvir, n0, r0, Smin, S0;

  // calculate a bunch of things based on user inputs
  Tvir = GalaxySimulationGasHaloTemperature;  // in Kelvin
  n0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert from density to electron number density (cm^-3)
  r0 = GalaxySimulationGasHaloScaleRadius*Mpc;  // scale radius in CGS
  Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8;  // given in keV cm^2, converted to Kelvin cm^2
  S0 = Tvir / POW(n0,Gamma-1);  // entropy at scale radius, in units of Kelvin cm^2

  if(GalaxySimulationGasHalo == 4){

    return S0*POW(r/r0,GalaxySimulationGasHaloAlpha);  // has units of Kelvin cm^2

  } else if (GalaxySimulationGasHalo == 5){

    return Smin + S0*POW(r/r0,GalaxySimulationGasHaloAlpha);  // has units of Kelvin cm^2

  } else if (GalaxySimulationGasHalo == 6){

    double vcirc2 = GravConst * halo_galmass_at_r(r) / r;
    double Tgrav = mu*mh * vcirc2 / (2*kboltz); // gravitational "temperature"
    double therm = 2*Tgrav / kboltz; // ?
    double dens = n*mu*mh; // current density
    double Lambda, vx=0, vy=0, vz=0;
    double hi, hii, hei, heii, heiii, de, hm, h2i, h2ii, di, dii, hdi, metal; // species
    int dim=1;

    setup_chem(dens, 2*Tgrav, EquilibrateChem, de, hi, hii, hei, heii, heiii, hm, h2i, h2ii, di, dii, hdi);
    metal = GalaxySimulationGasHaloMetallicity * CoolData.SolarMetalFractionByMass * dens;
    
    Grid->GrackleCustomCoolRate(1, &dim, &Lambda, &dens, &therm, &vx, &vy, &vz, &hi, &hii,
			  &hei, &heii, &heiii, &hm, &h2i, &h2ii, &di, &dii, &hdi, &metal);
    
  } else {
    ENZO_FAIL("halo_S_of_r: GalaxySimulationGasHalo set incorrectly.");
  }

}

/* dEntropy/dr as a function of radius for the user-specified CGM types that require numerical
   integration. 

   Input is radius in CGS units; output is entropy gradient in CGS units (Kelvin*cm) */
double halo_dSdr(double r){

  double Tvir, alpha, n0, r0, Smin, S0;

  // calculate a bunch of things based on user inputs
  Tvir = GalaxySimulationGasHaloTemperature;  // in Kelvin
  n0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert from density to electron number density (cm^-3)
  r0 = GalaxySimulationGasHaloScaleRadius*Mpc;  // scale radius in CGS
  Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8;  // given in keV cm^2, converted to Kelvin cm^2
  S0 = Tvir / POW(n0,Gamma-1);  // entropy at scale radius, in units of Kelvin cm^2

  if(GalaxySimulationGasHalo == 4 || GalaxySimulationGasHalo == 5){

    // has units of Kelvin*cm, same deriv for both halo types (since constant drops out)
    return S0*GalaxySimulationGasHaloAlpha*
      POW(r/r0,GalaxySimulationGasHaloAlpha-1.0)/r0;

  } else if (GalaxySimulationGasHalo == 6){
    
    return 0;//SOMETHING;

  } else {
    ENZO_FAIL("halo_dSdr: GalaxySimulationGasHalo set incorrectly.");
  }
}

/* dn/dr as a function of radius and halo electron number density.  This quantity is calculated
   by assuming that gravity and pressure are in hydrostatic equilibrium in a halo with a specified 
   entropy profile S(r).

   Input is radius in CGM units and electron number density in units 
   of particles per cm^-3.  Output is dn/dr in CGS units, so particles per cm^4. */
double halo_dn_dr(double r, double n, grid* Grid){
  
  return -1.0*( n*1.22*mh*halo_g_of_r(r) + kboltz*POW(n,Gamma)*halo_dSdr(r) ) /
    ( Gamma * kboltz * halo_S_of_r(r,n,Grid) * POW(n, Gamma-1));
  
}

/* halo gravitational acceleration as a function of radius.

   Input is the radius in CGS units and returns the MAGNITUDE of the 
   acceleration in CGS units.  */
double halo_g_of_r(double r){
  return GravConst*halo_galmass_at_r(r)/(r*r); 
}

/* halo galaxy mass at a given radius, using user-defined global parameters for galaxy
   quantities and assuming that all halo mass is in an NFW halo.  This is not totally
   correct near the center of the halo, but since we're using it for the CGM initialization 
   and are dealing with radii that aren't particularly near the center of the halo, this 
   approximation is probably fine. 

   Input is the radius in CGS units; output is the enclosed mass at that radius in CGS units.
*/
double halo_galmass_at_r(double r){

  double M, C, Rvir, rho_0, Rs, M_within_r;
  double rho_crit = 1.8788e-29*0.49;
  
  M = GalaxySimulationGasHaloGalaxyMass * SolarMass;  // halo total mass in CGS
  C = GalaxySimulationGasHaloDMConcentration;  // concentration parameter for NFW halo
  
  Rvir = POW(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  Rs = Rvir/C;  // scale radius of NFW halo in CGS
  rho_0 = 200.0*POW(C,3)/3.0/(log(1.0+C) - C/(1.0+C))*rho_crit;  // rho_0 for NFW halo in CGS

  // mass w/in radius R
  M_within_r = 4.0*3.14159*rho_0*POW(Rs,3.0)*(log((Rs+r)/Rs) - r/(Rs+r));

  return M_within_r;
  
}

/* -------------------- END OF Routines used for initializing the circumgalactic medium -------------------- */


float DiskPotentialGasDensity(FLOAT r,FLOAT z){
/*
 *  computes gas density within galaxy disk, according to eq
 *
 *    (Mgas/8*pi*a^2*b)*sech(r/a)*sech*(z/b)
 *
 *  Smoothed by a cosine fcn beyond SmoothRadius
 *
 *  Parameteres:
 *  ------------
 *    r - cylindrical radius (code units)
 *    z - cylindrical height (code units)
 *
 *  Returns: density (in grams/cm^3)
 *
 */
  double density = MgasScale*SolarMass/(8.0*pi*POW(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc);
  density /= (cosh(r*LengthUnits/gScaleHeightR/Mpc)*cosh(z*LengthUnits/gScaleHeightz/Mpc));

  if(fabs(r*LengthUnits/Mpc) > SmoothRadius && fabs(r*LengthUnits/Mpc) <= TruncRadius)
    density *= 0.5*(1.0+cos(pi*(r*LengthUnits-SmoothRadius*Mpc)/(SmoothLength*Mpc)));
  return density;
} // end DiskPotentialGasDensity




double findZicm(FLOAT r){
  /*  
   *  Finds the height above the disk plane where the disk gas density
   *  matches the halo's gas density (using bisection)
   *
   *  Parameters:
   *  -----------
   *    r - cylindrical radius (code units)
   *
   *  Returns: zicm, edge of disk, (code units)
   */

  static const double X_TOL = 1e-7*Mpc/LengthUnits; // sub pc resolution
  static const int MAX_ITERS = 50; int iters=0;

  double z_lo = 0.0,z_hi = 0.01*Mpc/LengthUnits,z_new,f_lo,f_hi,f_new;
  f_hi = DiskPotentialGasDensity(r,z_hi) - HaloGasDensity(sqrt(r*r+z_hi*z_hi)); // -ve
  f_lo = DiskPotentialGasDensity(r,z_lo) - HaloGasDensity(sqrt(r*r+z_lo*z_lo)); // +ve

  if(f_lo < 0.0) return 0.0; // beyond the disk
  if(f_hi > 0.0) ENZO_FAIL("ERROR IN GALAXY INITIALIZE: HALO IS UNDER-PRESSURIZED");

  while(iters++ < MAX_ITERS ){

    z_new = (z_hi+z_lo)/2.0;
    f_new = DiskPotentialGasDensity(r,z_new)
            - HaloGasDensity(sqrt(r*r+z_new*z_new));

    if( fabs(f_new) == 0.0 ) return z_new;
    if( f_new*f_lo > 0.0 ){
      z_lo = z_new; f_lo = f_new;
    }
    else{
      z_hi = z_new; f_hi = f_new;
    }
    if( fabs(z_hi - z_lo) <= X_TOL ) return z_new;
  }

  ENZO_FAIL("ERROR IN GALAXY INITIALIZE: findZicm FAILED TO CONVERGE");
  return -1.0;
}


/* 
 *  DISK POTENTIAL CIRCULAR VELOCITY
 *
 *      Returns disk circular velocity (in code units) given height z
 *      and drcyl (radius in plane) in disk.  This includes the effect of
 *      thermal pressure and so is not just the sqrt(G M/drcyl).
 *      Note that for historical reasons drcyl is an external. *
 */
float DiskPotentialCircularVelocity(FLOAT cellwidth, FLOAT z, FLOAT density, 
            FLOAT &temperature)
{

  extern double drcyl;
  double PbulgeComp1(double zint);       // (density times Stellar bulge force)
  double PbulgeComp2(double zint);       // same but for r2 (3D distance)
  double PstellarComp1(double zint);     // (density times stellar disk force)
  double PstellarComp2(double zint);     // same but for r2 (3D distance plane)
  double PDMComp1(double zint);          // (density times dark matter halo force)
  double PDMComp2(double zint);          // same but for r2 (3D distance plane)

  double Pressure,Pressure2,zicm,zicm2,zicmf=0.0,zsmall=0.0,
    zicm2f=0.0,zint,FdPdR,FtotR,denuse,rsph,vrot,bulgeComp,rsph_icm;

  r2 = (drcyl+0.01*cellwidth)*LengthUnits;  // in plane radius
  rsph = sqrt(POW(drcyl*LengthUnits,2)+POW(z,2)); // 3D radius

  /*  Determine zicm: the height above the disk where rho -> rho_ICM,
   *  use this to find P_icm and dP_icm  */

  if (fabs(drcyl*LengthUnits/Mpc) <= SmoothRadius) {

    zicm  = findZicm(drcyl)*LengthUnits;
    zicm2 = findZicm(r2/LengthUnits)*LengthUnits;

    if( fabs(z) < fabs(zicm) ){

      /* Integrate the density times force to get pressure.  Do this
   at two different locations to get a numerical gradient. */

      bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ? 
       0.0 : qromb(PbulgeComp1, fabs(zicm), fabs(z)));
      Pressure  = bulgeComp + qromb(PstellarComp1, fabs(zicm), fabs(z));
      Pressure += qromb(PDMComp1, fabs(zicm), fabs(z));

      bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ? 
       0.0 : qromb(PbulgeComp2, fabs(zicm2), fabs(z)));
      Pressure2  = bulgeComp + qromb(PstellarComp2, fabs(zicm2), fabs(z));
      Pressure2 += qromb(PDMComp2, fabs(zicm2), fabs(z));

    }  // end |z| < |zicm| if

  }  else {

    if (fabs(drcyl*LengthUnits/Mpc) <= TruncRadius ) {

      zicm  = findZicm(drcyl)*LengthUnits;
      zicm2 = findZicm(r2/LengthUnits)*LengthUnits;

      /* This checks to see if the returned density is smaller than the halo
   density and issues warning. */

#ifdef UNUSED
      if ( HaloGasDensity(sqrt(drcyl*drcyl+z*z)) >= DiskPotentialGasDensity(drcyl,z)
     && fabs(z) < zicm) {
  printf("warning: small density zicm = %g, z = %g\n", zicm/Mpc, z/Mpc);
      } // end small density if
#endif /* UNUSED */

      if (fabs(z) < fabs(zicm)) {
        
  bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ?
         0.0 : qromb(PbulgeComp1, fabs(zicm), fabs(z)));
  Pressure  = (bulgeComp + qromb(PDMComp1, fabs(zicm), fabs(z)) + qromb(PstellarComp1, fabs(zicm), fabs(z)))
    *(0.5*(1.0+cos(pi*(drcyl*LengthUnits-SmoothRadius*Mpc)/
       (SmoothLength*Mpc))));

  bulgeComp = (DiskGravityStellarBulgeMass == 0.0 ?
         0.0 : qromb(PbulgeComp2, fabs(zicm2), fabs(z)));
  Pressure2 = (bulgeComp + qromb(PDMComp2, fabs(zicm2), fabs(z)) + qromb(PstellarComp2, fabs(zicm2), fabs(z)))
    *(0.5*(1.0+cos(pi*(r2-SmoothRadius*Mpc)/(SmoothLength*Mpc))));

      } // end |z| < |zicm| if

    } // end r_cyle < TruncRadius if

  } // end r_cyl < SmoothRadius if/else

  denuse = density*DensityUnits; 

  if (Pressure < 0.0 && fabs(drcyl)*LengthUnits/Mpc <= TruncRadius && fabs(z) <= fabs(zicm)) {
    fprintf(stderr,"neg pressure:  P = %"FSYM", z = %"FSYM", r = %"FSYM"\n", Pressure, z/Mpc, drcyl*LengthUnits/Mpc);
  }
  if (fabs(drcyl)*LengthUnits/Mpc >= TruncRadius || fabs(zicm) <= fabs(z)){
    Pressure = 0.0;
    Pressure2 = 0.0;
    denuse = HaloGasDensity(rsph);
  }
  if (Pressure2 <= 0.0 && Pressure <= 0.0){
    Pressure = 0.0;
    Pressure2 = 0.0;
    denuse = HaloGasDensity(rsph);
  }
  if (Pressure <= 0.0) {
    Pressure = 0.0;
    Pressure2 = 0.0;
    denuse = HaloGasDensity(rsph);
  }
  if (denuse < HaloGasDensity(rsph)) {
    fprintf(stderr,"denuse small:  %"FSYM"\n", denuse);
  }
  rsph_icm = sqrt(drcyl*drcyl+POW(zicm/LengthUnits,2));
  Picm = HaloGasDensity(rsph_icm)*kboltz*HaloGasTemperature(rsph_icm)/(mu*mh);
  temperature=mu*mh*(Picm+Pressure)/(kboltz*denuse);

  /* Calculate pressure gradient */

  FdPdR = (Pressure2 - Pressure)/(r2-drcyl*LengthUnits)/density; 

  /* Calculate Gravity = Fg_DM + Fg_StellarDisk + Fg_StellaDiskGravityStellarBulgeR */
  
  FtotR  = (-pi)*GravConst*DiskGravityDarkMatterDensity*
          POW(DiskGravityDarkMatterR*Mpc,3)/POW(rsph,3)*drcyl*LengthUnits
    *(-2.0*atan(rsph/DiskGravityDarkMatterR/Mpc) + 
      2.0*log(1.0+rsph/DiskGravityDarkMatterR/Mpc) +
      log(1.0+POW(rsph/DiskGravityDarkMatterR/Mpc,2)));
  FtotR += -GravConst*DiskGravityStellarDiskMass*SolarMass*drcyl*LengthUnits
    /sqrt(POW(POW(drcyl*LengthUnits,2) + 
        POW(DiskGravityStellarDiskScaleHeightR*Mpc +
      sqrt(POW(z,2) + 
           POW(DiskGravityStellarDiskScaleHeightz*Mpc,2)),
      2),
        3));
  FtotR += -GravConst*DiskGravityStellarBulgeMass*SolarMass
           /POW(sqrt(POW(z,2) + POW(drcyl*LengthUnits,2)) + 
          DiskGravityStellarBulgeR*Mpc,2)*drcyl*LengthUnits/
           sqrt(POW(z,2) +POW(drcyl*LengthUnits,2));

  /* Some error checking. */

  if (temperature < 0.0) 
    fprintf(stderr,"G_GSIG: temp = %"FSYM", P = %"FSYM", z = %"FSYM", zicm = %"FSYM", zicmf=%"FSYM", zsmall=%"FSYM", drcyl = %"FSYM"\n", 
      temperature, Pressure, z/Mpc, zicm/Mpc, zicmf, zsmall, drcyl*LengthUnits/Mpc);
  if ((FtotR - FdPdR) > 0.0) { 
    fprintf(stderr,"G_GSIG: FtotR = %"FSYM", FdPdR = %"FSYM", P = %"FSYM",P2 = %"FSYM", Picm = %"FSYM", dr = %"FSYM", drcyl = %"FSYM", z = %"FSYM"\n", 
      FtotR, FdPdR, Pressure, Pressure2, Picm, r2-drcyl*LengthUnits, drcyl*LengthUnits/Mpc, z/Mpc);
    FdPdR = 0.0;
  } // end FtotR - FdPdr > 0.0 if

  /* Find circular velocity by balancing FG and dPdR against centrifugal force */

  vrot=sqrt(-drcyl*LengthUnits*(FtotR-FdPdR));

  if (denuse == densicm) vrot = 0.0;

  return (vrot/VelocityUnits); //code units

} // end DiskPotentialCircularVelocity


// *************************************************************
// The functions integrated by qromb (parameter must be external)

// Given the in place radius and height, this returns density times
//    the stellar bulge force

double PbulgeComp_general(double rvalue, double zint)
{
  return (-MgasScale*SolarMass/
    (2*pi*POW(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/
    cosh(rvalue/gScaleHeightR/Mpc) / cosh(fabs(zint)/gScaleHeightz/Mpc)*
    GravConst*DiskGravityStellarBulgeMass*SolarMass/
    POW((sqrt(POW(zint,2) + POW(rvalue,2)) + 
         DiskGravityStellarBulgeR*Mpc),2)*fabs(zint)/
    sqrt(POW(zint,2)+POW(rvalue,2)));
}

// Stellar Bulge functions

double PbulgeComp1(double zint)
{
  extern double drcyl;
  return PbulgeComp_general(drcyl*LengthUnits, zint);
}

double PbulgeComp2(double zint)
{
  extern double r2;
  return PbulgeComp_general(r2, zint);
}

double PstellarComp_general(double rvalue, double zint)
{
  return (-MgasScale*SolarMass/
    (2*pi*POW(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/
    cosh(rvalue/gScaleHeightR/Mpc) / cosh(fabs(zint)/gScaleHeightz/Mpc)*
    GravConst*DiskGravityStellarDiskMass*SolarMass*
    (DiskGravityStellarDiskScaleHeightR*Mpc + 
     sqrt(POW(zint,2) + POW(DiskGravityStellarDiskScaleHeightz*Mpc,2)))*
    fabs(zint)/
    sqrt(POW(POW(rvalue,2) + 
       POW((DiskGravityStellarDiskScaleHeightR*Mpc + 
      sqrt(POW(zint,2) + 
           POW(DiskGravityStellarDiskScaleHeightz*Mpc,2)))
           ,2)
       ,3))/
    sqrt(POW(zint,2)+POW(DiskGravityStellarDiskScaleHeightz*Mpc,2)));
}

// Stellar Disk functions

double PstellarComp1(double zint)
{
  extern double drcyl;
  return PstellarComp_general(drcyl*LengthUnits, zint);
}

double PstellarComp2(double zint)
{
  extern double r2;
  return PstellarComp_general(r2, zint);
}

double PDMComp_general(double rvalue, double zint){
/* --------------------------------------------------------
 * PDMComp_general
 * --------------------------------------------------------
 * General function for computing the DM contribution to
 * local (vertical) pressure on the gas in the galaxy's disk.
 * This returns the gas density at a given position
 * times the force on the gas due to the dark matter
 * --------------------------------------------------------- */

  float gas_density;
  float F;             // dark matter force
  float rsph; // 3D, spherical radius


  /* compute gas density */
  gas_density  = MgasScale*SolarMass / (8.0 * pi * POW(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc);
  gas_density /= (cosh(rvalue/gScaleHeightR/Mpc)*cosh(fabs(zint)/gScaleHeightz/Mpc));

  rsph = sqrt(rvalue*rvalue + zint*zint);

  /* fabs(zint) is because this is the force in the direction downward */
  F    = (-pi)*GravConst*DiskGravityDarkMatterDensity*
             POW(DiskGravityDarkMatterR*Mpc,3)/POW(rsph,3) * fabs(zint)
            *(-2.0*atan(rsph/DiskGravityDarkMatterR/Mpc) + 
            2.0*log(1.0+rsph/DiskGravityDarkMatterR/Mpc) +
            log(1.0+POW(rsph/DiskGravityDarkMatterR/Mpc,2)));

  return gas_density * F;
}

/* DM pressure integration */
double PDMComp1(double zint){
  extern double drcyl;
  return PDMComp_general(drcyl*LengthUnits, zint);
}

double PDMComp2(double zint){
  extern double r2;
  return PDMComp_general(r2, zint);
}

// Will be called by qromb to find the pressure at every point in disk.

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n)
{
  static double s;
  static int it;
  int j;
  double del, sum, tnm, x;

  if (n == 1){
    it = 1;
    return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
  }

  tnm = it;
  del = (b-a)/tnm;
  x = a+0.5*del;
  for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
  it *= 2;
  s = 0.5*(s+(b-a)*sum/tnm);
  return s;
} // end trapezoid

#define K 7  // FIXME
FLOAT polint_c[K+1];
FLOAT polint_d[K+1];

/* also called by qromb */
void polint(double xa[],double ya[],int n,double x,double *y,double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  void nrerror(char *);

  dif=fabs(x-xa[1]);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    polint_c[i]=ya[i];
    polint_d[i]=ya[i];
  } // end i for
  
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=polint_c[i+1]-polint_d[i];
      if ( (den=ho-hp) == 0.0 ) fprintf(stderr,"Error in routine POLINT\n");
      den = w/den;
      polint_d[i]=hp*den;
      polint_c[i]=ho*den;
    } // end i for
    *dy=(2*ns < (n-m) ? polint_c[ns+1] : polint_d[ns--]);
    *y += (*dy);
  } // end m for
} // end polint

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP JMAX+1

/* Main integration routine called by DiskPotentialCircularVelocity to find Pressure */
double qromb(double (*func)(double), double a, double b)
{
  if( a == b ) return 0.0;
  double ss,dss,trapzd(double (*func)(double), double a, double b, int n);
  int j;
  double h[JMAXP+1], s[JMAXP+1];
  void polint(double xa[],double ya[],int n,double x,double *y,double *dy),nrerror(char *);

  h[1] = 1.0;
  for (j=1;j<=JMAX;j++){
    s[j] = trapzd(func,a,b,j);
    if( isnan(s[j]) ) ENZO_FAIL("NaN's during pressure integration in GalaxySimulationInitialize");
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j]; 
  }
  /* Print bug report and exit */
  fprintf(stderr,"Too many steps in routine QROMB\n");
  fprintf(stderr,"\t>> drcyl = %"FSYM", z = %"FSYM", z_icm = %"FSYM"\n", drcyl*LengthUnits/Mpc, a/Mpc, b/Mpc);
  fprintf(stderr,"\t>> ss = %"FSYM", dss = %"FSYM"\n", ss, dss);
  ENZO_FAIL("FAILED IN QROMB IN GRID_GALAXYSIMULATIONINIALIZE\n");
  return -1.0;
}

double bilinear_interp(double x, double y, 
                       double x1, double x2, double y1, double y2,
                       double f_x1y1, double f_x1y2, 
                       double f_x2y1, double f_x2y2) {
double interp;

interp = f_x1y1*(x2-x)*(y2-y) + f_x2y1*(x-x1)*(y2-y) 
       + f_x1y2*(x2-x)*(y-y1) + f_x2y2*(x-x1)*(y-y1);

interp *= 1/( (x2-x1)*(y2-y1) );

return interp;

}
