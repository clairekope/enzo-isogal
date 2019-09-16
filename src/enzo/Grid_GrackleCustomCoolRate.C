/***********************************************************************
/
/  GRID CLASS (COMPUTE THE COOLING TIME FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Elizabeth Harper-Clark, August 2009
/              added in CoolingModel parameter
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the cooling time

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
 
/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */
 
extern int RadiationFieldRecomputeMetalRates;
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
int GadgetCoolingTime(float *d, float *e, float *ge, 
		      float *u, float *v, float *w,
		      float *cooltime,
		      int *in, int *jn, int *kn, int *iexpand, 
		      hydro_method *imethod, int *idual, int *idim,
		      int *is, int *js, int *ks, int *ie, int *je, 
		      int *ke, float *dt, float *aye,
		      float *fh, float *utem, float *uxyz, 
		      float *uaye, float *urho, float *utim,
		      float *gamma);

extern "C" void FORTRAN_NAME(cool_multi_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
	float *cooltime, int *coolonly,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
	hydro_method *imethod,
        int *idual, int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co,
	int *ipiht, int *igammah,
	float *dt, float *aye, float *redshift, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *eta1, float *eta2, float *gamma, float *z_solar,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, float *ciHeIa,
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa, float *gammaha,
	float *comp_xraya, float *comp_temp, float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI, float *metal,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
	float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
	float *gasgra, float *metala, int *n_xe, float *xe_start, float *xe_end,
	float *inutot, int *iradfield, int *nfreq, int *imetalregen,
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
	int *iradtrans, float *photogamma,
	int *ih2optical, int *iciecool, float *ciecoa,
 	int *icmbTfloor, int *iClHeat,
 	float *clEleFra, int *clGridRank, int *clGridDim,
 	float *clPar1, float *clPar2, float *clPar3, float *clPar4, float *clPar5,
 	int *clDataSize, float *clCooling, float *clHeating);

extern "C" void FORTRAN_NAME(cool_time)(
	float *d, float *e, float *ge, float *u, float *v, float *w,
           float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
	hydro_method *imethod, int *idual, int *idim, int *igammah,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	float *dt, float *aye, float *temstart, float *temend,
	float *fh, float *utem, float *urho, 
	float *eta1, float *eta2, float *gamma, float *coola, float *gammaha, float *mu);
 
int grid::GrackleCustomCoolRate(int rank, int *dim, float *cool_rate,
				float *dens, float *thrmeng,
				float *velx, float *vely, float *velz,
				float *HIdens, float *HIIdens,
				float *HeIdens, float *HeIIdens, float *HeIIIdens,
				float *edens,
				float *HMdens, float *H2Idens, float *H2IIdens,
				float *DIdens, float *DIIdens, float *HDIdens,
				float *metaldens,
				float *kphHI, float *kphHeI, float *kphHeII,
				float *kdissH2I, float *gamma)
{

  // All passed fields should be in code units

#ifdef USE_GRACKLE 

  /* Return if this doesn't concern us. */

  if (grackle_data->use_grackle == FALSE) return SUCCESS;
  
  if (RadiativeCooling == 0) return SUCCESS;
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i;
  int size = 1;
  for (int d = 0; d < rank; d++)
    size *= dim[d];
 
  /* Compute the cooling time. */
 
  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
 
    aUnits = 1.0/(1.0 + InitialRedshift);
  } else if (RadiationFieldRedshift > -1){
    a       = 1.0 / (1.0 + RadiationFieldRedshift);
    aUnits  = 1.0;
  }
  float afloat = float(a);

  float *volumetric_heating_rate = NULL;
  float *specific_heating_rate   = NULL;

  float *cooling_time = new float[size];
  
  // Double check if there's a metal field when we have metal cooling
  int metal_cooling = MetalCooling;
  if (metal_cooling && !metaldens) {
    if (debug)
      fprintf(stderr, "Warning: No metal field passed to GrackleCustomCoolRate. Not using metal cooling.\n");
    metal_cooling = FALSE;
  }
  printf("Setup Grackle\n"); 
  Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
  g_grid_dimension = new Eint32[rank];
  g_grid_start = new Eint32[rank];
  g_grid_end = new Eint32[rank];
  for (i = 0; i < GridRank; i++) {
    g_grid_dimension[i] = (Eint32) dim[i];
    g_grid_start[i] = (Eint32) 1;
    g_grid_end[i] = (Eint32) dim[i]+1;
  }

  /* Update units. */

  code_units grackle_units;
  grackle_units.comoving_coordinates = (Eint32) ComovingCoordinates;
  grackle_units.density_units        = (double) DensityUnits;
  grackle_units.length_units         = (double) LengthUnits;
  grackle_units.time_units           = (double) TimeUnits;
  grackle_units.velocity_units       = (double) VelocityUnits;
  grackle_units.a_units              = (double) aUnits;
  grackle_units.a_value              = (double) a;

  /* set up the my_fields */
  grackle_field_data my_fields;

  my_fields.grid_rank = (Eint32) GridRank;
  my_fields.grid_dimension = g_grid_dimension;
  my_fields.grid_start     = g_grid_start;
  my_fields.grid_end       = g_grid_end;
  my_fields.grid_dx        = this->CellWidth[0][0]; // CHANGE

  /* now add in the baryon fields */
  my_fields.density         = (gr_float*)dens;
  my_fields.internal_energy = thrmeng;
  my_fields.x_velocity      = velx;
  my_fields.y_velocity      = vely;
  my_fields.z_velocity      = velz;

  if (MultiSpecies) {
    my_fields.HI_density      = HIdens;
    my_fields.HII_density     = HIIdens;
    my_fields.HeI_density     = HeIdens;
    my_fields.HeII_density    = HeIIdens;
    my_fields.HeIII_density   = HeIIIdens;
    my_fields.e_density       = edens;

    if (MultiSpecies > 1) {
      my_fields.HM_density      = HMdens;
      my_fields.H2I_density     = H2Idens;
      my_fields.H2II_density    = H2IIdens;
    
      if (MultiSpecies > 2) {
	my_fields.DI_density      = DIdens;
	my_fields.DII_density     = DIIdens;
	my_fields.HDI_density     = HDIdens;
      }
    }
  }
  
  //my_fields.metal_density   = metaldens;
  
  my_fields.volumetric_heating_rate  = volumetric_heating_rate;
  my_fields.specific_heating_rate    = specific_heating_rate;

#ifdef TRANSFER

  /* unit conversion from Enzo RT units to CGS */
  const float ev2erg = 1.60217653E-12;
  float rtunits = ev2erg / TimeUnits;

  if ( RadiativeTransfer ){
    my_fields.RT_HI_ionization_rate = kphHI;

    if (RadiativeTransferHydrogenOnly == FALSE){
      my_fields.RT_HeI_ionization_rate  = kphHeI;
      my_fields.RT_HeII_ionization_rate = kphHeII;
    }

    if (MultiSpecies > 1)
      my_fields.RT_H2_dissociation_rate = kdissH2I;

    /* need to convert to CGS units */
    for( i = 0; i < size; i++) gamma[i] *= rtunits;

    my_fields.RT_heating_rate = gamma;

  }
#endif // TRANSFER
  printf("Calc cool time\n");
  if (calculate_cooling_time(&grackle_units, &my_fields, cooling_time) == FAIL) {
    ENZO_FAIL("Error in Grackle calculate_cooling_time.\n");
  }
  printf("Calc cool rate\n");
  // Code units
  for (i = 0; i < size; i++) {
    cool_rate[i] = thrmeng[i] / fabs(cooling_time[i]) / dens[i];
  }
    
#ifdef TRANSFER
  if (RadiativeTransfer){
    /* convert the RT units back to Enzo */
    for(i = 0; i < size; i ++) gamma[i] /= rtunits;

  }
#endif // TRANSFER
  printf("Cleanup\n");
  delete [] cooling_time;
  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

#else

    printf("WARNING: Calling GrackleCustomCoolRate but USE_GRACKLE is False!\n");
    
#endif // USE_GRACKLE
    
    return SUCCESS;
}



