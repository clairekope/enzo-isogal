/*------------------------------------------------------------------------
  SET EVOLVING STELLAR MASS THRESHOLD
  By Brian O'Shea

  History:
     23 April 2019 : BWO -- Created

  Note:
     StarMakerThermalFeedbackRamp = 0 is off
     StarMakerThermalFeedbackRamp = 1 is linear evolution of mass in time
     StarMakerThermalFeedbackRamp = 2 is linear evolution of mass in redshift
     StarMakerThermalFeedbackRamp = 3 is exponential evolution of mass in time
     StarMakerThermalFeedbackRamp = 4 is exponential evolution of mass in redshift

  If StarMakerThermalFeedbackRamp > 0, all values of the ramp parameters (starting and
  ending times and masses) are set in ReadParameterFile.C, and tests are made there 
  to ensure that the user has set them.

------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"

void my_exit(int status);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int SetStellarFeedbackEfficiency(FLOAT time)
{

  int timestep, i;
  FLOAT a, dadt, redshift=0.0;
  float early_fbeff, late_fbeff, current_fbeff, float_time=0.0, float_redshift=0.0;

  /* Return if not used */
  if (StarMakerThermalFeedbackRamp == 0)
    return SUCCESS;

  /* Calculate redshift */
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    redshift = (1 + InitialRedshift)/a - 1;
  }

  /* recast redshift and time to be the same precision as everything else; only important if
     the two different floating-point precisions are different. */
  float_redshift = (float) redshift;
  float_time = (float) time;

  if(StarMakerThermalFeedbackRamp == 1 || StarMakerThermalFeedbackRamp == 3){  // interpolation in time

    /* Set early and late efficiencies in linear or log */
    if(StarMakerThermalFeedbackRamp == 1){ // mass evolution linear in time
      early_fbeff = StarMakerThermalFeedbackRampStartValue;
      late_fbeff = StarMakerThermalFeedbackRampEndValue;
    } else { // mass evolution exponential in time
      early_fbeff = log10(StarMakerThermalFeedbackRampStartValue);
      late_fbeff = log10(StarMakerThermalFeedbackRampEndValue);
    }

    /* set current stellar feedback efficiency */
    if(float_time <= StarMakerThermalFeedbackRampStartTime){ // if time is before ramp start time, use early mass
      current_fbeff = early_fbeff;
    } else if (float_time >= StarMakerThermalFeedbackRampEndTime){ // if time is after ramp end time, use late mass
      current_fbeff = late_fbeff;
    } else {  // otherwise, linearly interpolate between start and end
      current_fbeff = early_fbeff + (float_time - StarMakerThermalFeedbackRampStartTime)
	* (late_fbeff-early_fbeff)/(StarMakerThermalFeedbackRampEndTime-StarMakerThermalFeedbackRampStartTime);  
    }

    /* set StarEnergyToThermalFeedback correctly */
    if(StarMakerThermalFeedbackRamp == 1){
      StarEnergyToThermalFeedback = current_fbeff;
    } else {
      StarEnergyToThermalFeedback = POW(10.0,current_fbeff);
    }

  } else if(StarMakerThermalFeedbackRamp == 2 || StarMakerThermalFeedbackRamp == 4){  // interpolation in redshift

    /* set early and late efficiencies in linear or log */
    if(StarMakerThermalFeedbackRamp == 2){ // mass evolution linear in redshift
      early_fbeff = StarMakerThermalFeedbackRampStartValue;
      late_fbeff = StarMakerThermalFeedbackRampEndValue;
    } else { // efficiency evolution exponential in time
      early_fbeff = log10(StarMakerThermalFeedbackRampStartValue);
      late_fbeff = log10(StarMakerThermalFeedbackRampEndValue);
    }

    /* set current stellar feedback efficiency */
    if(float_redshift >= StarMakerThermalFeedbackRampStartTime){ // if redshift is prior to ramp start redshift, use early mass
      current_fbeff = early_fbeff;
    } else if (float_redshift <= StarMakerThermalFeedbackRampEndTime){ // if redshift is after ramp end redshift, use late mass
      current_fbeff = late_fbeff;
    } else {  // otherwise, linearly interpolate between start and end
      current_fbeff = early_fbeff + (float_redshift - StarMakerThermalFeedbackRampStartTime)
	* (late_fbeff-early_fbeff)/(StarMakerThermalFeedbackRampEndTime-StarMakerThermalFeedbackRampStartTime);
    }

    /* set StarEnergyToThermalFeedback correctly */
    if(StarMakerThermalFeedbackRamp == 2){
      StarEnergyToThermalFeedback = current_fbeff;
    } else {
      StarEnergyToThermalFeedback = POW(10.0,current_fbeff);
    }

  } else {  // user has made a poor choice
    fprintf(stderr,"SetStellarFeedbackEfficiency:  StarMakerThermalFeedbackRamp improperly set!\n");
    my_exit(EXIT_FAILURE);
  }
  
  if(debug){
    printf("SetStellarFeedbackEfficiency:  StarEnergyToThermalFeedback set to %"FSYM" at time %"PSYM" (redshift %"PSYM")\n",
	   StarEnergyToThermalFeedback, time, redshift);
  }
  
  return SUCCESS;

}
