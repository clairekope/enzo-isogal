/*------------------------------------------------------------------------
  SET EVOLVING STELLAR MASS THRESHOLD
  By Brian O'Shea

  History:
     23 April 2019 : BWO -- Created

  Note:
     StarFeedbackThermalEfficiencyRamp = 0 is off
     StarFeedbackThermalEfficiencyRamp = 1 is linear evolution of mass in time
     StarFeedbackThermalEfficiencyRamp = 2 is linear evolution of mass in redshift
     StarFeedbackThermalEfficiencyRamp = 3 is exponential evolution of mass in time
     StarFeedbackThermalEfficiencyRamp = 4 is exponential evolution of mass in redshift

  If StarFeedbackThermalEfficiencyRamp > 0, all values of the ramp parameters (starting and
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
  if (StarFeedbackThermalEfficiencyRamp == 0)
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

  if(StarFeedbackThermalEfficiencyRamp == 1 || StarFeedbackThermalEfficiencyRamp == 3){  // interpolation in time

    /* Set early and late efficiencies in linear or log */
    if(StarFeedbackThermalEfficiencyRamp == 1){ // mass evolution linear in time
      early_fbeff = StarFeedbackThermalEfficiencyRampStartValue;
      late_fbeff = StarFeedbackThermalEfficiencyRampEndValue;
    } else { // mass evolution exponential in time
      early_fbeff = log10(StarFeedbackThermalEfficiencyRampStartValue);
      late_fbeff = log10(StarFeedbackThermalEfficiencyRampEndValue);
    }

    /* set current stellar feedback efficiency */
    if(float_time <= StarFeedbackThermalEfficiencyRampStartTime){ // if time is before ramp start time, use early mass
      current_fbeff = early_fbeff;
    } else if (float_time >= StarFeedbackThermalEfficiencyRampEndTime){ // if time is after ramp end time, use late mass
      current_fbeff = late_fbeff;
    } else {  // otherwise, linearly interpolate between start and end
      current_fbeff = early_fbeff + (float_time - StarFeedbackThermalEfficiencyRampStartTime)
	* (late_fbeff-early_fbeff)/(StarFeedbackThermalEfficiencyRampEndTime-StarFeedbackThermalEfficiencyRampStartTime);  
    }

    /* set StarEnergyToThermalFeedback correctly */
    if(StarFeedbackThermalEfficiencyRamp == 1){
      StarEnergyToThermalFeedback = current_fbeff;
    } else {
      StarEnergyToThermalFeedback = POW(10.0,current_fbeff);
    }

  } else if(StarFeedbackThermalEfficiencyRamp == 2 || StarFeedbackThermalEfficiencyRamp == 4){  // interpolation in redshift

    /* set early and late efficiencies in linear or log */
    if(StarFeedbackThermalEfficiencyRamp == 2){ // mass evolution linear in redshift
      early_fbeff = StarFeedbackThermalEfficiencyRampStartValue;
      late_fbeff = StarFeedbackThermalEfficiencyRampEndValue;
    } else { // efficiency evolution exponential in time
      early_fbeff = log10(StarFeedbackThermalEfficiencyRampStartValue);
      late_fbeff = log10(StarFeedbackThermalEfficiencyRampEndValue);
    }

    /* set current stellar feedback efficiency */
    if(float_redshift >= StarFeedbackThermalEfficiencyRampStartTime){ // if redshift is prior to ramp start redshift, use early mass
      current_fbeff = early_fbeff;
    } else if (float_redshift <= StarFeedbackThermalEfficiencyRampEndTime){ // if redshift is after ramp end redshift, use late mass
      current_fbeff = late_fbeff;
    } else {  // otherwise, linearly interpolate between start and end
      current_fbeff = early_fbeff + (float_redshift - StarFeedbackThermalEfficiencyRampStartTime)
	* (late_fbeff-early_fbeff)/(StarFeedbackThermalEfficiencyRampEndTime-StarFeedbackThermalEfficiencyRampStartTime);
    }

    /* set StarEnergyToThermalFeedback correctly */
    if(StarFeedbackThermalEfficiencyRamp == 2){
      StarEnergyToThermalFeedback = current_fbeff;
    } else {
      StarEnergyToThermalFeedback = POW(10.0,current_fbeff);
    }

  } else {  // user has made a poor choice
    fprintf(stderr,"SetStellarFeedbackEfficiency:  StarFeedbackThermalEfficiencyRamp improperly set!\n");
    my_exit(EXIT_FAILURE);
  }
  
  if(debug){
    printf("SetStellarFeedbackEfficiency:  StarEnergyToThermalFeedback set to %"FSYM" at time %"PSYM" (redshift %"PSYM")\n",
	   StarEnergyToThermalFeedback, time, redshift);
  }
  
  return SUCCESS;

}
