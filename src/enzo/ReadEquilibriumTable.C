/***********************************************************************
/
/
************************************************************************/

#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/**************************** Functions Prototypes ******************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

//herr_t read_chem_dsets(hid_t orig_id, const char* orig_name, 
//                        const H5O_info_t* orig_info, void* chemistry_storage);

// Read Equilibrium Table
int ReadEquilibriumTable(const char * name, FLOAT Time)
{

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	            &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
  /* Read data in from hdf5 file.*/

  hid_t       file_id, grp_id, dset_id, attr_id; 
  herr_t      status;
  herr_t      h5_error = -1;   

  if (debug) fprintf(stderr,"Reading from %s.\n",name);
  file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
  grp_id = H5Gopen(file_id, "data");
  if (grp_id == h5_error) {
    fprintf(stderr, "Can't open data group in %s.\n",name);
  }

  /* get the size (in one dimension) of the data */

  long_int* size = new long_int[1];
  attr_id = H5Aopen_name(grp_id, "num_elements");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open num_elements attribute in %s.\n",name);
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8, size);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to read Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  EquilibriumTable.dim_size = *size; // pointer to int
  if (debug) fprintf(stderr,"Equilibrium table size (num_elements) is %d\n",
                      EquilibriumTable.dim_size);
  delete [] size;

  /* what chem species are in the group? */

  /*H5Ovisit(grp_id, H5_INDEX_NAME, H5_ITER_NATIVE, read_chem_dsets, 
                      (void *) &EquilibriumTable);*/

  status = H5Gclose(grp_id);
  if (status == h5_error) {
    fprintf(stderr,"Failed to close group in %s.\n",name);
    return FAIL;
  }

  if (MultiSpecies) {
    EquilibriumTable.HI = new double[EquilibriumTable.dim_size
                                      * EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/data/HI");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /data/HI in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.HI);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /data/HI in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /data/HI in %s.\n",name);
      return FAIL;
    }

    EquilibriumTable.HII = new double[EquilibriumTable.dim_size
                                      * EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/data/HII");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /data/HII in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.HII);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /data/HII in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /data/HII in %s.\n",name);
      return FAIL;
    }

    EquilibriumTable.HeI = new double[EquilibriumTable.dim_size
                                      * EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/data/HeI");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /data/HeI in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.HeI);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /data/HeI in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /data/HeI in %s.\n",name);
      return FAIL;
    }

    EquilibriumTable.HeII = new double[EquilibriumTable.dim_size
                                      * EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/data/HeII");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /data/HeII in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.HeII);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /data/HeII in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /data/HeII in %s.\n",name);
      return FAIL;
    }

    EquilibriumTable.HeIII = new double[EquilibriumTable.dim_size
                                      * EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/data/HeIII");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /data/HeIII in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.HeIII);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /data/HeIII in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /data/HeIII in %s.\n",name);
      return FAIL;
    }

    EquilibriumTable.de = new double[EquilibriumTable.dim_size
                                      * EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/data/de");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /data/de in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.de);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /data/de in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /data/de in %s.\n",name);
      return FAIL;
    }

    if (MultiSpecies > 1) {
      EquilibriumTable.HM = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/data/HM");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /data/HM in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.HM);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /data/HM in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /data/HM in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.H2I = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/data/H2I");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /data/H2I in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.H2I);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /data/H2I in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /data/H2I in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.H2II = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/data/H2II");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /data/H2II in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.H2II);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /data/H2II in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /data/H2II in %s.\n",name);
        return FAIL;
      }
    }
    if (MultiSpecies > 2) {
      EquilibriumTable.DI = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/data/DI");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /data/DI in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.DI);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /data/DI in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /data/DI in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.DII = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/data/DII");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /data/DII in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.DII);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /data/DII in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /data/DII in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.HDI = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/data/HDI");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /data/HDI in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.HDI);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /data/HDI in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /data/HDI in %s.\n",name);
        return FAIL;
      }
    }
  }

  /* Read in the density and temperature */

  EquilibriumTable.density = new double[EquilibriumTable.dim_size];
  dset_id = H5Dopen(file_id, "/data/density");
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open /data/density in %s.\n", name);
    return FAIL;
  }
  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                    H5S_ALL, H5P_DEFAULT, EquilibriumTable.density);
  if (status == h5_error) {
    fprintf(stderr, "Failed to read /data/density in %s.\n",name);
    return FAIL;
  }
  status = H5Dclose(dset_id);
  if (status == h5_error) {
    fprintf(stderr,"Failed to close /data/density in %s.\n",name);
    return FAIL;
  }

  EquilibriumTable.temperature = new double[EquilibriumTable.dim_size];
  dset_id = H5Dopen(file_id, "/data/temperature");
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open /data/temperature in %s.\n", name);
    return FAIL;
  }
  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                    H5S_ALL, H5P_DEFAULT, EquilibriumTable.temperature);
  if (status == h5_error) {
    fprintf(stderr, "Failed to read /data/temperature in %s.\n",name);
    return FAIL;
  }
  status = H5Dclose(dset_id);
  if (status == h5_error) {
    fprintf(stderr,"Failed to close /data/temperature in %s.\n",name);
    return FAIL;
  }

  /* convert from cgs to code units */

  for (int i=0; i<EquilibriumTable.dim_size; ++i){
    EquilibriumTable.density[i] /= DensityUnits;
    //EquilibriumTable.temperature[i] /= TemperatureUnits; // keep in K
    if (MultiSpecies) {
      EquilibriumTable.HI[i] /= DensityUnits;
      EquilibriumTable.HII[i] /= DensityUnits;
      EquilibriumTable.HeI[i] /= DensityUnits;
      EquilibriumTable.HeII[i] /= DensityUnits;
      EquilibriumTable.HeIII[i] /= DensityUnits;
      EquilibriumTable.de[i] /= DensityUnits;
      if (MultiSpecies > 1) {
        EquilibriumTable.HM[i] /= DensityUnits;
        EquilibriumTable.H2I[i] /= DensityUnits;
        EquilibriumTable.H2II[i] /= DensityUnits;
      }
      if (MultiSpecies > 2) {
        EquilibriumTable.DI[i] /= DensityUnits;
        EquilibriumTable.DII[i] /= DensityUnits;
        EquilibriumTable.HDI[i] /= DensityUnits;
      }
    }
  }

  status = H5Fclose (file_id);
  
  return SUCCESS;
}

//herr_t read_chem_dsets(hid_t locm_id, const char* name, 
//                        const H5O_info_t* info, void* chemistry_storage){
    //printf ("/");               /* Print root group in object path */

    /*
     * Check if the current object is the root group, and if not print
     * the full path name and type.
     */
//    if (name[0] == '.')         /* Root group, do not print '.' */
/*        printf ("  (Group)\n");
    else
        switch (info->type) {
            case H5O_TYPE_GROUP:
                printf ("%s  (Group)\n", name);
                break;
            case H5O_TYPE_DATASET:
                printf ("%s  (Dataset)\n", name);
                break;
            case H5O_TYPE_NAMED_DATATYPE:
                printf ("%s  (Datatype)\n", name);
                break;
            default:
                printf ("%s  (Unknown)\n", name);
        }

    return 0;
}*/
