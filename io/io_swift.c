/*
 * Swift I/O for Rockstar
 * Fridolin Glatter (fridolin.glatter@univie.ac.at) and Lukas Winkler (winklerl23@univie.ac.at)
 */
#ifdef ENABLE_HDF5

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <hdf5.h> /* HDF5 required */
#include "io_hdf5.h"
#include "io_swift.h"
#include "io_util.h"
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"

#define SWIFT_NTYPES 7

void swift_read_dataset(hid_t HDF_FileID, char *filename, char *gid, char *dataid, struct particle *p, int64_t to_read, int64_t offset, int64_t stride, hid_t type) {
  int64_t width = (type == H5T_NATIVE_LLONG) ? 8 : 4;
  void *buffer = check_malloc_s(buffer, to_read, width*stride);
  int64_t *ibuffer = buffer;
  float *fbuffer = buffer;

  hid_t HDF_GroupID = check_H5Gopen(HDF_FileID, gid, filename);
  hid_t HDF_DatasetID = check_H5Dopen(HDF_GroupID, dataid, gid, filename);
  hid_t HDF_DataspaceID = check_H5Dget_space(HDF_DatasetID);

  check_H5Sselect_all(HDF_DataspaceID);
  hssize_t npoints = H5Sget_select_npoints(HDF_DataspaceID);

  if (npoints != to_read*stride) {
    fprintf(stderr, "[Error] dataspace %s/%s in HDF5 file %s not expected size!\n  (Actual size = %"PRId64" elements; expected size = %"PRId64" elements\n", 
	    gid, dataid, filename, (int64_t)(npoints), stride*to_read);
    exit(1);
  }

  check_H5Dread(HDF_DatasetID, type, buffer, dataid, gid, filename);
  
  H5Sclose(HDF_DataspaceID);
  H5Dclose(HDF_DatasetID);
  H5Gclose(HDF_GroupID);

  if (width == 8)
    for (int64_t i=0; i<to_read; i++)
      p[i].id = ibuffer[i];
  else
    for (int64_t i=0; i<to_read; i++)
      memcpy(((char *)&(p[i]))+offset, fbuffer+(i*stride), stride*width);

  free(buffer);
}

float swift_readheader_float(hid_t HDF_GroupID, char *filename, char *objName)
{
  char *gid = "Header";
  hid_t HDF_Type = H5T_NATIVE_FLOAT;
  hid_t HDF_AttrID = check_H5Aopen_name(HDF_GroupID, objName, gid, filename);
  hid_t HDF_DataspaceID = check_H5Aget_space(HDF_AttrID);

  check_H5Sselect_all(HDF_DataspaceID);
  
  float data = 0.0;
  check_H5Aread( HDF_AttrID, HDF_Type, &data, objName, gid, filename);

  H5Sclose(HDF_DataspaceID);
  H5Aclose(HDF_AttrID);
  return data;
}


void swift_readheader_array(hid_t HDF_GroupID, char *filename, char *objName, hid_t type, void *data)
{
  char *gid = "Header";
  hid_t HDF_AttrID = check_H5Aopen_name(HDF_GroupID, objName, gid, filename);
  hid_t HDF_DataspaceID = check_H5Aget_space(HDF_AttrID);
  check_H5Sselect_all(HDF_DataspaceID);

  int64_t ndims = check_H5Sget_simple_extent_ndims( HDF_DataspaceID );
  assert(ndims == 1);
  hsize_t dimsize = 0;
  check_H5Sget_simple_extent_dims(HDF_DataspaceID, &dimsize);
  assert(dimsize == SWIFT_NTYPES);
  
  check_H5Aread(HDF_AttrID, type, data, objName, gid, filename);

  H5Sclose(HDF_DataspaceID);
  H5Aclose(HDF_AttrID);
}

void swift_rescale_particles(struct particle *p, int64_t p_start, int64_t nelems) {
  for (int64_t i=0; i<nelems; i++) {
    for (int64_t j=0; j<3; j++) {
      p[p_start+i].pos[j]   *= h0;
    }
  }
}

void load_particles_swift(char *filename, struct particle **p, int64_t *num_p)
{	
  hid_t HDF_FileID = check_H5Fopen(filename, H5F_ACC_RDONLY);
  hid_t HDF_Header = check_H5Gopen(HDF_FileID, "Header", filename);
  hid_t HDF_Cosmology = check_H5Gopen(HDF_FileID, "Cosmology", filename);
  
  Ol = swift_readheader_float(HDF_Cosmology, filename, "Omega_lambda");
  Om = swift_readheader_float(HDF_Cosmology, filename, "Omega_b") + swift_readheader_float(HDF_Cosmology, filename, "Omega_cdm");
  h0 = swift_readheader_float(HDF_Cosmology, filename, "h");
  SCALE_NOW = swift_readheader_float(HDF_Cosmology, filename, "Scale-factor");
  BOX_SIZE = swift_readheader_float(HDF_Header, filename, "BoxSize");
  BOX_SIZE *= SWIFT_LENGTH_CONVERSION;
  BOX_SIZE *= h0;  //Rockstar wants Mpc/h
  
  uint32_t npart_low[SWIFT_NTYPES], npart_high[SWIFT_NTYPES] = {0};
  int64_t npart[SWIFT_NTYPES];
  float massTable[SWIFT_NTYPES];
  
  swift_readheader_array(HDF_Header, filename, "NumPart_ThisFile", H5T_NATIVE_UINT64, npart);
  swift_readheader_array(HDF_Header, filename, "NumPart_Total_HighWord", H5T_NATIVE_UINT32, npart_high);
  swift_readheader_array(HDF_Header, filename, "NumPart_Total", H5T_NATIVE_UINT32, npart_low);
  swift_readheader_array(HDF_Header, filename, "InitialMassTable", H5T_NATIVE_FLOAT, massTable);

   for (int i = 0; i < SWIFT_NTYPES; i++) {
       if (npart_low[i] > 0) {
           char parttype[10];
           sprintf(parttype, "PartType%d", i);
           printf("AAAAAAAA: %s\n", parttype);
           hid_t HDF_GroupID = check_H5Gopen(HDF_FileID, parttype, filename);
           hid_t HDF_DatasetID = check_H5Dopen(HDF_GroupID, "Masses", parttype, filename);
           printf("opened\n");
           hsize_t coord[1];
           coord[0] = 0; //we should be able to select the first mass since all are equal
           hid_t dataspace_masses = check_H5Dget_space(HDF_DatasetID);
           herr_t status = H5Sselect_elements(dataspace_masses, H5S_SELECT_SET, 1, (const hsize_t *)&coord);

           hsize_t number_of_points = 1;
           float mass[number_of_points];
           hid_t helpMemSpace = H5Screate_simple(1, &number_of_points, NULL);
           herr_t status = H5Dread(HDF_Dataset_ID, H5T_NATIVE_FLOAT, helpMemSpace, dataspace_masses, H5P_DEFAULT, mass)

// FIXME: the next line just hangs. Probably reading just the first value doesn't work
          //  check_H5Dread(HDF_DatasetID, H5T_NATIVE_FLOAT, mass, "Masses", parttype, filename);
           printf("loaded\n");

           massTable[i] = mass[0] * h0; //This will give Msun/h in the end
           H5Gclose(HDF_GroupID);
           H5Dclose(HDF_DatasetID);
           printf("particles finished\n");

       } else {
           massTable[i] = 0.0;
       }
   }


  TOTAL_PARTICLES = ( ((int64_t)npart_high[SWIFT_DM_PARTTYPE]) << 32 )
    + (int64_t)npart_low[SWIFT_DM_PARTTYPE];
  
  H5Gclose(HDF_Header);
  H5Gclose(HDF_Cosmology);
    
  PARTICLE_MASS   = massTable[SWIFT_DM_PARTTYPE] * SWIFT_MASS_CONVERSION;
  AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));
	
  if(RESCALE_PARTICLE_MASS)
    PARTICLE_MASS = Om*CRITICAL_DENSITY * pow(BOX_SIZE, 3) / TOTAL_PARTICLES;
 
  printf("SWIFT: filename:       %s\n", filename);
  printf("SWIFT: box size:       %g Mpc/h\n", BOX_SIZE);
  printf("SWIFT: h0:             %g\n", h0);
  printf("SWIFT: scale factor:   %g\n", SCALE_NOW);
  printf("SWIFT: Total DM Part:  %" PRIu64 "\n", TOTAL_PARTICLES);
  printf("SWIFT: ThisFile DM Part: %" PRIu64 "\n", npart[SWIFT_DM_PARTTYPE]);
  printf("SWIFT: DM Part Mass:   %g Msun/h\n", PARTICLE_MASS);
  printf("SWIFT: avgPartSpacing: %g Mpc/h\n\n", AVG_PARTICLE_SPACING);
  
  if (!npart[SWIFT_DM_PARTTYPE]) {
    H5Fclose(HDF_FileID);
    printf("   SKIPPING FILE, PARTICLE COUNT ZERO.\n");
    return;
  }

  int64_t to_read = npart[SWIFT_DM_PARTTYPE];
  check_realloc_s(*p, ((*num_p)+to_read), sizeof(struct particle));

  // read IDs, pos, vel
  char buffer[100];
  snprintf(buffer, 100, "PartType%"PRId64, SWIFT_DM_PARTTYPE);
  swift_read_dataset(HDF_FileID, filename, buffer, "ParticleIDs", *p + (*num_p),
	 to_read, (char *)&(p[0][0].id)-(char*)(p[0]), 1, H5T_NATIVE_LLONG);
  swift_read_dataset(HDF_FileID, filename, buffer, "Coordinates", *p + (*num_p),
	 to_read, (char *)&(p[0][0].pos[0])-(char*)(p[0]), 3, H5T_NATIVE_FLOAT);
  swift_read_dataset(HDF_FileID, filename, buffer, "Velocities", *p + (*num_p),
	 to_read, (char *)&(p[0][0].pos[3])-(char*)(p[0]), 3, H5T_NATIVE_FLOAT);

  H5Fclose(HDF_FileID);
  
  swift_rescale_particles(*p, *num_p, to_read); //only necessary to convert to Mpc/h
  
  *num_p += npart[SWIFT_DM_PARTTYPE];
}

#endif /* ENABLE_HDF5 */
