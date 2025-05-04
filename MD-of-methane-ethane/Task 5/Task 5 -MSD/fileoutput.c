#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

// Write the particle positions to a pdf file
// The filename (without extension) is given by p_parameters->filename_pdb.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_pdb(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time) {
    FILE *fp_traj;
    char filename[1024];
    double rs = p_parameters->rescale_output;

    snprintf(filename, 1024, "%s%s", p_parameters->filename_pdb, ".pdb");
    if (reset == 1) {
        fp_traj = fopen(filename, "w");
    } else {
        fp_traj = fopen(filename, "a");
    }

    fprintf(fp_traj, "MODEL\n");
    fprintf(fp_traj, "REMARK TIME = %f\n", time);
    fprintf(fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n", rs*p_parameters->L.x, rs*p_parameters->L.y, rs*p_parameters->L.z, 90.0, 90.0, 90.0, "P 1", "1");
    for (size_t i = 0; i < p_parameters->num_part; i++) {
        char atom_type = (p_vectors->type[i] == 0) ? 'M' : 'E'; // M for Methane, E for Ethane
        fprintf(fp_traj, "HETATM %5u  %c 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", (unsigned int)i % 100000, atom_type, rs*p_vectors->r[i].x, rs*p_vectors->r[i].y, rs*p_vectors->r[i].z);
    }
    fprintf(fp_traj, "ENDMDL\n");

    fclose(fp_traj);
}

// Write the particle positions to a xyz file
// The filename (without extension) is given by p_parameters->filename_xyz.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1024];
  double rs = p_parameters->rescale_output;

  snprintf(filename, 1024, "%s%s", p_parameters->filename_xyz, ".xyz");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "%lu\n", p_parameters->num_part);
  fprintf(fp_traj, "time = %f\n", time);
  struct Vec3D *r = p_vectors->r;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "  C        %10.5f %10.5f %10.5f\n", rs*r[i].x, rs*r[i].y, rs*r[i].z);
  }

  fclose(fp_traj);
}

// save arrays in vectors to binary file
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
  FILE* p_file = fopen( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);

  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

// load arrays in vectors from binary file
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
  FILE* p_file = fopen( p_parameters->restart_in_filename, "rb" );
  size_t num_part;
  fread(&num_part, sizeof(size_t), 1, p_file);
  size_t sz = num_part*sizeof(struct Vec3D);
  alloc_vectors(p_vectors,num_part);
  p_parameters->num_part = num_part;
  fread(p_vectors->r, sz, 1, p_file);
  fread(p_vectors->v, sz, 1, p_file);
  fread(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

// Function to write energy data to a file
void write_energy_to_file(const char *filename, unsigned long step, double kinetic_energy, double potential_energy, double total_energy, double temp) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening energy file");
        exit(EXIT_FAILURE);
    }
    
    if (step == 1) {
    fprintf(file, "%20s %20s %20s %20s %20s\n",
            "         Time Step", "    Kinetic Energy eV", "    Potential Energy eV", "      Total Energy eV", "        Current Temp");
    }

    fprintf(file, "%20lu %20.7f %20.7f %20.7f %20.7f\n",
                    step, kinetic_energy, potential_energy, total_energy, temp);

    fclose(file);
}

// Function to write MSD data to a file
void write_msd_data(const char *msd_output, unsigned long step, double msd[], double diffusion_coefficient[]) {
    // Open MSD data file, write headers if file is being created
    FILE *msd_file = fopen(msd_output, "a");
    if (msd_file != NULL) {
        // Check if the file is empty to write headers
        fseek(msd_file, 0, SEEK_END);
        long size = ftell(msd_file);
        
        if (size == 0) {
            // If the file is empty, write the headers
            fprintf(msd_file, "Time Step, Methane MSD (Å^2), Methane Diffusion Coefficient (Å^2/ps), Ethane MSD (Å^2), Ethane Diffusion Coefficient (Å^2/ps)\n");
        }

        // Return to the end of the file to append data
        fseek(msd_file, 0, SEEK_END);
        fprintf(msd_file, "%lu, %f, %f, %f, %f\n", step, msd[0], diffusion_coefficient[0], msd[1], diffusion_coefficient[1]);

        // Close the file
        fclose(msd_file);
    } else {
        // Error handling if the file could not be opened
        perror("Failed to open MSD output file");
    }
}

void removeFile(const char *fileName) {
    if (remove(fileName) == 0) {
        printf("%s was there and it has been removed\n", fileName);
    } else {
        printf("There was no %s file\n", fileName);
    }
}