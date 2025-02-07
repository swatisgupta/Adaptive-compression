/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write a global array from N processors with gwrite
 *
 * How to run: mpirun -np <N> global_array_aggregate_by_color_C
 * Output: global_array_aggregate_by_color_C.bp
 * ADIOS config file: None
 *
*/

/* This example will write out 2 sub blocks of the variable temperature
   and place these in the global array.
   This example illustrates both the use of sub blocks in writing, and
   the usage of the ADIOS non-xml API's
*/

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "public/adios.h"
#include "public/adios_types.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

int main (int argc, char ** argv) 
{
	char        filename [256];
        char        color_str[256];
	int         rank, size, i, color;
	int         NX = 100, Global_bounds, Offsets; 
	double      t[NX];
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t    adios_groupsize, adios_totalsize;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

        Global_bounds = NX * size;

	strcpy (filename, "global_array_aggregate_by_color_C.bp");

	adios_init_noxml (comm);
        adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);

        int64_t       m_adios_group;
        int64_t       m_adios_file;

        adios_declare_group (&m_adios_group, "restart", "iter", adios_flag_yes);

        // split into 2 groups 
        color = (rank % 2 == 0 ? 0 : 1);
        sprintf (color_str, "color=%d", color);
        adios_select_method (m_adios_group, "MPI", color_str, "");

        adios_define_var (m_adios_group, "NX"
			,"", adios_integer
			,0, 0, 0);
   
	adios_define_var (m_adios_group, "Global_bounds"
			,"", adios_integer
			,0, 0, 0);

        adios_define_var (m_adios_group, "Offsets"
                         ,"", adios_integer
                         ,0, 0, 0);
   
       adios_define_var (m_adios_group, "temperature"
                        ,"", adios_double
                        ,"NX", "Global_bounds", "Offsets");
   
        adios_open (&m_adios_file, "restart", filename, "w", comm);

        adios_groupsize = 4 + 4 + 4 + NX * 8;

        adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);
	adios_write(m_adios_file, "NX", (void *) &NX);
	adios_write(m_adios_file, "Global_bounds", (void *) &Global_bounds);

        Offsets = rank * NX;
        adios_write(m_adios_file, "Offsets", (void *) &Offsets);

        for (i = 0; i < NX; i++)
            t[i] = Offsets + i;

        adios_write(m_adios_file, "temperature", t);

        adios_close (m_adios_file);

        MPI_Barrier (comm);

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
