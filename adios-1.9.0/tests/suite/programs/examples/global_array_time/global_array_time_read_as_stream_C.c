/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 * which has multiple timesteps,
 * reading step by step
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"
#include "adios_error.h"
#include "core/adios_logger.h"

int main (int argc, char ** argv) 
{
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    ADIOS_FILE * f;
    ADIOS_VARINFO * v;
    ADIOS_SELECTION * sel;
    int steps = 0;
    int retval = 0;
    float timeout_sec = 0.0; 

    void * data = NULL;
    uint64_t start[2], count[2];

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=3");
    adios_logger_open ("log_read_as_stream_C", rank);

    f = adios_read_open ("global_array_time_C.bp", ADIOS_READ_METHOD_BP,
                          comm, ADIOS_LOCKMODE_NONE, timeout_sec);
    if (adios_errno == err_file_not_found)
    {
        log_error ("rank %d: Stream not found after waiting %f seconds: %s\n",
                rank, timeout_sec, adios_errmsg());
        retval = adios_errno;
    }
    else if (adios_errno == err_end_of_stream)
    {
        log_error ("rank %d: Stream terminated before open. %s\n", rank, adios_errmsg());
        retval = adios_errno;
    }
    else if (f == NULL) {
        error ("rank %d: Error at opening stream: %s\n", rank, adios_errmsg());
        retval = adios_errno;
    }
    else
    {
        /* process file here... */
        v = adios_inq_var (f, "temperature");
        adios_inq_var_blockinfo (f, v);

        log_test ("rank %d: ndim = %d\n", rank, v->ndim);
        //log_test ("nsteps = %d\n",  v->nsteps);
        log_test ("rank %d: dims[%llu][%llu]\n", rank, v->dims[0], v->dims[1]);

        uint64_t slice_size = v->dims[0]/size;
        if (rank == size-1)
            slice_size = slice_size + v->dims[0]%size;

        start[0] = rank * slice_size;
        count[0] = slice_size;
        start[1] = 0;
        count[1] = v->dims[1];

        data = malloc (slice_size * v->dims[1] * 8);

        /* Processing loop over the steps (we are already in the first one) */
        while (adios_errno != err_end_of_stream) {
            steps++; // steps start counting from 1

            sel = adios_selection_boundingbox (v->ndim, start, count);
            adios_schedule_read (f, sel, "temperature", 0, 1, data);
            adios_perform_reads (f, 1);

            log_test("rank=%d: step %d: [0:%lld,0:%lld] = [", rank, f->current_step, v->dims[0], v->dims[1]);
            for (i = 0; i < slice_size; i++) {
                log_test (" [");
                for (j = 0; j < v->dims[1]; j++) {
                    log_test ("%g ", *((double *)data + i * v->dims[1] + j));
                }
                log_test ("]");
            }
            log_test (" ]\n");

            // advance to 1) next available step with 2) blocking wait
            adios_advance_step (f, 0, timeout_sec);
            if (adios_errno == err_step_notready)
            {
                log_test ("rank %d: No new step arrived within the timeout. Quit. %s\n",
                        rank, adios_errmsg());
                break; // quit while loop
            }

        }

        adios_read_close (f);
    }

    log_test ("rank: %d: We have processed %d steps\n", rank, steps);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    adios_logger_close();
    free (data);
    MPI_Finalize ();

    return retval;
}


