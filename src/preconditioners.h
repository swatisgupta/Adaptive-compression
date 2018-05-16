#ifndef __ACOMP_PRECONDITIONERS_H__
#define __ACOMP_PRECONDITIONERS_H__

#include "util.h"

#define T_PREPROCESS 3
typedef 
enum {
  BSEG =  0,
  BWSEG,
  BWXOR,
} preconditioner_t;


#define PRECOND(x)   (x) == BSEG   ? "BSEG":             \
		     (x) == BWSEG  ? "BWSEG":            \
		     (x) == BWXOR  ? "BWXOR":            \
		      "NONE"

void preprocess_data (ubyte_t byteArray, byte_t *in_buffer, byte_t *comp_buffer, byte_t *uncomp_buffer,
                   int ele_size, acomp_size_t data_size, preconditioner_t transform_type );

void postprocess_data (ubyte_t byteArray, byte_t *in_buffer, byte_t *comp_buffer, byte_t *uncomp_buffer,
                   int ele_size, acomp_size_t data_size, preconditioner_t transform_type ); 

#endif 
