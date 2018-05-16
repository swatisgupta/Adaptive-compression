#ifndef __ACOMP_COMPRESSORS_H__
#define __ACOMP_COMPRESSORS_H__

#include "util.h"
#define T_COMP_ALGOS 3

typedef 
enum {
	ZLIB = 0,
	LZO,
	BZIP2,
} compressor_t;

#define COMPRESSOR(x)   (x) == ZLIB  ? "ZLIB":   \
			(x) == LZO   ?  "LZO":   \
			(x) == BZIP2 ? "BZIP2":  \
			"NONE" 

acomp_size_t compress_prep (byte_t *in_data, byte_t *comp_buffer, 
		acomp_size_t data_size, compressor_t compression_type );

void decompress_prep (byte_t *out_data, byte_t *comp_buffer,
   		       acomp_size_t comp_size, acomp_size_t data_size,
					compressor_t compression_type );


acomp_size_t compress_prep_by_bytecol( byte_t *in_data, byte_t *comp_buffer,
                                       int n_bytecols, acomp_size_t bytecol_size,
                                              compressor_t compression_type );

void decompress_prep_by_bytecol( byte_t *in_data, byte_t *comp_buffer,
                                       int n_bytecols, acomp_size_t bytecol_size,
                                           compressor_t compression_type );

#endif
