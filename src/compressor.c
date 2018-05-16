#include "compressor.h"
#include <zlib.h> 
#include <bzlib.h>
#include <lzo/lzoconf.h>
#include <lzo/lzo1x.h>

static int init_lzo = 0;

void bzip_err (int err);
void zlib_err (int err);
void lzo_err (int err);


acomp_size_t _zlib_compression( byte_t *in_data, byte_t *comp_buffer, acomp_size_t data_size ) {
	
    uLongf csize  =  data_size > 4* 1024? 4*data_size:4* 1024;
    DB_PRINTF("Size on entry %llu\n",csize);	
    int err  =  compress2((Bytef*)comp_buffer, (uLongf *)&csize, (Bytef*)in_data, 
			(uLongf)data_size, Z_DEFAULT_COMPRESSION );
    if(err !=  Z_OK) {
      zlib_err(err); 
      //exit(1); 
    }
    return (acomp_size_t)csize;
}


acomp_size_t _lzo2_compression( byte_t *in_data, byte_t *comp_buffer, acomp_size_t data_size ) {
	
    int err;
    lzo_uint size = (lzo_uint)data_size;
    lzo_uint input_size = (lzo_uint)data_size;	    
    if (!init_lzo) {
      err  =  lzo_init (  );
      if (err !=  LZO_E_OK) {
	printf("Couldn't initialize lzo library\n"); 
	//exit(1);
      }
      init_lzo  =  1;
    }

    char* wrkmem  =  (char*) malloc(LZO1X_1_MEM_COMPRESS * 40);

    err  =  lzo1x_1_compress ( (lzo_bytep) in_data, (lzo_uint) input_size, 
			       (lzo_bytep) comp_buffer, 
			       (lzo_uintp) &size, (lzo_voidp) wrkmem );

    if (err !=  LZO_E_OK) { 
      lzo_err(err); 
      //exit(1);
    }   

    free(wrkmem);
    return (acomp_size_t) size;
}     


acomp_size_t _bzip2_compression( byte_t *in_data, byte_t *comp_buffer, 
			acomp_size_t data_size ) {
	
    acomp_size_t size = data_size;
    unsigned int comp_size = data_size; 
    int err  =  BZ2_bzBuffToBuffCompress( comp_buffer, &comp_size, in_data,
				 	      (unsigned int)data_size, 
					   4 /*can be 1-9 */, 0, 30 ); 

    if(err !=  BZ_OK) {
      bzip_err(err); 
      //exit(1);
    }

    size = comp_size;
    DB_PRINTF("size %lu, data size %lu, err %d\n", size, data_size, err);
    return size;
}




void _zlib_decompression( byte_t *out_data, byte_t *comp_buffer, 
			acomp_size_t comp_size, acomp_size_t data_size ) {
   
    uLong uncomp_len  =  4*data_size; 
    int err  =  uncompress((Byte*) out_data, &uncomp_len, 
					    (Byte*)comp_buffer, comp_size); 
 
    if(err !=  Z_OK) {
      zlib_err(err); //fprintf(stderr, "Zlib uncompress error: %d\n", err);
      //exit(1); 
   }
}

void _lzo2_decompression( byte_t *out_data, byte_t *comp_buffer, 
			acomp_size_t comp_size, acomp_size_t data_size ) {

    int err;
    acomp_size_t uncompress_len  =  4*data_size;

    if (!init_lzo) {
      err  =  lzo_init ( );
      if (err !=  LZO_E_OK) { 
	printf("Couldn't initialize lao library\n"); exit(1);
      }
      init_lzo  =  1;
    }
 
    void* wrkmem  =  malloc(LZO1X_1_MEM_COMPRESS);


    err  =  lzo1x_decompress ( (lzo_bytep) comp_buffer, 
			      (lzo_uint) comp_size, (lzo_bytep) out_data, 
			     (lzo_uintp) &uncompress_len, (lzo_voidp) wrkmem);
    if (err !=  LZO_E_OK) { 
      lzo_err(err); 
      //exit(1); 
   } 

   free(wrkmem);

}

void _bzip2_decompression( byte_t *out_data, byte_t *comp_buffer, 
	   	    acomp_size_t comp_size, acomp_size_t data_size ) {
 
   unsigned int uncompress_len  =  4 * data_size;
   DB_PRINTF("compressed size %lu, uncompressed size %u\n", 
   					comp_size, uncompress_len);

   int err  =  BZ2_bzBuffToBuffDecompress(out_data, &uncompress_len,
	                comp_buffer, (unsigned int)comp_size, 0, 0 );

   DB_PRINTF("compressed size %lu, uncompressed size %u\n", comp_size,
						     uncompress_len);
   
   if(err !=  BZ_OK) {
      bzip_err(err); 
      //exit(1);
    } 
  
}

acomp_size_t compress_prep (byte_t *in_data, byte_t *comp_buffer,
   		       acomp_size_t data_size, compressor_t compression_type ) {

  acomp_size_t size = data_size;

  DB_PRINTF("Compress preprocess with type %d\n", compression_type); 

  switch(compression_type) {
   case ZLIB: 
         size = _zlib_compression(in_data, comp_buffer, data_size);
         break;
   case LZO:	  	    
         size = _lzo2_compression(in_data, comp_buffer, data_size);
         break;
   case BZIP2: 
         size = _bzip2_compression(in_data, comp_buffer, data_size);
         break;
        	
   default:
  	fprintf(stderr, "Invalid compressor type\n");
  } 

  return size;
}

acomp_size_t compress_prep_by_bytecol( byte_t *in_data, byte_t *comp_buffer,
			               int n_bytecols, acomp_size_t bytecol_size,
						    compressor_t compression_type ) {
   acomp_size_t total_sz = 0, *metadata = (acomp_size_t *)comp_buffer; 
   byte_t* out_buffer = comp_buffer + n_bytecols * sizeof(acomp_size_t);	
   int comp_col = 0;

   for(comp_col = 0; comp_col < n_bytecols; comp_col ++ ) {
	metadata[comp_col] = compress_prep (in_data, out_buffer, 
   					bytecol_size, compression_type );
	DB_PRINTF("DEBUG: Size after compression[Bytecol %d] = %lu\n", comp_col,
							metadata[comp_col])		
	out_buffer += metadata[comp_col];
	total_sz += metadata[comp_col]; 
	in_data += bytecol_size;
   }
  total_sz += (sizeof(acomp_size_t) * n_bytecols);
  return total_sz;	
}



void decompress_prep (byte_t *out_data, byte_t *comp_buffer,
   		       acomp_size_t comp_size, acomp_size_t data_size,
					compressor_t compression_type ) {

  DB_PRINTF("Compress preprocess with type %d\n", compression_type); 

  switch(compression_type) {
   case ZLIB: 
          _zlib_decompression(out_data, comp_buffer, comp_size, 
						     data_size);
         break;
   case LZO:	  	    
         _lzo2_decompression(out_data, comp_buffer, comp_size,
						     data_size);
         break;
   case BZIP2: 
         _bzip2_decompression(out_data, comp_buffer, comp_size,
						    data_size);
         break;
        	
   default:
  	fprintf(stderr, "Invalid compressor type\n");
  } 

}

void decompress_prep_by_bytecol( byte_t *out_data, byte_t *comp_buffer,
                               int n_bytecols, acomp_size_t bytecol_size,
                                          compressor_t compression_type ) { 

   acomp_size_t *metadata = (acomp_size_t *)comp_buffer; 
   byte_t* tmp_buffer = comp_buffer + n_bytecols * sizeof(acomp_size_t);        
   int comp_col = 0;

   for(comp_col = 0; comp_col < n_bytecols; comp_col ++ ) { 
	/* DB_PRINTF("DEBUG: Compressed size[Bytecol %d] = %lu\n", comp_col,
							metadata[comp_col])*/		
        decompress_prep (out_data, tmp_buffer, metadata[comp_col],
				   bytecol_size, compression_type );

        tmp_buffer += metadata[comp_col];
        out_data += bytecol_size;
   }   
}



void zlib_err(int err)
{
  switch (err) {
  case Z_ERRNO:
    if (ferror(stdin))
      fprintf(stdout, "ACOMPS ERROR: zlib : Error reading stdin\n");
 
    if (ferror(stdout))
      fprintf(stdout, "ACOMPS ERROR: zlib: Error writing stdout\n");
    break;

  case Z_STREAM_ERROR:
    fprintf(stdout, "ACOMPS ERROR: zlib: Invalid compression level\n");
    break;

  case Z_MEM_ERROR:
    fprintf(stdout, "ACOMPS ERROR: zlib: Insufficient memory\n");
    break;

  case Z_BUF_ERROR:
    fprintf(stdout, "ACOMPS ERROR: zlib: The buffer dest was not large enough\n");
    break;

  case Z_DATA_ERROR:
    fprintf(stdout, "ACOMPS ERROR: zlib: Invalid or incomplete deflate data\n");
    break;

  case Z_VERSION_ERROR:
    fprintf(stdout, "ACOMPS ERROR: zlib: Version mismatch!\n"); 
    break;

  default : fprintf(stdout, "ACOMPS ERROR: zlib: unknown error!\n"); 
  }
}

void lzo_err (int err) {

  switch ( err) {
  case LZO_E_INPUT_NOT_CONSUMED: 
    fprintf(stdout, " ACOMPS LZO: 'src_len' is too large \n"); break;
  case LZO_E_INPUT_OVERRUN: 
    fprintf(stdout, " Your data is corrupted (or 'src_len' is too small)\n"); 
    break;
  case LZO_E_OUTPUT_OVERRUN: 
    fprintf(stdout, "Either your data is corrupted, or you should increase the number \
					 of available bytes passed in the variable pointed by 'dst_len'\n"); 
    break;
  case LZO_E_LOOKBEHIND_OVERRUN:
    fprintf(stdout," Your data is corrupted.\n"); break;

  case LZO_E_EOF_NOT_FOUND: 
    fprintf(stdout, "No EOF code was found in the compressed block.\
				 Your data is corrupted (or 'src_len' is too small).\n"); break;
  case LZO_E_ERROR: fprintf( stdout, "Any other error (data corrupted).\n");           
  }
}


void bzip_err (int err) {
  switch(err) {
  case BZ_CONFIG_ERROR: 
    fprintf(stdout, " ACOMPS BZIP: library has been mis-compiled\n");
    break;	
  case BZ_PARAM_ERROR: 
    fprintf(stdout, " ACOMPS BZIP : Either dest is NULL or destLen is NULL or \
		small !=  0 && small !=  1 verbosity < 0 or verbosity > 4\n");
    break;
  case BZ_MEM_ERROR: 
    fprintf(stdout, "ACOMPS BZIP: insufficient memory is available \n");
    break;
  case BZ_OUTBUFF_FULL : 
    fprintf(stdout, "ACOMPS BZIP: The size of the compressed data exceeds *destLen\n"); 
    break;
  case BZ_DATA_ERROR : 
    fprintf(stdout, "ACOMPS BZIP : data integrity error was detected in the compressed data\n"); 
    break;
  case BZ_DATA_ERROR_MAGIC : 
    fprintf(stdout, "ACOMPS BZIP: the compressed data doesn't begin with the right magic bytes\n"); 
    break;
  case BZ_UNEXPECTED_EOF : 
    fprintf(stdout, "ACOMPS BZIP: if the compressed data ends unexpectedly\n"); break;
    default : fprintf(stdout, "ACOMPS BZIP: unknown error!\n"); 
  }
}



