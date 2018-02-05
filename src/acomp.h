#ifndef ACOMP_H
#define ACOMP_H

#define DEFAULT_W_CS 0.0f
#define DEFAULT_W_CR 1.0f
#define DEFAULT_DELTA 1.00f

int IS_REAL=1;

enum { 
	NAIVE_ZLIB = 0, 
	NAIVE_LZO, 
	NAIVE_BZIP,  
	BYTE_ZLIB, 
	BYTE_LZO,
	BYTE_BZIP, 
	BYTEXOR_ZLIB, 	
	BYTEXOR_LZO, 
	BYTEXOR_BZIP, 
	ZLIB_N, 
	LZO_N, 
	BZIP_N  
};

int acomp_decompress(int t_size, int ndims, unsigned char* metadata, char* compressed_data, 
	                char* output_data, 
	                unsigned int* dimLen, unsigned int* req_start, unsigned int* req_count );

unsigned int acomp_get_max_metadata_size();

unsigned int acomp_compress( char * varname, int ndims, unsigned int* spatial_dim_len, int t_size, int is_real, 
	                unsigned int input_size, unsigned int spatial_size,  
	                char* output_buff,  char* input_buff, 
	                unsigned char * metadata);


void acomp_get_actual_size( unsigned char* metadata, unsigned int *uncompressed_size_meta );

void acomps_conf(double w_cs, double w_cr, double delta);

int acomp_get_nchunks(unsigned char* metadata, unsigned int** chunk_meta);

unsigned char acomp_get_byteflags(unsigned char* metadata);

int acomp_is_compressed(unsigned char* metadata, int* compression_method);

#endif
