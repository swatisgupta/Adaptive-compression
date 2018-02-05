/**
@author : Swati Singhal ( University of Maryland )
@contact : swati@cs.umd.edu

**/


#include <sys/time.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>

/** Currently supported compression methods **/
#include <zlib.h> 
#include <bzlib.h>
#include <lzo/lzoconf.h>
#include <lzo/lzo1x.h>
#include "acomp.h"
#include "util.h"

#define FLOAT_EQ(a,b) fabs((a) - (b)) < (DBL_EPSILON * fabs((a) + (b) ))
/*-- Global declarations ------*/
#define MAX_CHUNKS 6
#define MAX_CMP_ALGO 15
#define T_PREPROCESS 3
#define T_COMP_ALGOS 3
#define T_METHODS (T_PREPROCESS* T_COMP_ALGOS) 

char C_ALGO[T_METHODS][20] = { "NAIVE_ZLIB",
			       "NAIVE_LZO",
			       "NAIVE_BZIP",
			       "BYTE_ZLIB",
			       "BYTE_LZO",
			       "BYTE_BZIP",
			       "BYTEXOR_ZLIB",
			       "BYTEXOR_LZO",
			       "BYTEXOR_BZIP",
};

enum {
  BSEG =  0,
  BWSEG,
  BWXOR,
};

struct chunk_info_t 
{
  int  chunk_id;
  char     *data;
  unsigned int  chunk_size;
  unsigned int  chunk_offset; 
  unsigned int  compSize; 

};

double weight_cr =  DEFAULT_W_CR;
double weight_cs =  DEFAULT_W_CS;
double delta =  DEFAULT_DELTA;

int init_lzo = 0;

int CHUNK_UNIT = 2097152;
char *varname = NULL;


unsigned int compress_prep ( char *in_data, char *comp_buffer, unsigned int data_size, int compression_type );

unsigned char get_byte_flags(int dtype, int t_size, unsigned int data_size, 
				   unsigned char* data, int *ncompbytes) 
{
  unsigned char byte_flags = 0;
  int *byte_freq =  NULL;
  int *max_freq = NULL;
  int i, j, measure = 1;
  unsigned int data_nbytes = data_size/t_size;
  unsigned char tmp; 
  double skewness = 0, std_devitaion = 0, kurtosis = 0, S_thres = 0, K_thres = 0;
  double entropy[UCHAR_MAX], t_entropy = 0.0f, E_thres = CHAR_BIT - 1, val;
  printf("Inside get_byte_flags_dist \n");fflush(stdout);	
  if ( dtype !=  IS_REAL ) {
    *ncompbytes = t_size;
    return UCHAR_MAX - 1;	
  }

  byte_freq = calloc(t_size*UCHAR_MAX, sizeof(int));
  max_freq = calloc(t_size, sizeof(int));

  if(byte_freq ==  NULL || max_freq == NULL)  {
    fprintf(stdout, "unable to allocate memory for computing compressible bytes");
    return UCHAR_MAX - 1;
  }
	
  *ncompbytes = t_size;

  for( i = 0; i < t_size; i++) {

    for( j = 0; j< data_nbytes; j++) {
      tmp = data[j * t_size  + i];
      byte_freq[ i * UCHAR_MAX + (unsigned int)tmp]++;
      //printf("tmp = %u, i * UCHAR_MAX + (unsigned int)tmp - %u\n", tmp, i * UCHAR_MAX + (unsigned int)tmp ); 
      if(max_freq[i] < byte_freq[ i * UCHAR_MAX + (unsigned int)tmp]) 
	max_freq[i] = byte_freq[ i * UCHAR_MAX + (unsigned int)tmp];
    }
    
    double _mean = 0;  
    for ( j = 0; j < UCHAR_MAX; j ++ ) {
	_mean +=  j * byte_freq[ i * UCHAR_MAX + j ];
    }	 
    _mean /=  UCHAR_MAX;
                
    switch (measure) {
    case 1: /* Entropy */

      t_entropy  =  0;
      for ( j =  0; j < UCHAR_MAX; j++) {
	val  =  (double)byte_freq[ i * UCHAR_MAX + j]/(double)data_nbytes;

	if( FLOAT_EQ(val, 0) )
		 entropy[j]  =  0;
	else 
	  entropy[j]  =  -1 * val * log2(val);
	  //entropy[j]  =  -1 * log2(val);
        
	//printf("Entropy of byte column %i: state %d (%d), (-1* val* log2(val)  =  -1 * %g * %g)  =  %g \n", i, j, byte_freq[ i * UCHAR_MAX + j],  val, log2(val), entropy[j]);
	//printf("Entropy of byte column %i: state %d (%d), (-1* log2(%g)  =  -1 * %g)  =  %g \n", i, j, byte_freq[ i * UCHAR_MAX + j], val, log2(val), entropy[j]);                           
	t_entropy +=  entropy[j];                                                          
      }
    
      if(t_entropy > E_thres || i  ==  1 )   {    /* required more than 6 bits to encode a character
					*/
          printf("Entropy of the byte column %d is %g high\n", i, t_entropy);					
          (*ncompbytes) --;
      } else {
	 byte_flags |=  1<<i;
      }
   
      printf("Entropy of the byte column %d is %g\n", i, t_entropy);					
      break;

    case 2: /* Skewness */
      skewness  =  0;
      std_devitaion  =  0;
      for( j  =  0; j < UCHAR_MAX; j ++) { 	
      	skewness +=  pow(byte_freq[ i * UCHAR_MAX + j] - _mean, 3)/UCHAR_MAX;
      	std_devitaion +=  pow(byte_freq[ i * UCHAR_MAX + j] - _mean, 2);
      }
      std_devitaion  =  sqrt(std_devitaion); 
      skewness /=  pow(std_devitaion, 3); 	
      printf(" Skewness of byte column %d is %g\n", i, skewness );
      
      if ( skewness > S_thres )   {    
          printf("Skewness of the byte column %d is %g high\n", i, skewness);					
          (*ncompbytes) --;
      } else {
	 byte_flags |=  1<<i;
      }
 
      break;

    case 3: /* Kurtosis */
      kurtosis  =  0;
      std_devitaion  =  0;
      for( j  =  0; j < UCHAR_MAX; j ++) { 	
      	kurtosis +=  pow(byte_freq[ i * UCHAR_MAX + j] - _mean, 4)/UCHAR_MAX;
      	std_devitaion +=  pow(byte_freq[ i * UCHAR_MAX + j] - _mean, 2);
      }
      std_devitaion  =  sqrt(std_devitaion); 
      kurtosis /=  pow(std_devitaion, 4); 	
      // kurtosis /=  pow(std_devitaion, 4) - 3; 	
      printf(" Kurtosis of byte column %d is %g\n", i, kurtosis );

      if(kurtosis > K_thres )   {  
          printf("kurtosis of the byte column %d is %g high\n", i, kurtosis);					
          (*ncompbytes) --;
      } else {
	 byte_flags |=  1<<i;
      }

      break;

    default: 
      if(max_freq[i] < 0.58 * data_nbytes)   {    
	(*ncompbytes) --;
      } else {
	byte_flags |=  1<<i;
      }
    } 
  }
	
  if( byte_flags  ==  0) {
    int sz  =  t_size/2; 
    for( i  =  0 ; i < sz; i ++ ) {
      byte_flags |=  1<<(t_size - i - 1);
    }
    *ncompbytes  =  sz; 
  }
  /**
     for(i = 0; i < t_size; i++) {
     printf("Byte %d, max freq %d\n", i, max_freq[i]);
     }
  **/
  //	printf("Byte FLags %u, Datatype size: %d, number of compressible bytes: %d\n", byte_flags, t_size , *ncompbytes);
  free(byte_freq);
  free(max_freq);

  return byte_flags;
}	

int get_chunks_1d (struct chunk_info_t **dchunks, unsigned int* diml, int ndims, int t_size) {

  int i  =  ndims -1, nchunks, total_size  =  1;
  unsigned int chunk_sz  =  CHUNK_UNIT/ t_size;
  struct chunk_info_t * chunks  =  NULL;	

  for(; i >=  0 ; i--) {
    total_size *=  diml[i];
  }
	
  nchunks  =  total_size / chunk_sz;
       
  if(total_size % chunk_sz) 
    nchunks ++;
 
  if( nchunks > MAX_CHUNKS) { 
    chunk_sz  =  total_size / MAX_CHUNKS;
    nchunks  =   total_size / chunk_sz;
  }

  if( nchunks  ==  0) 
    nchunks  =  1;
	
  *dchunks  =  (struct chunk_info_t *) malloc(sizeof(struct chunk_info_t) * nchunks);

  assert(*dchunks !=  NULL );
	
  chunks  =  *dchunks;

  chunk_sz *=  t_size;
	
  printf("chunk_size %u, total_size %u, nchunks %d ndims %d\n", chunk_sz, total_size * t_size, nchunks, ndims);
   
  for (i  =  0; i< nchunks ; i++) {
    chunks[i].chunk_offset  =  i * chunk_sz;
    chunks[i].chunk_id  =  i;
    if( i  ==  nchunks - 1) {
      chunks[i].chunk_size  =  (total_size * t_size) - i * chunk_sz;
      chunks[i].compSize  =  (total_size * t_size) - i * chunk_sz;
    } else  {
      chunks[i].chunk_size  =  chunk_sz;
      chunks[i].compSize  =  chunk_sz;
    }
    chunks[i].data  =  NULL;
    printf("chunk %d  => offset %u, size %u\n",chunks[i].chunk_id  =  i ,chunks[i].chunk_offset, chunks[i].chunk_size);
  }

  return nchunks;
}


int get_overlap_chunks_1d( struct chunk_info_t **dchunks, int ndims, unsigned int* dimLen, int nchunks, 
			   unsigned int* subreq_start, unsigned int* subreq_count, int t_size) {

  int i =  ndims - 1, treq_size  =  t_size, total_size  =  1, offset  = 0;
  int chunk_start, chunk_end;
  unsigned int chunk_sz; 
  struct chunk_info_t * chunks  =  NULL;

  for(; i>= 0; i--) { 
    offset *=  (subreq_start[i]* t_size);
    treq_size *=  subreq_count[i];   
    total_size *=  dimLen[i];   
  }
	
  chunk_sz  =  (total_size * t_size) / nchunks; 

  if( total_size * t_size != chunk_sz * nchunks) {
    chunk_sz = (CHUNK_UNIT / t_size) * t_size;
  } 

  chunk_start  =  offset / chunk_sz;

  chunk_end  =  (offset + treq_size) / chunk_sz; 
	
  if( chunk_end >=  nchunks) 
    chunk_end  =  nchunks -1;

  *dchunks  =  malloc(sizeof(struct chunk_info_t)*(chunk_end - chunk_start + 1));

  assert(*dchunks !=  NULL);
   
  chunks  =  *dchunks;
  printf("Noverlap %d, nchunks %d,  chunk_start %d, chunk_end %d, treq_size %u, total_size %u\n",
	 (chunk_end - chunk_start + 1), nchunks, chunk_start, chunk_end, treq_size, total_size * t_size);
        

  for( i  =  0 ; chunk_start <=  chunk_end; chunk_start++) {
    chunks[i].chunk_id  =  chunk_start;
    chunks[i].chunk_offset  =  chunk_start * chunk_sz;
    if(chunk_start  ==  chunk_end) { 
      chunks[i].chunk_size  =  (total_size * t_size) - chunk_start * chunk_sz;
      chunks[i].compSize  =  chunks[i].chunk_size;
    } else {
      chunks[i].chunk_size  =  chunk_sz;
      chunks[i].compSize  =  chunk_sz;
    }
    chunks[i].data  =  NULL; 
    printf("chunk %d  => chunk_offset %u, chunk_size %u\n", chunks[i].chunk_id, chunks[i].chunk_offset, chunks[i].chunk_size);
    i++;
  }
	
  return i; 
}

#if 0
int recreate_chunk_1d(unsigned char* input_data, unsigned char* output_data, unsigned int offset, 
		      unsigned int count, int t_size) {
    
  unsigned char* out_data  =  output_data;
  unsigned char* in_data  =  input_data;
  int datalen  =  0; 

  memcpy(out_data + offset * t_size, in_data, count * t_size);
  datalen  =  count * t_size;
  return datalen;

}


int extract_chunk_data_1d(unsigned char* input_data, unsigned char* output_data, unsigned int offset, 
			  unsigned int count, int ndims, int t_size) {
 
  unsigned char* out_data  =  output_data;
  unsigned char* in_data  =  input_data;
  int datalen  =  0; 
    
  memcpy(out_data, in_data + offset * t_size, count * t_size);
  datalen  =  count * t_size;
  return datalen;
}
#endif


/***  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  PREPROCESSING METHODS  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  =  ***/
void preprocess_bSeg(unsigned char byteArray, int t_size, unsigned int data_size, char* comp_buffer, 
		     char* uncomp_buffer, char* in_buffer) {

  int startOffset  =  0, byteCount  =  0, ubyteCount  =  0;
  unsigned int varid  =  0, byte;

  for ( varid  =  0; varid < data_size ; ++varid ) {
    unsigned int idx  =  startOffset + varid * t_size;
    for( byte  =  0; byte < t_size; ++byte) {
      if (( byteArray >> byte) & 1) {
	comp_buffer[byteCount]  =  in_buffer[idx + byte]; 
	++byteCount;
      }   else {
	uncomp_buffer[ubyteCount]  =  in_buffer[idx + byte]; 
	ubyteCount++;
      }
    }
  }
}

void preprocess_bWiseSegXor(unsigned char byteArray, int t_size, unsigned int data_size, char* comp_buffer, 
			    char* uncomp_buffer, char* in_buffer) {

  unsigned int varid  =  0, byte_pos[8], byte  =  0; 
  int nextOffset  =  0;
  unsigned int ncb_count  =  0, bcount  =  0;
       
  for ( byte  =  0; byte < t_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount;
      bcount ++;	 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    printf("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
  }
  printf("Non compressible recount %d\n", ncb_count);	

  for( varid  =  0; varid < data_size; varid ++) {
    for( byte  =  0; byte < t_size; ++ byte) {
      unsigned int curr_index  =  nextOffset + varid * t_size + byte;
      unsigned int next_index  =  nextOffset + (varid + 1) * t_size + byte;
      unsigned int bindex  =  byte_pos[byte] * data_size + varid;
      if(( byteArray >> byte) & 1) {
	if( varid < data_size - 1 ) 
	  comp_buffer[bindex]  =   in_buffer[next_index] ^ in_buffer[curr_index];
	else 
	  comp_buffer[bindex]  =   in_buffer[curr_index];
      } else {
	uncomp_buffer[bindex]  =  in_buffer[curr_index];
      }
      // *data_ptr >>=  8; 	
    }
    //data_ptr ++;	
  }
}

void preprocess_bWiseSeg(unsigned char byteArray, int t_size, unsigned int data_size, char* comp_buffer, 
			 char* uncomp_buffer, char* in_buffer) {
  unsigned int varid  =  0, byte_pos[8], byte  =  0; 
  unsigned int ncb_count  =  0, bcount  =  0;
     
  int nextOffset  =  0;
  nextOffset  =  0;

    

  for ( byte  =  0; byte < t_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount; ++bcount; 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    printf("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
  }

  for( varid  =  0; varid < data_size; varid ++) {
    for( byte  =  0; byte < t_size; ++ byte) {
      unsigned int curr_index  =  nextOffset + varid * t_size + byte;
      unsigned int bindex  =  byte_pos[byte] * data_size + varid;
      if(( byteArray >> byte) & 1) { 
	comp_buffer[bindex]  =   in_buffer[curr_index];
      } else {
	uncomp_buffer[bindex]  =  in_buffer[curr_index];
      }
    }
  }

}

void postprocess_bSeg(unsigned char byteArray, int t_size, unsigned int data_size, char* comp_buffer, 
		      char* uncomp_buffer, char* in_buffer) {
  int startOffset  =  0, byteCount  =  0, ubyteCount  =  0;
  unsigned int varid  =  0, byte  =  0; 
	
  for ( varid  =  0; varid < data_size; ++varid) {
    unsigned int idx  =  startOffset + varid * t_size;
    for( byte  =  0; byte < t_size; ++byte) {
      if (( byteArray >> byte) & 1) {
	in_buffer[idx + byte]  =  comp_buffer[byteCount];
	++byteCount;
      } else {
	in_buffer[idx + byte] =  uncomp_buffer[ubyteCount];
	ubyteCount++;
      }
    }
  }
}


void postprocess_bWiseSeg(unsigned char byteArray, int t_size, unsigned int data_size, char* comp_buffer, 
			  char* uncomp_buffer, char* in_buffer) {

  unsigned int varid  =  0, byte_pos[8], byte  =  0; 
  int nextOffset  =  0; 
  unsigned int ncb_count  =  0, bcount  =  0;


  for ( byte  =  0; byte < t_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount; ++bcount; 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    printf("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
  }  


	
  for( varid  =  0; varid < data_size; varid ++) {
    for( byte  =  0; byte < t_size; ++ byte) {
      unsigned int curr_index  =  nextOffset + varid * t_size + byte;
      unsigned int bindex  =  byte_pos[byte] * data_size + varid;
      if(( byteArray >> byte) & 1) {
	in_buffer[curr_index]  =  comp_buffer[bindex];
      } else {
	in_buffer[curr_index]  =  uncomp_buffer[bindex];
      }
    }
  }
}

void postprocess_bWiseSegXor(unsigned char byteArray, int t_size, unsigned int data_size, char* comp_buffer, 
			     char* uncomp_buffer, char* in_buffer) {
  unsigned int varid  =  0, byte_pos[8], byte  =  0; 
  unsigned int ncb_count  =  0, bcount  = 0;
  int nextOffset  =  0; 

  for ( byte  =  0; byte < t_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount; ++bcount; 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    printf("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
  }

  for( varid  =  data_size - 1; varid >=  0; varid --) {
     
    for( byte  =  0; byte < t_size; ++ byte) {
      unsigned int curr_index  =  nextOffset + varid * t_size + byte;
      unsigned int next_index  =  nextOffset + (varid + 1) * t_size + byte;
      unsigned int bindex  =  byte_pos[byte] * data_size + varid;
      if(( byteArray >> byte) & 1) {
	//printf("varid %lu bindex %lu\n", varid, bindex);
	if (varid < data_size -1 ) 
	  in_buffer[curr_index]  =  comp_buffer[bindex] ^ in_buffer[next_index];
	else 
	  in_buffer[curr_index]  =  comp_buffer[bindex];	
      } else {
	//printf("Uncompress bindex %lu\n", bindex);
	in_buffer[curr_index]  =  uncomp_buffer[bindex];
      }
    }
 		
    if (varid  ==  0) 
      break;
  }

}
/*****  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ****/


void process_data (unsigned char byteArray, char *in_buffer, char *comp_buffer, char *uncomp_buffer, 
		   int t_size, unsigned int data_size, int transform_type ) {

  switch ( transform_type ) {
  case BSEG:    /**---collecting compressible bytes ---***/
    preprocess_bSeg(byteArray, t_size, data_size, comp_buffer, 
		    uncomp_buffer, in_buffer);	
    break;

  case BWSEG: /** -- collecting compressible bytes seperately -- **/
    preprocess_bWiseSeg(byteArray, t_size, data_size, comp_buffer, 
			uncomp_buffer, in_buffer);
			    
    break;        

  case BWXOR:  /** --- collecting compressible bytes seperately + Xor ---**/
    preprocess_bWiseSegXor(byteArray, t_size, data_size, comp_buffer, 
			   uncomp_buffer, in_buffer);
    break;

  case (T_PREPROCESS + BSEG):    /**--- uncollecting compressible data ---***/
    postprocess_bSeg(byteArray, t_size, data_size, comp_buffer, 
		     uncomp_buffer, in_buffer);	
    break;
        
  case T_PREPROCESS + BWSEG: /**-----  uncollecting compressible bytes seperately ---***/
    postprocess_bWiseSeg(byteArray, t_size, data_size, comp_buffer, 
			 uncomp_buffer, in_buffer);
                        
    break;
  case T_PREPROCESS + BWXOR:  /** --- uncollecting compressible bytes seperately + Xor ---**/
    printf("Inside BWXOR post : size %u\n", data_size);
    postprocess_bWiseSegXor(byteArray, t_size, data_size, comp_buffer, 
			    uncomp_buffer, in_buffer);
    break;

  default:
    printf("Error: unsupported preprocess method\n");        
    //printf("Here...data_size %u t_size %d \n", data_size, t_size);
                         

  }     
}
 
int analyzeData(char * varname, char *input_data, unsigned int total_bytes, int dtype, int t_size, 
		unsigned int data_size, double weight_cr, double weight_ct, int* compress_algo, double* cPerform, unsigned char* byteFlags, int *n_comp_bytes) {

  int i, j, index, idx;
     
  char *comp_buffer, *uncomp_buffer;
  char *prep_data;
  double analysis_data[MAX_CMP_ALGO];  

  int bAcnt  =  0;

  index  =  find_index(varname, &idx);

  if( index  ==  -1) {
    index  =  idx; 
    assert(index !=  -1); 
    strcpy(hash_table[index].varname, varname);
    hash_table[index].is_empty  =  1;
  }

  printf("Variable : %s found at index %d \n", varname, index); fflush(stdout);
  *compress_algo  =  hash_table[index].compress_algo;
  
  unsigned char byteArray  =  get_byte_flags(dtype , t_size,  total_bytes, (unsigned char*) input_data, &bAcnt);
  *n_comp_bytes  =  bAcnt;
  *byteFlags  =  byteArray;

  if( ( *compress_algo !=  -1 ) && 
	 ~( hash_table[index].comp_perform < hash_table[index].avg_comp_perform + delta &&
	    hash_table[index].comp_perform > hash_table[index].avg_comp_perform - delta ) ) {
    return 0;
  }    

    
  //printf("analyze data: input_data %p, total_bytes %llu, data_size %llu type_size %llu\n", input_data, total_bytes, data_size, t_size);
  //printf("analyze data: byteArray %d\n", byteArray);
  //printf("weight cr %g, weight cs %g\n", weight_cr, weight_ct);
  int cSize;

  //printf("analyze data : bAcnt %d\n", bAcnt); 
  unsigned int size  =  data_size*t_size >1024? 4*data_size*t_size: 4*1024;
  comp_buffer  =  malloc(size * sizeof(char));
  uncomp_buffer  =  malloc(size * sizeof(char));
  prep_data  =  malloc(size * sizeof(char));
  //unsigned char* output_data  =  malloc(size * sizeof(unsigned char));

  if(comp_buffer  ==  NULL || uncomp_buffer  ==  NULL || prep_data  ==  NULL) { 
    fprintf (stderr, "ADIOS ERROR: failed while allocating  \
				 comp_buffer and uncomp_buffer in analyzeData() \n");
    exit(-1);
  }

  struct timeval start_time, end_time;
  double cTime, cSpeed, cRatio, res;  

  memcpy(prep_data, input_data, data_size);
  *cPerform  =  -0.0f;
  *compress_algo  =  -1;

  for ( i  =  0; i < T_PREPROCESS; i++ ) {
    process_data(byteArray, input_data, prep_data, uncomp_buffer, t_size, data_size, i);
    for ( j  =  0; j < T_COMP_ALGOS; j++) {
      gettimeofday(&start_time, NULL);    
      cSize  =  compress_prep(prep_data, comp_buffer, data_size * t_size, j);
      gettimeofday(&end_time, NULL); 
      cTime  =  (double)(end_time.tv_usec - start_time.tv_usec) / 1000000 + (double)(end_time.tv_sec - start_time.tv_sec);
      cSpeed  =  (data_size*t_size)/cTime;	
      cSpeed  =  cSpeed/(1024*1024); 
      cRatio  =  (double)(data_size*bAcnt)/cSize;
      res  =  cRatio * weight_cr + cSpeed * weight_ct;
      analysis_data[i*T_COMP_ALGOS + j]  =  res; 
            
      if( *cPerform  < res )  { 
	*cPerform  =  analysis_data[ i * T_COMP_ALGOS + j ];
	*compress_algo  =  i * T_COMP_ALGOS + j;
      }
      printf("%s, %s[%d]: compress_ratio : %g, compress_speed %g byte/second, compress_size %d perform %g\n", varname, C_ALGO[ i * T_COMP_ALGOS+ j ], (i * T_COMP_ALGOS+ j), cRatio, cSpeed, cSize, analysis_data[ i * T_COMP_ALGOS+ j ]);
    }

  }
  //*compress_algo  =  6;
  printf("%s: compression algo %d, compression perform %g\n", varname, *compress_algo, *cPerform);      
  return 1;
}


static void zerr(int ret)
{
  switch (ret) {
  case Z_ERRNO:
    if (ferror(stdin))
      fprintf(stdout, "ADIOS ERROR: zlib : Error reading stdin\n");
 
    if (ferror(stdout))
      fprintf(stdout, "ADIOS ERROR: zlib: Error writing stdout\n");
    break;

  case Z_STREAM_ERROR:
    fprintf(stdout, "ADIOS ERROR: zlib: Invalid compression level\n");
    break;

  case Z_MEM_ERROR:
    fprintf(stdout, "ADIOS ERROR: zlib: Insufficient memory\n");
    break;

  case Z_BUF_ERROR:
    fprintf(stdout, "ADIOS ERROR: zlib: The buffer dest was not large enough\n");
    break;

  case Z_DATA_ERROR:
    fprintf(stdout, "ADIOS ERROR: zlib: Invalid or incomplete deflate data\n");
    break;

  case Z_VERSION_ERROR:
    fprintf(stdout, "ADIOS ERROR: zlib: Version mismatch!\n"); 
    break;

  default : fprintf(stdout, "ADIOS ERROR: zlib: unknown error!\n"); 
  }
}

static void lzo_err (int ret) {

  switch ( ret) {
  case LZO_E_INPUT_NOT_CONSUMED: 
    fprintf(stdout, " ADIOS LZO: 'src_len' is too large \n"); break;
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


static void bzip_err (int ret) {
  switch(ret) {
  case BZ_CONFIG_ERROR: 
    fprintf(stdout, " ADIOS BZIP: library has been mis-compiled\n");
    break;	
  case BZ_PARAM_ERROR: 
    fprintf(stdout, " ADIOS BZIP : Either dest is NULL or destLen is NULL or \
						small !=  0 && small !=  1 verbosity < 0 or verbosity > 4\n");
    break;
  case BZ_MEM_ERROR: 
    fprintf(stdout, "ADIOS BZIP: insufficient memory is available \n");
    break;
  case BZ_OUTBUFF_FULL : 
    fprintf(stdout, "ADIOS BZIP: The size of the compressed data exceeds *destLen\n"); 
    break;
  case BZ_DATA_ERROR : 
    fprintf(stdout, "ADIOS BZIP : data integrity error was detected in the compressed data\n"); 
    break;
  case BZ_DATA_ERROR_MAGIC : 
    fprintf(stdout, "ADIOS BZIP: the compressed data doesn't begin with the right magic bytes\n"); 
    break;
  case BZ_UNEXPECTED_EOF : 
    fprintf(stdout, "ADIOS BZIP: if the compressed data ends unexpectedly\n"); break;
  default : fprintf(stdout, "ADIOS BZIP: unknown error!\n"); 
  }
}

void decompress_prep(char *out_data, char *comp_buffer, unsigned int comp_size, unsigned int data_size, int decompression_type) {
     
  unsigned int uncompress_len  =  4*data_size;
  if ( decompression_type  ==  0) {
    uLong uncomp_len  =  4*data_size; 
    int err  =  uncompress((Byte*) out_data, &uncomp_len, (Byte*)comp_buffer, comp_size); 
 
    if(err !=  Z_OK) {
      zerr(err); //fprintf(stderr, "Zlib uncompress error: %d\n", err);
      exit(1);
    }
  } else if ( decompression_type  ==  1 ) {
    int ret;

    if (!init_lzo) {
      ret  =  lzo_init ( );
      if (ret !=  LZO_E_OK) { 
	printf("Couldn't initialize lao library\n"); exit(1);
      }
      init_lzo  =  1;
    }
 
    void* wrkmem  =  malloc(LZO1X_1_MEM_COMPRESS);
    //printf("decompress: %p \n", wrkmem); //fflush(stdout);
    ret  =  lzo1x_decompress ( (lzo_bytep) comp_buffer, (lzo_uint) comp_size, (lzo_bytep) out_data, 
			     (lzo_uintp) &uncompress_len, (lzo_voidp) wrkmem);
    if (ret !=  LZO_E_OK) { 
      lzo_err(ret); exit(1);
    } 
    //printf("decompress: %p \n", wrkmem); //fflush(stdout);
    free(wrkmem);

  } else {
    printf("compressed size %u, uncompressed size %u\n", comp_size, uncompress_len);
    int err  =  BZ2_bzBuffToBuffDecompress (out_data, &uncompress_len, comp_buffer, comp_size, 0, 0 );
    printf("compressed size %u, uncompressed size %u\n", comp_size, uncompress_len);
    if(err !=  BZ_OK) {
      bzip_err(err); //fprintf(stderr, "Zlib uncompress error: %d\n", err);
      exit(1);
    } 
  }
 
}


unsigned int compress_prep (char *in_data, char *comp_buffer, unsigned int data_size, int compression_type ) {

  unsigned int size  =  data_size > 4 * 1024 ? 4 * data_size : 4 * 1024; 

  //printf("Compress preprocess with type %d\n", compression_type); 
  if( compression_type  ==  0 ) {    
    uLongf csize  =  data_size > 4* 1024? 4*data_size:4* 1024;
    int err  =  compress2((Bytef*)comp_buffer, (uLongf *)&csize, (Bytef*)in_data, 
			(uLongf)data_size, Z_DEFAULT_COMPRESSION );
    size  =  csize;
    if(err !=  Z_OK) {
      zerr(err); exit(1); //fprintf(stderr, "Zlib compress error: %d\n", err); exit(1);
    }
   
  } else if (compression_type  ==  1) {

    int ret;
       
    if (!init_lzo) {
      ret  =  lzo_init (  );
      if (ret !=  LZO_E_OK) {
	printf("Couldn't initialize lao library\n"); exit(1);
      }
      init_lzo  =  1;
    }

    char* wrkmem  =  (char*) malloc(LZO1X_1_MEM_COMPRESS);
    ret  =  lzo1x_1_compress ( (lzo_bytep) in_data, (lzo_uint) data_size, (lzo_bytep) comp_buffer, 
			     (lzo_uintp) &size, (lzo_voidp) wrkmem);

    if (ret !=  LZO_E_OK) { 
      lzo_err(ret); exit(1);
    }   

    free(wrkmem);
     
  } else {
    int err  =  BZ2_bzBuffToBuffCompress( comp_buffer, &size, in_data, data_size, 4 /*can be 1-9 */, 0, 30 ); 
    if(err !=  BZ_OK) {
      bzip_err(err); exit(1); //fprintf(stderr, "Zlib compress error: %d\n", err); exit(1);
    }
    printf("size %u, data size %u, err %u\n", size, data_size, err);
    assert(size < data_size);	
  }
  return size;
}


void verify(unsigned char byteArray, unsigned char** bSeg, unsigned char** dbSeg, int t_size, unsigned int data_size) {
  /* verify */

  unsigned char* dbSegment, *bSegment;

  int byte, varIdx;
   
  for ( varIdx  =  0; varIdx < data_size; ++varIdx) {
    int idx  =  varIdx;
    for ( byte = 0; byte < t_size; ++byte) {
      dbSegment  =  dbSeg[byte];
      bSegment  =  bSeg[byte];

      if ( ( byteArray >> byte) & 1) {
	if (bSegment[idx] !=  dbSegment[idx]) {
	  fprintf(stderr, "Mismatch byte : %d varIdx: %d\n", byte, varIdx);
	  exit(1);
	}
      }
    }
  }
   
}


void verify_compression(unsigned char byteArray, unsigned char* input_data, 
			unsigned char* output_data, int t_size, 
			unsigned int data_size) {
  /* verify */

  unsigned char* dbSegment, *bSegment;

  int byte, varIdx;
  dbSegment  =  output_data;
  bSegment  =  input_data;

  for ( varIdx  =  0; varIdx < data_size; ++varIdx) {
    int idx  =  varIdx * t_size;
    for ( byte = 0; byte < t_size; ++byte) {
      if ( ( byteArray >> byte) & 1) {
	if (bSegment[idx + byte] !=  dbSegment[idx + byte]) {
	  fprintf(stderr, "Mismatch at byte: %d varIdx: %d\n", byte, varIdx);
	  exit(1);
	}
      }
    }
  }
   
}

void acomps_conf(double w_cs, double w_cr, double del) {
  weight_cr  =  w_cr;
  weight_cs  =  w_cs;
  //printf("weight cr %f \n", weight_cr); fflush(stdout);
  //printf("weight cs %f \n", weight_cs); fflush(stdout);
  delta  =  del;

}

unsigned int acomp_compress( char * varname, int ndims, unsigned int * spatial_dim_len, int t_size, int dtype, 
			     unsigned int input_size, unsigned int spatial_size,  
			     char* output_buff,  char* input_buff, 
			     unsigned char * metadata) {

  int chunk_id  =  0, nchunk  =  0, idx, index;
  struct chunk_info_t* chunks;
  unsigned char byteFlags;
  int  n_comp_bytes;
  double compress_perform  =  0;
  // unsigned int test_size  =  input_size > 1024? 1024:input_size; 
  int update_log  =  0, compress_algo  =  -1;
  char *comp_buff  =  NULL, *uncomp_buff =  NULL;
  unsigned char compress_ok  =  1;
  unsigned int actual_output_size  =  0; //output_size;


  if(read_analysis_log()  ==  -1) 
    initialize_hash();

  printf("Inside heap\n"); fflush(stdout);

  update_log  =  analyzeData(varname, input_buff, input_size, dtype, t_size, 
			   spatial_size, weight_cr, weight_cs, &compress_algo, 
			   &compress_perform, &byteFlags, &n_comp_bytes);
        
  //byteFlags  =  get_byte_flags(dtype, t_size, input_size,  input_buff, &n_comp_bytes); 

  printf("N compressible bytes %d of %d bytes",n_comp_bytes, t_size);
 
  comp_buff  =  output_buff + (t_size - n_comp_bytes) * spatial_size; 
      
  uncomp_buff  =  output_buff; 
  nchunk  =  get_chunks_1d(&chunks, spatial_dim_len, ndims, t_size);

  struct timeval start_time, end_time;
  double cTime, total_comp_time = 0, compression_speed  =  0, compression_ratio;
     
  while ( chunk_id < nchunk && compress_ok) {
 
    gettimeofday(&start_time, NULL); 
    chunks[chunk_id].data  =  input_buff + chunks[chunk_id].chunk_offset; 
    printf("Processing chunk %d[of %d] offset %u, buffer pointer %p\n", chunk_id, nchunk, chunks[chunk_id].chunk_offset, chunks[chunk_id].data);
           

    unsigned int d_size  =  chunks[chunk_id].chunk_size / t_size; 
    unsigned int compressed_size;
    unsigned int num_comps  =  n_comp_bytes;

    char* prep_data  =  (char*) malloc(d_size * num_comps* 10);     
      
    if(compress_algo <=  T_METHODS )
      process_data(byteFlags, chunks[chunk_id].data, prep_data, uncomp_buff,
		   t_size, d_size, compress_algo % T_PREPROCESS);
    else {
      memcpy(prep_data, chunks[chunk_id].data, d_size*t_size);
      num_comps  =  t_size;
    } 	   
    printf("Data size %u total_dsize %u\n", d_size, d_size * num_comps);
 
    compressed_size  =  compress_prep ( prep_data, comp_buff, d_size * num_comps, 
				      compress_algo % T_PREPROCESS);
  
    free(prep_data);
    if( compressed_size >=  d_size * num_comps) {
      //printf("compressed_size %u input_size %d \n", compressed_size, num_comps );
      compressed_size  =  -1;
      compress_ok  =  0;
 
    }
   
    chunks[chunk_id].compSize  =  compressed_size;

    gettimeofday(&end_time, NULL); 

    cTime  =  (double)(end_time.tv_usec - start_time.tv_usec) / 1000000 
      + (double)(end_time.tv_sec - start_time.tv_sec);

    total_comp_time +=  cTime;  
    printf("Done compression and compression size  =  %u  \n", compressed_size);
    actual_output_size +=  compressed_size;
    comp_buff +=  compressed_size;
    uncomp_buff +=  (chunks[chunk_id].chunk_size/t_size) *(t_size - n_comp_bytes);
    chunk_id ++;
  }

  compression_ratio  =  ((double)(spatial_size * n_comp_bytes))/actual_output_size;
  compression_speed  =  ((double) (spatial_size * n_comp_bytes))/(1024*1024);
  compression_speed  =  compression_speed/total_comp_time;

  compress_perform  =  weight_cr * compression_ratio + weight_cs * compression_speed;
  //printf("done compressing. Compressed size  =  %d\n", actual_output_size);

  if(compress_ok !=  1) {
    actual_output_size  =  input_size;
  }
  else {
    actual_output_size +=  (t_size - n_comp_bytes) * spatial_size;
  }
   
  compress_ok  =  compress_algo << 1  | compress_ok;

  // copy the metadata, simply the original size before compression
  unsigned char *metadata1  =  metadata; 
  if( metadata )
    {   int i;
      unsigned char nc  =  nchunk;
      unsigned int total_compress  =  0;	

      memcpy((unsigned char*)metadata, &input_size, sizeof(unsigned int));
      metadata +=  sizeof(unsigned int);

      memcpy((unsigned char*)metadata, &compress_ok, sizeof(unsigned char));
      metadata +=  sizeof(unsigned char);

      memcpy((unsigned char*)metadata, &byteFlags, sizeof(unsigned char));
      metadata +=  sizeof(unsigned char);

      memcpy((unsigned char*)metadata, &nc, sizeof(unsigned char));
      metadata +=  sizeof(unsigned char);

      for(i  =  0; i< nchunk; i++) {
	printf("chunk id %d, chunk offset %u, chunk size %u\n", i, total_compress, chunks[i].compSize);
	total_compress +=  chunks[i].compSize;
	memcpy((unsigned char*)metadata, &total_compress, sizeof(unsigned int));
	printf("written compression size %u\n", *(unsigned int*)metadata);
	metadata +=  sizeof(unsigned int);
      }
    }

    
  int j  =  0;
  char* uncompressed_data  =  malloc(4*input_size);
  unsigned int dimLen[5], req_start[5], req_count[5];
  for(j = 0; j<5; j++) {
    dimLen[j]  =  1;
    req_start[j]  =  0;
    req_count[j]  =  1;
  }			
  dimLen[ndims-1]  =  input_size/t_size; // + input_size%t_size;
  req_count[ndims-1]  =  input_size/t_size; // + input_size%t_size;
  printf("Ndims %d, %u, DimLen[%d] = %u, req_count[%d] = %u\n", ndims, input_size/t_size ,ndims-1, dimLen[ndims-1], ndims -1, req_count[ndims-1]);

  acomp_decompress(t_size, ndims, metadata1, output_buff,
		   uncompressed_data,
		   dimLen, req_start, req_count );

  for( j  =  0; j< input_size; j++) {
    if( input_buff[j] !=  uncompressed_data[j])  {
      printf("Data mismatches at %d \n", j); exit(1); 
    }
  }
    
  printf("Data is matched\n");


  printf("%s, compress_algo  =  %d, compress_perform  =  %f\n", varname, 
	 compress_algo, compress_perform);

  index  =  find_index(varname, &idx);
  update_hash_table(index, compress_algo, compress_perform);   

  if ( update_log) {	
    write_analysis_log(); 
  }

  free(chunks);
    
  return actual_output_size;

}

void acomp_get_actual_size( unsigned char* metadata, unsigned int *uncompressed_size_meta ) {
  *uncompressed_size_meta  =  *((unsigned int*)metadata);
}

int acomp_is_compressed(unsigned char* metadata, int* compression_method) {
  char compress_ok  =  *((unsigned char*)(metadata + sizeof(unsigned int)));  
  *compression_method  =  (compress_ok>>1) & 0x0f;
  compress_ok &=  1;
  return compress_ok;
}

unsigned int acomp_get_max_metadata_size() {
  return (sizeof(unsigned int) + 3 * sizeof(unsigned char) + sizeof(unsigned int) *  MAX_CHUNKS);    // metadata: original data size (unsigned int) + compression succ flag (char)

}

unsigned char acomp_get_byteflags(unsigned char* metadata) {
  return *((unsigned char *)(metadata + sizeof(unsigned int) + 1));
    
}

int acomp_get_nchunks(unsigned char* metadata, unsigned int** chunk_meta) {
  *chunk_meta  =  (unsigned int *)((unsigned char*) metadata + sizeof(unsigned int) + 3);
  return *((unsigned char*)(metadata + sizeof(unsigned int) +  2));  
}

int acomp_decompress(int t_size, int ndims, unsigned char* metadata, char* compressed_data, 
		     char* output_data, 
		     unsigned int* dimLen, unsigned int* req_start, unsigned int* req_count ) {

  int compression_method, compress_ok ;
  unsigned char byteFlags, blg;
  char *uncomp_data, *uncomp_cdata;
  unsigned int* chunk_meta  =  NULL;
  char *prep_data  =  NULL, *chunk_out_data  =  NULL;
  int i  =  0, nchunks, noverlap, num_comp_bytes  =  0;
  unsigned int uncompressed_size; 
  int chunk_index;
  unsigned int comp_size, chunk_offset, safe_size  =  4096;  
  struct chunk_info_t* chunks;
    
  compress_ok  =  acomp_is_compressed( metadata, &compression_method);
  printf("compress_ok %d compression_method %d \n", compress_ok, compression_method); fflush(stdout);
  if(compress_ok !=  0) {

    byteFlags  =  acomp_get_byteflags(metadata);
    nchunks  =  acomp_get_nchunks(metadata, &chunk_meta);
    acomp_get_actual_size( metadata, &uncompressed_size);

      
    
    blg  =  byteFlags;
      
    printf("uncompressed_size %u, compress_ok %d, byteflags %d, nchunks %d \
compression method %d num_comp_bytes %d\n", uncompressed_size, compress_ok, byteFlags, nchunks, 
	   compression_method, num_comp_bytes);  	    
     
          
    //output_data  =  malloc(uncompressed_size);
      
    assert(output_data !=  NULL);
	
    while(i < t_size) {
      if((blg >> i)& 1) num_comp_bytes++;
      i++;    
    }

    noverlap  =  get_overlap_chunks_1d(&chunks, ndims, dimLen, nchunks, 
				     req_start, req_count, t_size);

    //noverlap  =  1;

    printf("N compressible bytes %d, noverlap %d\n", num_comp_bytes, noverlap);

    uncomp_data  =  compressed_data;

    compressed_data +=  (uncompressed_size/t_size) * (t_size - num_comp_bytes);

    assert(noverlap > 0);

    chunk_index  =  chunks[0].chunk_id;
   
    if( chunk_index  ==  0)
      chunk_offset  =  0;
    else 
      chunk_offset  =  chunk_meta[chunk_index - 1];
	
    if(chunks[noverlap-1].chunk_size > chunks[0].chunk_size)
      prep_data  =  malloc(4* chunks[noverlap-1].chunk_size);	     
    else if(chunks[0].chunk_size > 1024)
      prep_data  =  malloc(4* chunks[0].chunk_size);
    else 	
      prep_data  =  malloc(4* safe_size);

    assert(prep_data !=  NULL);

    chunk_out_data  =  output_data;

    uncomp_cdata  =  uncomp_data;		     

    for( i  =  0; i < noverlap ; i++) {
      chunk_index  =  chunks[i].chunk_id;
	
      comp_size  =  chunk_meta[chunk_index] - chunk_offset;

      printf("chunk %d, chunkoffset %u, comp_size %u chunk_size %u\n", chunks[i].chunk_id, chunk_offset, comp_size, chunks[i].chunk_size); 
      printf("chunk %d, chunk output offset %u\n", chunks[i].chunk_id, chunks[i].chunk_offset); 

      	
      chunk_out_data  =  output_data + chunks[i].chunk_offset;

      decompress_prep ( prep_data, compressed_data + chunk_offset, comp_size, 
			(chunks[i].chunk_size<1024?1024:chunks[i].chunk_size), 
			compression_method % T_PREPROCESS); 
    
      printf("New pointer %u\n",(chunks[i].chunk_size/t_size)*(t_size - num_comp_bytes)); fflush(stdout);

         
      if( compression_method <=  T_METHODS)
	process_data( byteFlags, chunk_out_data, prep_data, uncomp_cdata, 
		      t_size, chunks[i].chunk_size/t_size, 
		      compression_method % T_PREPROCESS + T_PREPROCESS);
      else 
	memcpy(chunk_out_data, prep_data, chunks[i].chunk_size);

      if(num_comp_bytes) {
	uncomp_cdata +=  (chunks[i].chunk_size/t_size)*(t_size - num_comp_bytes); 
      }

      chunk_offset  =  chunk_meta[chunk_index];

    }

    free(prep_data);
    free(chunks);
    return 1;
  }  

  return 0;
}

