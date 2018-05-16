#include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "acomp.h"

/** Included files for hashing, analysis, compressing
chuncking etc **/
#include "hash_util.h"
#include "byte_analyzer.h"
#include "preconditioners.h"
#include "compressor.h"
#include  "chunking.h"

/*-- Global declarations ------*/
#define MAX_CMP_ALGO 15

#define MAX_STRING_LEN 50

#define COMPRESSION_SCHEME(i, j, z) { char* x = PRECOND(i);            \
				      char* y = COMPRESSOR(j);         \
				      sprintf(z, "%s_%s", x, y); } 

#define T_METHODS (T_PREPROCESS * T_COMP_ALGOS) 

#define COMP_BY_BYTECOL 0

double weight_cr =  DEFAULT_W_CR;
double weight_cs =  DEFAULT_W_CS;
double delta =  DEFAULT_DELTA;
char *c_var = NULL;
int acomp_init = 0;

bool_t has_changed = 0; 

char C_ALGO[T_METHODS][MAX_STRING_LEN];

/* ---------------------------- */

int analyzeData (char *varname, char *input_data, acomp_size_t data_size,
		int is_real_type, int ele_size, acomp_size_t total_bytes, int* compress_algo, 
		double* cPerform, unsigned char* byteFlags, int *n_comp_bytes) {

  int i, j, index, idx;
     
  byte_t *comp_buffer, *uncomp_buffer;
  byte_t *prep_data;
  double analysis_data[MAX_CMP_ALGO];  

  int bAcnt  =  0;

  index  =  find_index(varname, &idx);

  if( index  ==  -1) {
    index  =  idx; 
    assert(index !=  -1); 
    strcpy(hash_table[index].varname, varname);
    hash_table[index].is_empty  =  1;
  }
  
  c_var = varname;

  DB_PRINTF("Variable : %s found at index %d \n", varname, index); fflush(stdout);
  *compress_algo  =  hash_table[index].compress_algo;
  
  ubyte_t byteArray  =  get_compressibility_flags(is_real_type , ele_size, 
	         			 total_bytes, (byte_t*) input_data, &bAcnt);

  *n_comp_bytes  =  bAcnt;
  *byteFlags  =  byteArray;

#if 1
  if( ( *compress_algo !=  -1 ) && 
	 ~( hash_table[index].comp_perform < hash_table[index].avg_comp_perform + delta &&
	    hash_table[index].comp_perform > hash_table[index].avg_comp_perform - delta ) ) {
    return 0;
  }    
#endif
  
      
  DB_PRINTF("Analyze data: input_data %p, total_bytes %llu, data_size %llu type_size %llu\n", 
					       input_data, total_bytes, data_size, ele_size);
  DB_PRINTF("analyze data: byteArray %d\n", byteArray);
  DB_PRINTF("weight cr %g, weight cs %g\n", weight_cr, weight_cs);

  acomp_size_t cSize;

  /** Test size is 10-15% of the original size **/
  acomp_size_t testSize = data_size/8;

  DB_PRINTF("analyze data : bAcnt %d\n", bAcnt); 

  /** Is this required?? **/
  acomp_size_t size  =  total_bytes > 1024 ? 2 * total_bytes: 1024;

  comp_buffer  =  (byte_t *) malloc(size * sizeof(byte_t));
  uncomp_buffer  =  (byte_t*) malloc(size * sizeof(byte_t));
  prep_data  =  (byte_t*) malloc(size * sizeof(byte_t));

  //ubyte_t* output_data  =  malloc(size * sizeof(ubyte_t));

  if(comp_buffer  ==  NULL || uncomp_buffer  ==  NULL || prep_data  ==  NULL) { 
    fprintf (stderr, "ADIOS ERROR: failed while allocating  \
			 comp_buffer and uncomp_buffer in analyzeData() \n");
    exit(-1);
  }

  struct timeval start_time, end_time;
  double pTime, cTime, cSpeed, cRatio, res;  

  memcpy(prep_data, input_data, testSize);

  *cPerform  =  10000.0f;
  *compress_algo  =  -1;

#ifdef PRINT_STATS
  FILE *file = fopen("comp_indv_summary", "a+" ); 
  fprintf(file, "%s:", varname);
#endif 

  for ( i  =  0; i < T_PREPROCESS; i++ ) {

    gettimeofday(&start_time, NULL);
    preprocess_data(byteArray, input_data, prep_data, uncomp_buffer, 
					       ele_size, testSize, i);
    gettimeofday(&end_time, NULL);
    
    pTime  =  (double)(end_time.tv_usec - start_time.tv_usec) / 1000000 
		+ (double)(end_time.tv_sec - start_time.tv_sec);


    for ( j  =  0; j < T_COMP_ALGOS; j++) {

      gettimeofday(&start_time, NULL);    

#if COMP_BY_BYTECOL == 0
      cSize  =  compress_prep(prep_data, comp_buffer, testSize * bAcnt, j);
#else
      cSize  =  compress_prep_by_bytecol(prep_data, comp_buffer, bAcnt, 
							      testSize, j);
#endif

      gettimeofday(&end_time, NULL); 

      cTime  =  (double)(end_time.tv_usec - start_time.tv_usec) / 1000000 
		+ (double)(end_time.tv_sec - start_time.tv_sec);

      cTime += pTime;

      if( cSize >= testSize * bAcnt )
		cTime = 100000; /** Do not use this compression scheme if 
					it doesn't compress the data **/

      cSpeed  =  ((double) (testSize * (bAcnt)) )/cTime;	
      cSpeed  =  cSpeed/(1024*1024); 

      //cRatio  =  (double)( testSize * bAcnt )/cSize;
      cRatio  =  (double)(testSize * ele_size)/(cSize + testSize * (ele_size - bAcnt));

      res  =  ((cSize + testSize * (ele_size - bAcnt)) * IO_BANDWIDTH ) * weight_cr 
		+ cTime * weight_cs;

      analysis_data[ i * T_COMP_ALGOS + j]  =  res; 
            
      if( *cPerform  > res )  { 
	*cPerform  =  analysis_data[ i * T_COMP_ALGOS + j ];
	*compress_algo  =  i * T_COMP_ALGOS + j;
      }

#ifdef PRINT_STATS      
      fprintf(file, "%s[%d]: %g (%u), %g(%g), %g <==>", C_ALGO[ i * T_COMP_ALGOS+ j ], 
	                          (i * T_COMP_ALGOS+ j), cRatio, cSize, cSpeed, cTime,   
 		    			        analysis_data[ i * T_COMP_ALGOS+ j ]);
#endif
      
    }

  }
#ifdef PRINT_STATS      
  fprintf(file,"\n");
  fclose(file);
#endif
      

  if( *compress_algo != hash_table[index].compress_algo)
       has_changed = 1;


#ifdef PRINT_STATS
      
  FILE * comp_stats = fopen("comp_summary", "a+");
	fprintf(comp_stats,
	"(%s) => [%s] %d, %g\n", has_changed? "Changed": "", varname, 
			*compress_algo, analysis_data[*compress_algo]);
  fclose(comp_stats);
#endif

  DB_PRINTF("%s: compression algo %d, compression perform %g\n", varname, *compress_algo, *cPerform);      

  return 1;
}


void verify(ubyte_t byteArray, ubyte_t** bSeg, ubyte_t** dbSeg, int ele_size, acomp_size_t data_size) {
  /* verify */

  ubyte_t* dbSegment, *bSegment;

  int byte;
  acomp_size_t varIdx;
   
  for ( varIdx  =  0; varIdx < data_size; ++varIdx) {
    int idx  =  varIdx;
    for ( byte = 0; byte < ele_size; ++byte) {
      dbSegment  =  dbSeg[byte];
      bSegment  =  bSeg[byte];

      if ( ( byteArray >> byte) & 1) {
	if (bSegment[idx] !=  dbSegment[idx]) {
	  fprintf(stderr, "Mismatch byte : %d varIdx: %lu\n", byte, varIdx);
	  exit(1);
	}
      }
    }
  }
   
}


void verify_compression(ubyte_t byteArray, ubyte_t* input_data, 
			ubyte_t* output_data, int ele_size, 
			acomp_size_t data_size) {
  /* verify */

  ubyte_t* dbSegment, *bSegment;

  int byte;
  acomp_size_t varIdx;
  dbSegment  =  output_data;
  bSegment  =  input_data;

  for ( varIdx  =  0; varIdx < data_size; ++varIdx) {
    acomp_size_t idx  =  varIdx * ele_size;
    for ( byte = 0; byte < ele_size; ++byte) {
      if ( ( byteArray >> byte) & 1) {
	if (bSegment[idx + byte] !=  dbSegment[idx + byte]) {
	  fprintf(stderr, "Mismatch at byte: %d varIdx: %lu\n", byte, varIdx);
	  exit(1);
	}
      }
    }
  }
   
}

void acomp_conf(double w_cs, double w_cr, double del) {
  int i = 0, j = 0;
  weight_cr  =  w_cr;
  weight_cs  =  w_cs;
  DB_PRINTF("weight cr %f \n", weight_cr); fflush(stdout);
  DB_PRINTF("weight cs %f \n", weight_cs); fflush(stdout);
  delta  =  del;
  for ( i = 0; i < T_PREPROCESS; i++) {
    for ( j = 0; j < T_COMP_ALGOS; j++ ) {
        char * str = &C_ALGO[i * T_PREPROCESS + j ][0];
	COMPRESSION_SCHEME(i, j, str);
    }
  }
  acomp_init = 1;	  
}


acomp_size_t acomp_compress( char * varname, int ndims, acomp_size_t *input_dim_len, 
			     int ele_size, int is_real_type, 
			     acomp_size_t input_size_bytes, acomp_size_t n_elements,  
			     char* output_buff,  char* input_buff, 
			     unsigned char * metadata) {

  int chunk_id  =  0, nchunk  =  0, idx, index;
  struct chunk_info_t* chunks;
  ubyte_t byteFlags;
  int  n_comp_bytes;
  double compress_perform  =  0, total_time = 0;
  struct timeval start_time1, end_time1;

  int update_log  =  0, compress_algo  =  -1;
  char *comp_buff  =  NULL, *uncomp_buff =  NULL;
  ubyte_t compress_ok  =  1;
  size_t actual_output_size  =  0; //output_size;

  has_changed = 0;
  gettimeofday(&start_time1, NULL);
  
  if( acomp_init == 0) {
    acomp_conf(0, 1.0 , delta);
  }	 

  if(read_analysis_log()  ==  -1) 
    initialize_hash();

  STAT_PRINTF("********* COMPRESSING VAR %s ***************** \n", varname)

  update_log  =  analyzeData(varname, input_buff, n_elements, is_real_type, ele_size, 
			   input_size_bytes, &compress_algo, 
			   &compress_perform, &byteFlags, &n_comp_bytes);
        
  DB_PRINTF("N compressible bytes %d of %d bytes",n_comp_bytes, ele_size);
 
  comp_buff  =  output_buff + ((unsigned int)(ele_size - n_comp_bytes)) * n_elements; 
      
  uncomp_buff  =  output_buff; 
  nchunk  =  get_chunks_info_1D(&chunks, input_dim_len, ndims, ele_size);


  struct timeval start_time, end_time;
  double cTime, total_comp_time = 0, compression_speed  =  0, compression_ratio;
  byte_t* prep_data  =  (byte_t*) malloc(input_size_bytes * 4);     
     
  while ( chunk_id < nchunk && compress_ok) {
 
    chunks[chunk_id].data  =  input_buff + chunks[chunk_id].chunk_offset; 
    DB_PRINTF("Processing chunk %d[of %d] offset %lu, buffer pointer %p\n", chunk_id, nchunk,
				chunks[chunk_id].chunk_offset, chunks[chunk_id].data);
          

    acomp_size_t subset_size  =  chunks[chunk_id].chunk_size / ele_size; 
    acomp_size_t compressed_size;
    int num_comps  =  n_comp_bytes;

      
    gettimeofday(&start_time, NULL); 
    if(compress_algo <=  T_METHODS )
      preprocess_data(byteFlags, chunks[chunk_id].data, prep_data, uncomp_buff,
		   ele_size, subset_size, compress_algo / T_PREPROCESS);
    else {
      memcpy(prep_data, chunks[chunk_id].data, subset_size * ele_size);
      num_comps  =  ele_size;
    } 	   
    
    DB_PRINTF("Data size %lu total_dsize %lu\n", subset_size, subset_size * num_comps);

#if COMP_BY_BYTECOL == 0 
    compressed_size  =  compress_prep ( prep_data, comp_buff, subset_size * num_comps, 
				      compress_algo % T_PREPROCESS);
#else
    compressed_size  =  compress_prep_by_bytecol ( prep_data, comp_buff, num_comps, subset_size, 
				      compress_algo % T_PREPROCESS);

#endif  
    gettimeofday(&end_time, NULL); 


#if 1
    if( compressed_size >=  subset_size * num_comps) {

      compressed_size  =  0;
      compress_ok  =  0;
 
    }
#endif
   
    chunks[chunk_id].compSize  =  compressed_size;


    cTime  =  (double)(end_time.tv_usec - start_time.tv_usec) / 1000000 
      + (double)(end_time.tv_sec - start_time.tv_sec);

    total_comp_time +=  cTime;  
    DB_PRINTF("Done compression and compression size = %lu\n", compressed_size);
    actual_output_size +=  compressed_size;
    comp_buff +=  compressed_size;
    uncomp_buff +=  (chunks[chunk_id].chunk_size/ele_size) *(ele_size - n_comp_bytes);
    chunk_id ++;
  }
  free(prep_data);

  //compression_ratio  =  ((double)(n_elements * n_comp_bytes))/actual_output_size;
  compression_ratio  =  ((double)(n_elements * ele_size))/(actual_output_size + 
  			n_elements * ( ele_size - n_comp_bytes));

  compression_speed  =  ((double) (n_elements * n_comp_bytes))/(1024*1024);

  compression_speed  =  compression_speed/total_comp_time;

  //compress_perform  =  weight_cr * compression_ratio + weight_cs * compression_speed;
  compress_perform  =  weight_cr * IO_BANDWIDTH *
		       (actual_output_size + n_elements * ( ele_size - n_comp_bytes))  
			+ weight_cs * total_comp_time;

  DB_PRINTF("done compressing. Compressed size  =  %lu\n", actual_output_size);

  if(compress_ok !=  1) {
    actual_output_size  =  input_size_bytes;
  }
  else {
    actual_output_size +=  (ele_size - n_comp_bytes) * n_elements;
  }
   
  DB_PRINTF("DEBUG: Compression status[%d]\n ", compress_ok);
  
  compress_ok  =  compress_algo << 1  | compress_ok;


  // copy the metadata, simply the original size before compression
  ubyte_t *metadata1  =  metadata; 

  if( metadata )
    {   int i;
      ubyte_t nc  =  nchunk;
      acomp_size_t total_compress  =  0;	

      memcpy((ubyte_t*)metadata, &input_size_bytes, sizeof(acomp_size_t));
      metadata +=  sizeof(acomp_size_t);

      memcpy((ubyte_t*)metadata, &compress_ok, sizeof(ubyte_t));
      metadata +=  sizeof(ubyte_t);

      memcpy((ubyte_t*)metadata, &byteFlags, sizeof(ubyte_t));
      metadata +=  sizeof(ubyte_t);

      memcpy((ubyte_t*)metadata, &nc, sizeof(ubyte_t));
      metadata +=  sizeof(ubyte_t);

      for(i  =  0; i< nchunk; i++) {
	DB_PRINTF("chunk id %d, chunk offset %lu, chunk size %lu\n", i, 
				total_compress, chunks[i].compSize);

	total_compress +=  chunks[i].compSize;

	memcpy((ubyte_t*)metadata, &total_compress, sizeof(acomp_size_t));

	DB_PRINTF("written compression size %lu\n", *(acomp_size_t*)metadata);
	metadata +=  sizeof(acomp_size_t);

      }
    }

#ifdef VERIFY    
  acomp_size_t j  =  0;
  byte_t* uncompressed_data = calloc( 4 * input_size_bytes, sizeof(byte_t));
  acomp_size_t dimLen[5], req_start[5], req_count[5];

  for(j = 0; j < 5; j++) {
    dimLen[j] = 1;
    req_start[j] = 0;
    req_count[j] = 1;
  }			

  dimLen[ndims-1]  =  input_size_bytes/ele_size; 

  req_count[ndims-1]  =  input_size_bytes/ele_size; 

  DB_PRINTF("Ndims %d, %u, DimLen[%d] = %u, req_count[%d] = %u\n", ndims, 
				input_size_bytes/ele_size ,ndims-1, dimLen[ndims-1], 
						ndims -1, req_count[ndims-1]);
  acomp_decompress(ele_size, ndims, metadata1, output_buff,
		   uncompressed_data, dimLen, req_start, req_count );

  for( j  =  0; j< input_size_bytes; j++) {

    if( input_buff[j] !=  uncompressed_data[j])  {
      printf("Data mismatches at index %lu \n", j); 
      exit(1); 
    }
  }
    
  printf("Verification is successful!\n");
#endif

  STAT_PRINTF("Compress_algo used is %d, weight_cr %g, compression ratio %g, weight_cs %g, \
               compression_speed %g, compression time %g, compress_perform  =  %f\n", 
	       compress_algo, weight_cr, compression_ratio, weight_cs, compression_speed,
							 total_comp_time, compress_perform);




  index  =  find_index(varname, &idx);
  update_hash_table(index, compress_algo, compress_perform);   

  if ( update_log) {	
    write_analysis_log(); 
  }

  free(chunks);
  gettimeofday(&end_time1, NULL);
  total_time  =  (double)(end_time1.tv_usec - start_time1.tv_usec) / 1000000 + 
				(double)(end_time1.tv_sec - start_time1.tv_sec);
#ifdef PLOT
   FILE *statfile2 = fopen("stat.txt", "a+");
   fprintf(statfile2, "WCR %g WCS %g CA(%d) CR %g CS %g CP %g CT %g CSZ %u Total_T %g\n", weight_cr, weight_cs,compress_algo, compression_ratio, 
				compression_speed, compress_perform, total_comp_time, actual_output_size, total_time);
   fprintf(statfile2, "%g\n", compress_perform);
   fclose(statfile2);
#endif

  STAT_PRINTF("**********%s*******************************\n", "****" );
    
  return actual_output_size;

}

void acomp_get_actual_size( ubyte_t* metadata, acomp_size_t *uncompressed_size_meta ) {
  *uncompressed_size_meta  =  *((acomp_size_t*)metadata);
}

int acomp_is_compressed(ubyte_t* metadata, int* compression_method) {
  byte_t compress_ok  =  *((byte_t*)(metadata + sizeof(acomp_size_t)));  
  //*compression_method  =  (compress_ok>>1) & 0x0f;
  *compression_method  =  (compress_ok>>1);
  compress_ok &=  1;
  return compress_ok;
}

int acomp_get_max_metadata_size() {
  return (3 * sizeof(ubyte_t) + sizeof(acomp_size_t) * (MAX_CHUNKS + 1));
}

ubyte_t acomp_get_byteflags(ubyte_t* metadata) {
  return *((ubyte_t *)(metadata + sizeof(acomp_size_t) + 1));
    
}

int acomp_get_nchunks(ubyte_t* metadata, acomp_size_t** chunk_meta) {
  *chunk_meta  =  (acomp_size_t *)((ubyte_t*) metadata + sizeof(acomp_size_t) + 3);
  return *((ubyte_t*)(metadata + sizeof(acomp_size_t) +  2));  
}

int acomp_decompress(int ele_size, int ndims, unsigned char* metadata, char* compressed_data, 
		     char* output_data, acomp_size_t* dimLen, acomp_size_t* req_start, acomp_size_t* req_count ) {

  int compression_method, compress_ok ;
  ubyte_t byteFlags, blg;
  byte_t *uncomp_data, *uncomp_cdata;
  acomp_size_t *chunk_meta  =  NULL;
  byte_t *prep_data  =  NULL, *chunk_out_data  =  NULL;
  int i  =  0, nchunks, noverlap, num_comp_bytes  =  0;
  acomp_size_t uncompressed_size; 
  int chunk_index;
  acomp_size_t comp_size, chunk_offset, safe_size  =  4096;  
  struct chunk_info_t* chunks;
    
  compress_ok  =  acomp_is_compressed( metadata, &compression_method);

  DB_PRINTF("compress_ok %d compression_method %d \n", compress_ok, compression_method)

  if(compress_ok !=  0) {

    byteFlags  =  acomp_get_byteflags(metadata);
    nchunks  =  acomp_get_nchunks(metadata, &chunk_meta);
    acomp_get_actual_size( metadata, &uncompressed_size);

    blg  =  byteFlags;
      
    DB_PRINTF("uncompressed_size %lu, compress_ok %d, byteflags %d, nchunks %d \
               compression method %d num_comp_bytes %d\n", uncompressed_size, 
	       compress_ok, byteFlags, nchunks, compression_method, num_comp_bytes)  	    
          
      
    assert(output_data !=  NULL);
	
    while(i < ele_size) {
      if((blg >> i)& 1) num_comp_bytes++;
      i++;    
    }

    noverlap  =  get_overlap_chunks_1D(&chunks, ndims, dimLen, nchunks, 
				     req_start, req_count, ele_size);

    DB_PRINTF("N compressible bytes %d, noverlap %d\n", num_comp_bytes, noverlap);

    uncomp_data  =  compressed_data;

    compressed_data +=  (uncompressed_size/ele_size) * (ele_size - num_comp_bytes);

    assert(noverlap > 0);

    chunk_index  =  chunks[0].chunk_id;
   
    if( chunk_index  ==  0)
      chunk_offset  =  0;
    else 
      chunk_offset  =  chunk_meta[chunk_index - 1];

	
    if( chunks[ noverlap - 1].chunk_size > chunks[0].chunk_size)
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

      
      DB_PRINTF("chunk %d, chunkoffset %lu, comp_size %lu chunk_size %lu\n", 
	chunks[i].chunk_id, chunk_offset, comp_size, chunks[i].chunk_size); 

      DB_PRINTF("chunk %d, chunk output offset %lu\n", chunks[i].chunk_id, 
						   chunks[i].chunk_offset); 
      	
      chunk_out_data  =  output_data + chunks[i].chunk_offset;

#if COMP_BY_BYTECOL == 0
      decompress_prep ( prep_data, compressed_data + chunk_offset, comp_size, 
			(chunks[i].chunk_size<1024?1024:chunks[i].chunk_size), 
			compression_method % T_PREPROCESS); 
#else
      decompress_prep_by_bytecol ( prep_data, compressed_data + chunk_offset, 
			num_comp_bytes, (chunks[i].chunk_size)/ ele_size, 
			compression_method % T_PREPROCESS); 

#endif
    
      DB_PRINTF("New pointer %lu\n",
	   (chunks[i].chunk_size/ele_size)*(ele_size - num_comp_bytes)) 

         
      if( compression_method <=  T_METHODS)
	postprocess_data( byteFlags, chunk_out_data, prep_data, uncomp_cdata, 
		      ele_size, chunks[i].chunk_size/ele_size, 
		      compression_method / T_PREPROCESS);
      else 
	memcpy(chunk_out_data, prep_data, chunks[i].chunk_size);

      uncomp_cdata +=  (chunks[i].chunk_size/ele_size)*(ele_size - num_comp_bytes); 

      chunk_offset  =  chunk_meta[chunk_index];

    }

    free(prep_data);
    free(chunks);
    return 1;
  }  

  return 0;
}

