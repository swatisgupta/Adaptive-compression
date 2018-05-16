#include <zlib.h> 
#include <bzlib.h>
#include <lzo/lzoconf.h>
#include <lzo/lzo1x.h>
#include "preconditioners.h"


void preprocess_bSeg(ubyte_t byteArray, int ele_size, acomp_size_t data_size, byte_t* comp_buffer, 
		     byte_t* uncomp_buffer, byte_t* in_buffer) {

  int startOffset = 0, byte = 0, byteCount  =  0, ubyteCount  =  0;
  acomp_size_t varid  =  0;

  for ( varid  =  0; varid < data_size ; ++varid ) {
    acomp_size_t idx  =  startOffset + varid * ele_size;
    for( byte  =  0; byte < ele_size; ++byte) {
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

void preprocess_bWiseSegXor(ubyte_t byteArray, int ele_size, acomp_size_t data_size, byte_t* comp_buffer, 
			    byte_t* uncomp_buffer, byte_t* in_buffer) {

  acomp_size_t varid  =  0; 
  int byte_pos[8], byte, nextOffset  =  0;
  int ncb_count  =  0, bcount  =  0;
    
 //printf("Eleement size %d\n", ele_size); fflush(stdout); 
       
  for ( byte  =  0; byte < ele_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount;
      bcount ++;	 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    DB_PRINTF("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
  }
  DB_PRINTF("Non compressible recount %d\n", ncb_count);	

  for( varid  =  0; varid < data_size; varid ++) {
    for( byte  =  0; byte < ele_size; ++ byte) {
      acomp_size_t curr_index  =  nextOffset + varid * ele_size + byte;
      acomp_size_t next_index  =  nextOffset + (varid + 1) * ele_size + byte;
      acomp_size_t bindex  =  byte_pos[byte] * data_size + varid;
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

void preprocess_bWiseSeg(ubyte_t byteArray, int ele_size, acomp_size_t data_size, byte_t* comp_buffer, 
			 byte_t* uncomp_buffer, byte_t* in_buffer) {
  acomp_size_t varid  =  0;
  int byte_pos[8], byte  =  0; 
  int ncb_count  =  0, bcount  =  0;
     
  unsigned int nextOffset  =  0;
  nextOffset  =  0;

    
 //printf("Eleement size %d\n", ele_size); fflush(stdout); 
  for ( byte  =  0; byte < ele_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount; ++bcount; 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    DB_PRINTF("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
    fflush(stdout); 
  }
  for( varid  =  0; varid < data_size; varid ++) {
    for( byte  =  0; byte < ele_size; ++ byte) {
      acomp_size_t curr_index  =  nextOffset + varid * ele_size + byte;
      acomp_size_t bindex  =  byte_pos[byte] * data_size + varid;
      if(( byteArray >> byte) & 1) { 
	comp_buffer[bindex]  =   in_buffer[curr_index];
      } else {
	uncomp_buffer[bindex]  =  in_buffer[curr_index];
      }
    }
  }

}

void postprocess_bSeg(ubyte_t byteArray, int ele_size, acomp_size_t data_size, byte_t* comp_buffer, 
		      byte_t* uncomp_buffer, byte_t* in_buffer) {
  acomp_size_t startOffset  =  0, byteCount  =  0, ubyteCount  =  0;
  acomp_size_t varid  =  0;
  int byte  =  0; 
	
  for ( varid  =  0; varid < data_size; ++varid) {
    acomp_size_t idx  =  startOffset + varid * ele_size;
    for( byte  =  0; byte < ele_size; ++byte) {
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


void postprocess_bWiseSeg(ubyte_t byteArray, int ele_size, acomp_size_t data_size, byte_t* comp_buffer, 
			  byte_t* uncomp_buffer, byte_t* in_buffer) {

  acomp_size_t varid  =  0;
  int byte_pos[8], byte = 0, nextOffset  =  0; 
  int ncb_count  =  0, bcount  =  0;


  for ( byte  =  0; byte < ele_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount; ++bcount; 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    DB_PRINTF("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
  }  
	
  for( varid  =  0; varid < data_size; varid ++) {
    for( byte  =  0; byte < ele_size; ++ byte) {
      acomp_size_t curr_index  =  nextOffset + varid * ele_size + byte;
      acomp_size_t bindex  =  byte_pos[byte] * data_size + varid;
      if(( byteArray >> byte) & 1) {
	in_buffer[curr_index]  =  comp_buffer[bindex];
      } else {
	in_buffer[curr_index]  =  uncomp_buffer[bindex];
      }
    }
  }
}

void postprocess_bWiseSegXor(ubyte_t byteArray, int ele_size, acomp_size_t data_size, byte_t* comp_buffer, 
			     byte_t* uncomp_buffer, byte_t* in_buffer) {
  acomp_size_t varid  =  0;
  int byte_pos[8], byte = 0, ncb_count  =  0, bcount  = 0;
  int nextOffset  =  0; 

  for ( byte  =  0; byte < ele_size; ++byte) {
    if (( byteArray >> byte) & 1) { 
      byte_pos[byte]  =  bcount; ++bcount; 
    }
    else { 
      byte_pos[byte]  =  ncb_count; 
      ncb_count++; 
    }
    DB_PRINTF("Byte_pos[%d]  =  %d\n", byte, byte_pos[byte]); 
  }

  for( varid  =  data_size - 1; varid >=  0; varid --) {
     
    for( byte  =  0; byte < ele_size; ++ byte) {
      acomp_size_t curr_index  =  nextOffset + varid * ele_size + byte;
      acomp_size_t next_index  =  nextOffset + (varid + 1) * ele_size + byte;
      acomp_size_t bindex  =  byte_pos[byte] * data_size + varid;
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

void preprocess_data (ubyte_t byteArray, byte_t *in_buffer, byte_t *comp_buffer, byte_t *uncomp_buffer, 
		   int ele_size, acomp_size_t data_size,  preconditioner_t transform_type ) {

  DB_PRINTF("Compressible flags %u\n", byteArray)
 fflush(stdout);	

  switch ( transform_type ) {
  case BSEG:    /**---collecting compressible bytes ---***/
    preprocess_bSeg(byteArray, ele_size, data_size, comp_buffer, 
		    uncomp_buffer, in_buffer);	
    break;

  case BWSEG: /** -- collecting compressible bytes seperately -- **/
    preprocess_bWiseSeg(byteArray, ele_size, data_size, comp_buffer, 
			uncomp_buffer, in_buffer);
			    
    break;        

  case BWXOR:  /** --- collecting compressible bytes seperately + Xor ---**/
    preprocess_bWiseSegXor(byteArray, ele_size, data_size, comp_buffer, 
			   uncomp_buffer, in_buffer);
    break;

  default:
    DB_PRINTF("Error: unsupported preprocess method %s\n", "");        
    //printf("Here...data_size %u ele_size %d \n", data_size, ele_size);
 }
}

void postprocess_data (ubyte_t byteArray, byte_t *in_buffer, byte_t *comp_buffer, byte_t *uncomp_buffer,
                   int ele_size, acomp_size_t data_size, preconditioner_t transform_type ) {

  switch( transform_type) {
   case BSEG:    /**--- uncollecting compressible data ---***/
     postprocess_bSeg(byteArray, ele_size, data_size, comp_buffer, 
		     uncomp_buffer, in_buffer);	
     break;
        
   case BWSEG: /**-----  uncollecting compressible bytes seperately ---***/
    postprocess_bWiseSeg(byteArray, ele_size, data_size, comp_buffer, 
			 uncomp_buffer, in_buffer);
    break;

   case BWXOR:  /** --- uncollecting compressible bytes seperately + Xor ---**/
    DB_PRINTF("Inside BWXOR post : size %lu\n", data_size);
    postprocess_bWiseSegXor(byteArray, ele_size, data_size, comp_buffer, 
			    uncomp_buffer, in_buffer);
    break;

  default:
    DB_PRINTF("Error: unsupported postprocess method %s\n", "");        
    //printf("Here...data_size %u ele_size %d \n", data_size, ele_size);
  }     
}
 

