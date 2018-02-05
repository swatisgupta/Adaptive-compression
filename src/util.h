#ifndef ADAPTIVE_UTIL_H

#define ADAPTIVE_UTIL_H

#include <sys/stat.h>
#include <string.h>
#include <assert.h>

#define MAX_HASH_ENTRY 200


int inithash = 0;
int rank;

struct analysis_t {
   
   char varname[40];
   float avg_comp_perform;	
   float comp_perform;	
   int compress_algo;
   int tsteps;
   int is_empty;

} hash_table[MAX_HASH_ENTRY];


long hash(char* key) {
	long hashVal = 0;
	while (*key != '\0') {
		hashVal = (hashVal << 4) + *(key++);
		long g = hashVal & 0xF0000000L;
		if (g != 0) hashVal ^= g >> 24;
		hashVal &= ~g;
	}
	return hashVal;
}		


int read_analysis_log() {

  char log[100];
  sprintf(log, ".analysis%d.log", rank);	
  FILE * fp = NULL;
  struct stat buffer;
  int exist = stat(log, &buffer);

  if(exist != 0)
        return -1;

  fp = fopen(log, "rb");
	
  if(fp) {
      //printf("Reading from file %s\n", log);
      fread(hash_table, sizeof(struct analysis_t), MAX_HASH_ENTRY, fp);
      fclose(fp); 
      return 0;
  }

  return -1;

}


void write_analysis_log() {
  char log[100];
  sprintf(log, ".analysis%d.log", rank);	
  FILE * fp = fopen(log, "wb");
   if(fp) {
       fwrite(hash_table, sizeof(struct analysis_t), MAX_HASH_ENTRY, fp);
       fclose(fp); 
    }      	
}


void initialize_hash() {
    int j = 0 ; 
    for( j =0 ; j< MAX_HASH_ENTRY; j++ ) {
	hash_table[j].varname[0] = '\0';
        hash_table[j].avg_comp_perform = 0;    
        hash_table[j].comp_perform = -1;    
        hash_table[j].tsteps = 0;    
        hash_table[j].is_empty = 0;    
        hash_table[j].compress_algo = -1;    
    }
    //printf("done initialize the hash table %d \n", j);
}


int get_hash_index(char *key ) {
   
   long hashcode = hash(key);
   int index = hashcode % MAX_HASH_ENTRY;
   int counter = index;
   if(hash_table[index].is_empty == 0) 
        return index;
   
   do {
 		counter ++;
                counter = counter % MAX_HASH_ENTRY;
    } while(counter != index && hash_table[counter].is_empty == 0);

   if(counter != index)  return counter;

   return -1;      
} 

void update_hash_table(int index, int compress_algo, int compress_perform) {
        hash_table[index].avg_comp_perform = 
			hash_table[index].avg_comp_perform * hash_table[index].tsteps + compress_perform;    
        hash_table[index].tsteps ++;    
	hash_table[index].avg_comp_perform /= hash_table[index].tsteps;
        hash_table[index].comp_perform = compress_perform;    
        hash_table[index].compress_algo = compress_algo;    

}

int find_index(char* key, int *potential_index) {
   unsigned long hashcode = hash(key);
   int index = hashcode % MAX_HASH_ENTRY;
   int counter = index;
   //printf("Key %s Index %d\n", key, index);

   assert(potential_index != NULL); 
   *potential_index = -1;  

   if(strcmp(hash_table[index].varname, key) == 0)
	return index;

   if(hash_table[index].is_empty == 0) 
         *potential_index = index;
   
   do {
         counter ++;
         counter = counter % MAX_HASH_ENTRY;
  	 if(*potential_index == -1 && hash_table[counter].is_empty == 0) 
             *potential_index = counter;	
 
         //printf("Key %s Varname %s Counter %d\n", key, hash_table[counter].varname, counter);
   } while(counter != index && strcmp(hash_table[counter].varname, key) != 0);
        
   if(counter != index)  return counter;
    
   return -1;
}

#endif
