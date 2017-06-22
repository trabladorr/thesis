#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "x86intrin.h"
#include "mult_utils.h"

//helper printing function
void printNxMbytes(void *bytes, uint8_t n, uint8_t m){
    static uint8_t tmp[PRINT_BUFF];
    memcpy(tmp,bytes,n*m);
    int i=0,j=0;
    for (i=0;i<n;i++){
	    printf("%03d: ", i);
	    for (j=0;j<m;j++){
	        printf("%02x ", tmp[i*m+j]);
	    }
	    printf("\n");
	}
	printf("\n");
}

//helper printing function
void printbits(void *bytes, uint8_t num){
    static uint8_t tmp[PRINT_BUFF];
    memcpy(tmp,bytes,num);
    int i=0,j=0;
    for (i=0;i<num;i++){
    	for (j=7;j>=0;j--){
            char byte = (tmp[i] >> j) & 1;
            printf("%u", byte);
        }
        printf(" ");
    }
    printf("\n");
}

void create_16_vector_byteslice(__m128i *v){
    uint8_t i, values[] = {0x01,0x02,0x03,0x04};
    for (i=0;i<4;i++)
        memset(v+i,values[i],16);
}

void create_128_vector_bitslice(__m128i *v){
    int i, values[] = {0,0,0,0,0,0,0,-1,  0,0,0,0,0,0,-1,0,  0,0,0,0,0,0,-1,-1,  0,0,0,0,0,-1,0,0};
    for (i=0;i<32;i++)
        memset(v+i,values[i],16);
}

void create_128_vector_serial(__m128i *v){
    int j;
    v[0] = _mm_setr_epi8(0x01,0x02,0x03,0x04,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00);
    for (j=1;j<32;j++)
        v[j] = _mm_setr_epi8(0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00);
}

void create_16_vector_serial(__m128i *v){
    uint8_t input[64];
    int i,j;
    for (i=0;i<16;i++)
        for (j=0;j<4;j++)
            input[i*4+j]=j+1;
    memcpy(v,input,64);
}

//prevent compiler from over optimizing code in cycle test, by branching on variable
void use_m128i_variable(__m128i v,int zero){
    if (zero)
        printNxMbytes(&v,1,16);
}

//prevent compiler from over optimizing code in cycle test, by branching on variable
void use_m128i_variables(__m128i v[],int zero){
    if (zero)
        printNxMbytes(v,1,16);
}