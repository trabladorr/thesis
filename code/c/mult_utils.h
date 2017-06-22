#ifndef _MULT_UTILS_H
#define _MULT_UTILS_H

void printNxMbytes(void *bytes, uint8_t n, uint8_t m);
void printbits(void *bytes, uint8_t num);

void create_16_vector_byteslice(__m128i *v);
void create_128_vector_bitslice(__m128i *v);
void create_128_vector_serial(__m128i *v);
void create_16_vector_serial(__m128i *v);

void use_m128i_variable(__m128i v,int zero);
void use_m128i_variables(__m128i v[],int zero);

#define PRINT_BUFF 1000

#endif