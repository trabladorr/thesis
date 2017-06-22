#ifndef _MULT_IO_H
#define _MULT_IO_H

#define TABLE_ANUBIS 0
#define TABLE_OTHERH 1

#define MULT_SHUFFLE 0
#define MULT_X2 1
#define MULT_XOR 2

#define IMPL_INLINE 0
#define IMPL_MACRO 1
#define IMPL_INLINE_OPT 2
#define IMPL_INLINE_OPT2 3

//Represents the parameters needed to choose the right function to test
typedef struct parameters{
	int f; 		   //Galois Field exponent
	int n_vectors; //Number of input vectors in every test
	int n_tests;   //Number of tests to run
	int table;     //Indicates which table (and field irreducible polynomial) will be used
	int mult_type; //Indicates which multiplication algorithm will be used
	int impl_type; //Indicates whether the macro or the inline version will be executed
} parameters;

parameters init(int argc, char **argv);

void mult_io(parameters params, int verbose);

#endif