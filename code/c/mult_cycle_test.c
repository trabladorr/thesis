#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "x86intrin.h"

#include "mult_utils.h"
#include "mult_cycle_test.h"
#include "utils_inline.h"
#include "utils_macro.h"

#include "direct_shufflemult_inline.h"
#include "direct_shufflemult_macro.h"

#include "byteslice_shufflemult_inline.h"
#include "byteslice_x2mult_inline.h"

#include "bitslice_xormult_inline.h"
#include "bitslice_xormult_macro.h"

float tests = 1000000.0;

int main(int argc,char **argv){
	if (argc > 1)
		tests = atoi(argv[1]);

	init_macro_constants();

	uint64_t (*test[DIFFTESTS]) (int zero) = {\
		test_mult_4_shuffle,\
		test_mult_4_shuffle_macro,\
		test_mult_2_shuffle_gf16,\
		test_mult_2_shuffle_gf16_macro,\
		test_mult_16_shuffle_anubis,\
		test_mult_16_shuffle_other,\
		test_mult_16_shuffle_load_anubis,\
		test_mult_16_shuffle_load_other,\
		test_mult_16_x2_anubis,\
		test_mult_16_x2_other,\
		test_mult_16_x2_load_anubis,\
		test_mult_16_x2_load_other,\
		test_mult_128_xor_anubis,\
		test_mult_128_xor_anubis_opt_compiler,\
		test_mult_128_xor_anubis_opt_tmp,\
		test_mult_128_xor_anubis_macro,\
		test_mult_128_xor_load_anubis};	

	char *test_descriptions[] = {\
		"4   vectors, GF(2^8) , shuffle, Anubis, inline,     direct",\
		"4   vectors, GF(2^8) , shuffle, Anubis,  macro,     direct",\
		"2   vectors, GF(2^16), shuffle, Anubis, inline,     direct",\
		"2   vectors, GF(2^16), shuffle, Anubis,  macro,     direct",\
		"16  vectors, GF(2^8) , shuffle, Anubis, inline, bytesliced",\
		"16  vectors, GF(2^8) , shuffle, OtherH, inline, bytesliced",\
		"16  vectors, GF(2^8) , shuffle, Anubis, inline, w/ loading",\
		"16  vectors, GF(2^8) , shuffle, OtherH, inline, w/ loading",\
		"16  vectors, GF(2^8) , x2     , Anubis, inline, bytesliced",\
		"16  vectors, GF(2^8) , x2     , OtherH, inline, bytesliced",\
		"16  vectors, GF(2^8) , x2     , Anubis, inline, w/ loading",\
		"16  vectors, GF(2^8) , x2     , OtherH, inline, w/ loading",\
		"128 vectors, GF(2^8) , XOR    , Anubis, inline,  bitsliced",\
		"128 vectors, GF(2^8) , XOR    , Anubis, optCmp,  bitsliced",\
		"128 vectors, GF(2^8) , XOR    , Anubis, optTmp,  bitsliced",\
		"128 vectors, GF(2^8) , XOR    , Anubis,  macro,  bitsliced",\
		"128 vectors, GF(2^8) , XOR    , Anubis, inline, w/ loading"};

	int num_vectors[] = {4, 4, 2, 2, 16, 16, 16, 16, 16, 16, 16, 16, 128, 128, 128, 128, 128};

	uint64_t avg[DIFFTESTS];
	uint64_t tmp;
	int i,j;

	
	//Run tests tests times in silent, to acquire accurate averages
	for (j=0;j<DIFFTESTS;j++){
		avg[j] = test[j](0);
	}

	//Print stats
	for (j=0;j<DIFFTESTS;j++){
		printf("%s:\t%.2f\tCAvg,\t%.2f\tC/Vec\n",test_descriptions[j],avg[j]/tests,avg[j]/(tests*num_vectors[j	]));
	}
}






uint64_t test_mult_4_shuffle(int zero){

	__m128i x;
	uint64_t t;
	int i;

	t = __rdtsc();
	for(i=0;i<tests;i++){
		x = anubis_mult_4_shuffle(x);
	}
	t = __rdtsc() - t;

	use_m128i_variable(x,zero);

	return t;
}

uint64_t test_mult_4_shuffle_macro(int zero){

	__m128i res,tmp1,tmp2,x;
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_4_shuffle_macro(x,tmp1,tmp2,res);
		anubis_mult_4_shuffle_macro(res,tmp1,tmp2,x);
	}

	t = __rdtsc() - t;

	use_m128i_variable(x,zero);

	return t;
}


uint64_t test_mult_2_shuffle_gf16(int zero){

	__m128i x;
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests;i++){
		x = gf16_anubis_mult_2_shuffle(x);
	}

	t = __rdtsc() - t;

	use_m128i_variable(x,zero);

	return t;
}

uint64_t test_mult_2_shuffle_gf16_macro(int zero){

	__m128i res,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,x;
	uint64_t t;
	int i;

	//time cycles
	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		gf16_anubis_mult_2_shuffle_macro(x,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,res);
		gf16_anubis_mult_2_shuffle_macro(res,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,x);
	}

	t = __rdtsc() - t;

	use_m128i_variable(x,zero);

	return t;
}



uint64_t test_mult_16_shuffle_anubis(int zero){
	__m128i input[4],output[4];
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_16_shuffle(input,output);
		anubis_mult_16_shuffle(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}

uint64_t test_mult_16_shuffle_other(int zero){
	__m128i input[4],output[4];
	uint64_t t;
	int i;


	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		h2_mult_16_shuffle(input,output);
		h2_mult_16_shuffle(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}

uint64_t test_mult_16_shuffle_load_anubis(int zero){
	__m128i input[4],output[4];
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_16_shuffle_mem(input,output);
		anubis_mult_16_shuffle_mem(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}

uint64_t test_mult_16_shuffle_load_other(int zero){
	__m128i input[4],output[4];
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		h2_mult_16_shuffle_mem(input,output);
		h2_mult_16_shuffle_mem(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}


uint64_t test_mult_16_x2_anubis(int zero){
	__m128i ret[4],v[4];
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_16_x2(v,ret);
		anubis_mult_16_x2(ret,v);
	}

	t = __rdtsc() - t;

	use_m128i_variables(v,zero);

	return t;
}

uint64_t test_mult_16_x2_other(int zero){
	__m128i ret[4],v[4];
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		h2_mult_16_x2(v,ret);
		h2_mult_16_x2(ret,v);
	}

	t = __rdtsc() - t;

	use_m128i_variables(v,zero);

	return t;
}

uint64_t test_mult_16_x2_load_anubis(int zero){
	__m128i input[4],output[4];
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_16_x2_mem(input,output);
		anubis_mult_16_x2_mem(output,input);
	}
		
	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}

uint64_t test_mult_16_x2_load_other(int zero){
	__m128i input[4],output[4];
	uint64_t t;
	int i;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		h2_mult_16_x2_mem(input,output);
		h2_mult_16_x2_mem(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}


uint64_t test_mult_128_xor_anubis(int zero){
	__m128i input[32], output[32];
	int i;
	uint64_t t;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_128_xor(input,output);
		anubis_mult_128_xor(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}

uint64_t test_mult_128_xor_anubis_opt_compiler(int zero){
	__m128i input[32], output[32];
	int i;
	uint64_t t;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_128_xor_opt_compiler(input,output);
		anubis_mult_128_xor_opt_compiler(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}

uint64_t test_mult_128_xor_anubis_opt_tmp(int zero){
	__m128i input[32], output[32];
	int i;
	uint64_t t;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_128_xor_opt_tmp(input,output);
		anubis_mult_128_xor_opt_tmp(output,input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}

uint64_t test_mult_128_xor_anubis_macro(int zero){
	__m128i m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31,m128i_tmp0,m128i_tmp1;
	int i;
	uint64_t t;

	t = __rdtsc();

	for(i=0;i<tests/2;i++){
		anubis_mult_128_xor_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31,m128i_tmp0,m128i_tmp1);
		anubis_mult_128_xor_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31,m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_tmp0,m128i_tmp1);
	}

	t = __rdtsc() - t;

	use_m128i_variable(m128i_input0,zero);
	use_m128i_variable(m128i_input1,zero);
	use_m128i_variable(m128i_input2,zero);
	use_m128i_variable(m128i_input3,zero);
	use_m128i_variable(m128i_input4,zero);
	use_m128i_variable(m128i_input5,zero);
	use_m128i_variable(m128i_input6,zero);
	use_m128i_variable(m128i_input7,zero);
	use_m128i_variable(m128i_input8,zero);
	use_m128i_variable(m128i_input9,zero);
	use_m128i_variable(m128i_input10,zero);
	use_m128i_variable(m128i_input11,zero);
	use_m128i_variable(m128i_input12,zero);
	use_m128i_variable(m128i_input13,zero);
	use_m128i_variable(m128i_input14,zero);
	use_m128i_variable(m128i_input15,zero);
	use_m128i_variable(m128i_input16,zero);
	use_m128i_variable(m128i_input17,zero);
	use_m128i_variable(m128i_input18,zero);
	use_m128i_variable(m128i_input19,zero);
	use_m128i_variable(m128i_input20,zero);
	use_m128i_variable(m128i_input21,zero);
	use_m128i_variable(m128i_input22,zero);
	use_m128i_variable(m128i_input23,zero);
	use_m128i_variable(m128i_input24,zero);
	use_m128i_variable(m128i_input25,zero);
	use_m128i_variable(m128i_input26,zero);
	use_m128i_variable(m128i_input27,zero);
	use_m128i_variable(m128i_input28,zero);
	use_m128i_variable(m128i_input29,zero);
	use_m128i_variable(m128i_input30,zero);
	use_m128i_variable(m128i_input31,zero);

	return t;
}

uint64_t test_mult_128_xor_load_anubis(int zero){
	__m128i input[32];
	int i;
	uint64_t t;

	t = __rdtsc();

	for(i=0;i<tests;i++){
		anubis_mult_128_xor_mem(input);
	}

	t = __rdtsc() - t;

	use_m128i_variables(input,zero);

	return t;
}