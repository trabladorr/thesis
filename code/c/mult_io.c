#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "x86intrin.h"

#include "mult_utils.h"
#include "mult_io.h"
#include "utils_inline.h"
#include "utils_macro.h"

#include "direct_shufflemult_inline.h"
#include "direct_shufflemult_macro.h"

#include "byteslice_shufflemult_inline.h"
#include "byteslice_x2mult_inline.h"

#include "bitslice_xormult_inline.h"
#include "bitslice_xormult_macro.h"

#define DEBUG_VERBOSE 0

int main(int argc, char **argv){
	parameters params = init(argc,argv);
	init_macro_constants();
	

	int i = 0;
	for (;i<params.n_tests;++i){
		mult_io(params,DEBUG_VERBOSE);
	}
	
	return 0;
}

//parameters choose implementation, stdin accepts vectors as raw bytes, stdout returns multiplication results
//allows testing implementation for correctness
void mult_io(parameters params, int verbose){

	int no_function_flag = 0;
	int ret = 0;
	ssize_t input_size = params.n_vectors*4*params.f/8;

	void *buffer = malloc(input_size);
	if (!buffer){
		perror("mult_io:malloc");
		exit(EXIT_FAILURE);
	}
	

	ret = read(0, buffer, input_size);
	if (ret != input_size){
		if (ret < 0)
			perror("mult_io:read");
		else 
			fprintf(stderr, "mult_io:read: Read %d bytes instead of %ld\n",ret,input_size);
		exit(EXIT_FAILURE);
	}
	if(verbose){
		printf("Read %d bytes:\n",ret);
		printNxMbytes(buffer,input_size/16,16);
	}

	//ugly jungle of if else cases, essentially each case is choice of implementation
	//each case loads input from stdin, performs multiplication, writes output to stdout
	if (params.n_vectors == 2){
		if (params.table != TABLE_ANUBIS)
			no_function_flag = 1;
		else if (params.mult_type == MULT_SHUFFLE){
			if (params.impl_type == IMPL_INLINE){
				__m128i x;
				memcpy(&x,buffer,input_size);

				if (verbose){
					printf("Vec-Mat mult using shuffle, GF(16), 2 input vectors(concatenated):\n");
					printNxMbytes(&x,1,16);
				}

				x = gf16_anubis_mult_2_shuffle(x);

				if (verbose){
					printf("Result:\n");
					printNxMbytes(&x,1,16);
				}

				memcpy(buffer,&x,input_size);
			}
			else if (params.impl_type == IMPL_MACRO){
				__m128i res,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,x;

				memcpy(&x,buffer,input_size);

				if (verbose){
					printf("Vec-Mat mult using shuffle macro, GF(16), 2 input vectors(concatenated):\n");
					printNxMbytes(&x,1,16);
				}

				gf16_anubis_mult_2_shuffle_macro(x,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,res);
				
				if (verbose){
					printf("Result:\n");
					printNxMbytes(&res,1,16);
				}

				memcpy(buffer,&res,input_size);
			}
			else
				no_function_flag = 1;
		}
		else
			no_function_flag = 1;
	}
	else if (params.n_vectors == 4){
		if (params.table != TABLE_ANUBIS)
			no_function_flag = 1;
		else if (params.mult_type == MULT_SHUFFLE){
			if (params.impl_type == IMPL_INLINE){
				__m128i x;
				memcpy(&x,buffer,input_size);

				if (verbose){
					printf("Vec-Mat mult using shuffle, Anubis, 4 input vectors(concatenated):\n");
					printNxMbytes(&x,1,16);
				}

				x = anubis_mult_4_shuffle(x);

				if (verbose){
					printf("Result:\n");
					printNxMbytes(&x,1,16);
				}

				memcpy(buffer,&x,input_size);
			}
			else if (params.impl_type == IMPL_MACRO){
				__m128i res,tmp1,tmp2,x;

				memcpy(&x,buffer,input_size);

				if (verbose){
					printf("Vec-Mat mult using shuffle macro, Anubis, 4 input vectors(concatenated):\n");
					printNxMbytes(&x,1,16);
				}

				anubis_mult_4_shuffle_macro(x,tmp1,tmp2,res);

				if (verbose){
					printf("Result:\n");
					printNxMbytes(&res,1,16);
				}

				memcpy(buffer,&res,input_size);
			}
			else
				no_function_flag = 1;
		}
		else
			no_function_flag = 1;
	}
	else if (params.n_vectors == 16){
		if (params.table == TABLE_ANUBIS){
			if (params.mult_type == MULT_SHUFFLE){
				if (params.impl_type == IMPL_INLINE){
					__m128i input[4],output[4];
					memcpy(input,buffer,input_size);

					if (verbose){
						printf("Vec-Mat mult using shuffle, Anubis, 16 input vectors, with loading:\n");
						printNxMbytes(input,4,16);
					}

					anubis_mult_16_shuffle_mem(input,output);

					if (verbose){
						printf("Result:\n");
						printNxMbytes(output,4,16);
					}

					memcpy(buffer,output,input_size);
				}
				else
					no_function_flag = 1;
			}
			else if (params.mult_type == MULT_X2){
				if (params.impl_type == IMPL_INLINE){
					__m128i input[4],output[4];
					memcpy(input,buffer,input_size);

					if (verbose){
						printf("Vec-Mat mult using  x2, Anubis, 16 input vectors, with loading:\n");
						printNxMbytes(input,4,16);
					}

					anubis_mult_16_x2_mem(input,output);
						
					if (verbose){
						printf("Result:\n");
						printNxMbytes(output,4,16);
					}

					memcpy(buffer,output,input_size);
				}
				else
					no_function_flag = 1;
			}
			else
				no_function_flag = 1;
		}
		else if (params.table == TABLE_OTHERH){
			if (params.mult_type == MULT_SHUFFLE){
				if (params.impl_type == IMPL_INLINE){
					__m128i input[4],output[4];
					memcpy(input,buffer,input_size);

					if (verbose){
						printf("Vec-Mat mult using shuffle, Other Hadamard, 16 input vectors, with loading:\n");
						printNxMbytes(input,4,16);
					}

					h2_mult_16_shuffle_mem(input,output);

					if (verbose){
						printf("Result:\n");
						printNxMbytes(output,4,16);
					}

					memcpy(buffer,output,input_size);
				}
				else
					no_function_flag = 1;
			}
			else if (params.mult_type == MULT_X2){
				if (params.impl_type == IMPL_INLINE){
					__m128i input[4],output[4];
					memcpy(input,buffer,input_size);

					if (verbose){
						printf("Vec-Mat mult using  x2, Other Hadamard, 16 input vectors, with loading:\n");
						printNxMbytes(input,4,16);
					}

					h2_mult_16_x2_mem(input,output);
						
					if (verbose){
						printf("Result:\n");
						printNxMbytes(output,4,16);
					}

					memcpy(buffer,output,input_size);
				}
				else
					no_function_flag = 1;
			}
			else
				no_function_flag = 1;
		}
		else
			no_function_flag = 1;
	}
	else if (params.n_vectors == 128){
		if (params.table != TABLE_ANUBIS)
			no_function_flag = 1;
		else if (params.mult_type == MULT_XOR){
			if (params.impl_type == IMPL_INLINE){
				__m128i input[32];
				memcpy(input,buffer,input_size);

				if (verbose){
					printf("Vec-Mat mult using x2, Anubis, 128 input vectors, (bitsliced):\n");
					printNxMbytes(input,32,16);
				}

				anubis_mult_128_xor_mem(input);

				if (verbose){
					printf("Result:\n");
					printNxMbytes(input,32,16);
				}

				memcpy(buffer,input,input_size);
			}
			else if (params.impl_type == IMPL_INLINE_OPT){
				__m128i input[32], tmp[32];
				memcpy(input,buffer,input_size);

				if (verbose){
					printf("Vec-Mat mult using x2, Anubis, 128 input vectors, opt Compiler, (bitsliced):\n");
					printNxMbytes(input,32,16);
				}


				into_bitslice(input,tmp);
				anubis_mult_128_xor_opt_compiler(input,tmp);
				from_bitslice(input,tmp);

				if (verbose){
					printf("Result:\n");
					printNxMbytes(input,32,16);
				}

				memcpy(buffer,input,input_size);
			}
			else if (params.impl_type == IMPL_INLINE_OPT2){
				__m128i input[32], tmp[32];
				memcpy(input,buffer,input_size);

				if (verbose){
					printf("Vec-Mat mult using x2, Anubis, 128 input vectors, opt Tmp, (bitsliced):\n");
					printNxMbytes(input,32,16);
				}

				into_bitslice(input,tmp);
				anubis_mult_128_xor_opt_tmp(input,tmp);
				from_bitslice(input,tmp);

				if (verbose){
					printf("Result:\n");
					printNxMbytes(input,32,16);
				}

				memcpy(buffer,input,input_size);
			}
			else
				no_function_flag = 1;
		}
		else
			no_function_flag = 1;
	}
	else
		no_function_flag = 1;

	if (no_function_flag){
		fprintf(stderr, "No function found for f=%d, n_vectors=%d, table=%d, mult_type=%d, impl_type=%d\n", params.f, params.n_vectors, params.table, params.mult_type, params.impl_type);
	}


	ret = write(1, buffer, input_size);
	fflush(stdout);

	if (ret != input_size){
		if (ret < 0)
			perror("mult_io:write");
		else 
			fprintf(stderr, "mult_io:write: Wrote %d bytes instead of %ld\n",ret,input_size);
		exit(EXIT_FAILURE);
	}
	
	
	if(verbose){
		printf("\nWrote %d bytes:\n",ret);
		printNxMbytes(buffer,input_size/16,16);
	}

	free(buffer);
}

parameters init(int argc, char **argv){
	if (argc != 7){
		printf("Usage: [f] [n_vectors] [n_tests] [%d:TABLE_ANUBIS / %d:TABLE_OTHERH] [%d:MULT_SHUFFLE / %d:MULT_X2 / %d:MULT_XOR] [%d:IMPL_INLINE / %d:IMPL_MACRO]\n",TABLE_ANUBIS,TABLE_OTHERH,MULT_SHUFFLE,MULT_X2,MULT_XOR,IMPL_INLINE,IMPL_MACRO);
		exit(EXIT_FAILURE);
	}

	parameters params;
	params.f = atoi(argv[1]);
	params.n_vectors = atoi(argv[2]);
	params.n_tests = atoi(argv[3]);
	params.table = atoi(argv[4]);
	params.mult_type = atoi(argv[5]);
	params.impl_type = atoi(argv[6]);

	if (params.f < 1){
		fprintf(stderr, "Invalid field exponent:%d\n",params.f);
		exit(EXIT_FAILURE);
	}
	else if (params.n_vectors < 1){
		fprintf(stderr, "Invalid number of vectors:%d\n",params.n_vectors);
		exit(EXIT_FAILURE);
	}
	else if (params.n_tests < 1){
		fprintf(stderr, "Invalid number of tests:%d\n",params.n_tests);
		exit(EXIT_FAILURE);
	}


	return params;
}