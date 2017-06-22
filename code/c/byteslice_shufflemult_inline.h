#ifndef _MULT_BYTESLICE_SHUFFLEMULT_INLINE
#define _MULT_BYTESLICE_SHUFFLEMULT_INLINE

//calculate 16 multiplications of 8-bit values, using two shuffle multiplications
static inline __m128i shuffle_mult(__m128i x_low,__m128i x_high,__m128i in,__m128i in_high){
	__m128i  tmp = _mm_shuffle_epi8(x_low,in);
	return _mm_xor_si128(tmp, _mm_shuffle_epi8(x_high,in_high));
}

//input v contains 16 vectors to mutliply, bytesliced, output ret contains 16 vectors mutliplied, bytesliced
//multiply using two shuffle multiplications, anubis matrix/field
static inline void anubis_mult_16_shuffle(__m128i v[],__m128i ret[]){

	const __m128i midword_zeroes = _mm_setr_epi8(0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0);
	const __m128i shuffle_no_zeroing_mask = _mm_setr_epi8(0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f);

	const __m128i mul_x2_low = _mm_setr_epi8( 0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x10, 0x12, 0x14, 0x16, 0x18, 0x1a, 0x1c, 0x1e );
	const __m128i mul_x4_low = _mm_setr_epi8( 0x0, 0x4, 0x8, 0xc, 0x10, 0x14, 0x18, 0x1c, 0x20, 0x24, 0x28, 0x2c, 0x30, 0x34, 0x38, 0x3c );
	const __m128i mul_x6_low = _mm_setr_epi8( 0x0, 0x6, 0xc, 0xa, 0x18, 0x1e, 0x14, 0x12, 0x30, 0x36, 0x3c, 0x3a, 0x28, 0x2e, 0x24, 0x22 );
	const __m128i mul_x2_high = _mm_setr_epi8( 0x0, 0x20, 0x40, 0x60, 0x80, 0xa0, 0xc0, 0xe0, 0x1d, 0x3d, 0x5d, 0x7d, 0x9d, 0xbd, 0xdd, 0xfd );
	const __m128i mul_x4_high = _mm_setr_epi8( 0x0, 0x40, 0x80, 0xc0, 0x1d, 0x5d, 0x9d, 0xdd, 0x3a, 0x7a, 0xba, 0xfa, 0x27, 0x67, 0xa7, 0xe7 );
	const __m128i mul_x6_high = _mm_setr_epi8( 0x0, 0x60, 0xc0, 0xa0, 0x9d, 0xfd, 0x5d, 0x3d, 0x27, 0x47, 0xe7, 0x87, 0xba, 0xda, 0x7a, 0x1a );

	__m128i v_high[4];

	//diagonal of ones
	ret[0] = v[0];
	ret[1] = v[1];
	ret[2] = v[2];
	ret[3] = v[3];

	//high registers contain the low register bytes right-shifted by 4
	v_high[0] = _mm_srli_epi16(_mm_and_si128(v[0], midword_zeroes), 4);
	v_high[1] = _mm_srli_epi16(_mm_and_si128(v[1], midword_zeroes), 4);
	v_high[2] = _mm_srli_epi16(_mm_and_si128(v[2], midword_zeroes), 4);
	v_high[3] = _mm_srli_epi16(_mm_and_si128(v[3], midword_zeroes), 4);

	//mask to prevent shuffle zeroing based on highest order bit, ORed to low order multiplication
	v[0] = _mm_and_si128(v[0],shuffle_no_zeroing_mask);
	v[1] = _mm_and_si128(v[1],shuffle_no_zeroing_mask);
	v[2] = _mm_and_si128(v[2],shuffle_no_zeroing_mask);
	v[3] = _mm_and_si128(v[3],shuffle_no_zeroing_mask);
	

	//multiply by shuffling, xor for results
	ret[0] = _mm_xor_si128 (ret[0], shuffle_mult(mul_x2_low,mul_x2_high,v[1],v_high[1]));
	ret[0] = _mm_xor_si128 (ret[0], shuffle_mult(mul_x4_low,mul_x4_high,v[2],v_high[2]));
	ret[0] = _mm_xor_si128 (ret[0], shuffle_mult(mul_x6_low,mul_x6_high,v[3],v_high[3]));

	ret[1] = _mm_xor_si128 (ret[1], shuffle_mult(mul_x2_low,mul_x2_high,v[0],v_high[0]));
	ret[1] = _mm_xor_si128 (ret[1], shuffle_mult(mul_x6_low,mul_x6_high,v[2],v_high[2]));
	ret[1] = _mm_xor_si128 (ret[1], shuffle_mult(mul_x4_low,mul_x4_high,v[3],v_high[3]));

	ret[2] = _mm_xor_si128 (ret[2], shuffle_mult(mul_x4_low,mul_x4_high,v[0],v_high[0]));
	ret[2] = _mm_xor_si128 (ret[2], shuffle_mult(mul_x6_low,mul_x6_high,v[1],v_high[1]));
	ret[2] = _mm_xor_si128 (ret[2], shuffle_mult(mul_x2_low,mul_x2_high,v[3],v_high[3]));

	ret[3] = _mm_xor_si128 (ret[3], shuffle_mult(mul_x6_low,mul_x6_high,v[0],v_high[0]));
	ret[3] = _mm_xor_si128 (ret[3], shuffle_mult(mul_x4_low,mul_x4_high,v[1],v_high[1]));
	ret[3] = _mm_xor_si128 (ret[3], shuffle_mult(mul_x2_low,mul_x2_high,v[2],v_high[2]));
}

//input v contains 16 vectors to mutliply, bytesliced, output ret contains 16 vectors mutliplied, bytesliced
//multiply using two shuffle multiplications, other Hadamard matrix/field
static inline void h2_mult_16_shuffle(__m128i v[],__m128i ret[]){

	const __m128i midword_zeroes = _mm_setr_epi8(0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0);
	const __m128i shuffle_no_zeroing_mask = _mm_setr_epi8(0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f);

	const __m128i mul_x2_low = _mm_setr_epi8( 0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x10, 0x12, 0x14, 0x16, 0x18, 0x1a, 0x1c, 0x1e );
	const __m128i mul_xb0_low = _mm_setr_epi8( 0x0, 0xb0, 0x5, 0xb5, 0xa, 0xba, 0xf, 0xbf, 0x14, 0xa4, 0x11, 0xa1, 0x1e, 0xae, 0x1b, 0xab );
	const __m128i mul_xb2_low = _mm_setr_epi8( 0x0, 0xb2, 0x1, 0xb3, 0x2, 0xb0, 0x3, 0xb1, 0x4, 0xb6, 0x5, 0xb7, 0x6, 0xb4, 0x7, 0xb5 );
	
	const __m128i mul_x2_high = _mm_setr_epi8( 0x0, 0x20, 0x40, 0x60, 0x80, 0xa0, 0xc0, 0xe0, 0x65, 0x45, 0x25, 0x5, 0xe5, 0xc5, 0xa5, 0x85 );
	const __m128i mul_xb0_high = _mm_setr_epi8( 0x0, 0x28, 0x50, 0x78, 0xa0, 0x88, 0xf0, 0xd8, 0x25, 0xd, 0x75, 0x5d, 0x85, 0xad, 0xd5, 0xfd );
	const __m128i mul_xb2_high = _mm_setr_epi8( 0x0, 0x8, 0x10, 0x18, 0x20, 0x28, 0x30, 0x38, 0x40, 0x48, 0x50, 0x58, 0x60, 0x68, 0x70, 0x78 );
	
	__m128i v_high[4];

	//diagonal of ones
	ret[0] = v[0];
	ret[1] = v[1];
	ret[2] = v[2];
	ret[3] = v[3];

	//high registers contain the low register bytes right-shifted by 4
	v_high[0] = _mm_srli_epi16(_mm_and_si128(v[0], midword_zeroes), 4);
	v_high[1] = _mm_srli_epi16(_mm_and_si128(v[1], midword_zeroes), 4);
	v_high[2] = _mm_srli_epi16(_mm_and_si128(v[2], midword_zeroes), 4);
	v_high[3] = _mm_srli_epi16(_mm_and_si128(v[3], midword_zeroes), 4);

	//mask to prevent shuffle zeroing based on highest order bit, ORed to low order multiplication
	v[0] = _mm_and_si128(v[0],shuffle_no_zeroing_mask);
	v[1] = _mm_and_si128(v[1],shuffle_no_zeroing_mask);
	v[2] = _mm_and_si128(v[2],shuffle_no_zeroing_mask);
	v[3] = _mm_and_si128(v[3],shuffle_no_zeroing_mask);

	//multiply by shuffling, xor for results
	ret[0] = _mm_xor_si128 (ret[0], shuffle_mult(mul_x2_low,mul_x2_high,v[1],v_high[1]));
	ret[0] = _mm_xor_si128 (ret[0], shuffle_mult(mul_xb0_low,mul_xb0_high,v[2],v_high[2]));
	ret[0] = _mm_xor_si128 (ret[0], shuffle_mult(mul_xb2_low,mul_xb2_high,v[3],v_high[3]));

	ret[1] = _mm_xor_si128 (ret[1], shuffle_mult(mul_x2_low,mul_x2_high,v[0],v_high[0]));
	ret[1] = _mm_xor_si128 (ret[1], shuffle_mult(mul_xb2_low,mul_xb2_high,v[2],v_high[2]));
	ret[1] = _mm_xor_si128 (ret[1], shuffle_mult(mul_xb0_low,mul_xb0_high,v[3],v_high[3]));

	ret[2] = _mm_xor_si128 (ret[2], shuffle_mult(mul_xb0_low,mul_xb0_high,v[0],v_high[0]));
	ret[2] = _mm_xor_si128 (ret[2], shuffle_mult(mul_xb2_low,mul_xb2_high,v[1],v_high[1]));
	ret[2] = _mm_xor_si128 (ret[2], shuffle_mult(mul_x2_low,mul_x2_high,v[3],v_high[3]));

	ret[3] = _mm_xor_si128 (ret[3], shuffle_mult(mul_xb2_low,mul_xb2_high,v[0],v_high[0]));
	ret[3] = _mm_xor_si128 (ret[3], shuffle_mult(mul_xb0_low,mul_xb0_high,v[1],v_high[1]));
	ret[3] = _mm_xor_si128 (ret[3], shuffle_mult(mul_x2_low,mul_x2_high,v[2],v_high[2]));
}

//input/output pointer contains 16 vectors to multiply
//multiply using two shuffle multiplications, anubis matrix/field
static inline void anubis_mult_16_shuffle_mem(__m128i input[],__m128i output[]){

	into_byteslice_4x128(input);
	
	anubis_mult_16_shuffle(input,output);

	from_byteslice_4x128(output);
}


//input/output pointer contains 16 vectors to multiply
//multiply using two shuffle multiplications, other Hadamard matrix/field
static inline void h2_mult_16_shuffle_mem(__m128i input[],__m128i output[]){

	into_byteslice_4x128(input);
	
	h2_mult_16_shuffle(input,output);

	from_byteslice_4x128(output);
}


 #endif