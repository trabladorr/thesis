#ifndef _MULT_UTILS_INLINE
#define _MULT_UTILS_INLINE


//transpose 4x4, 32bit value matrix in input(4x128)
static inline __m128i transpose_4x4_8bit(__m128i input){
	const __m128i transpose_4x4_8bit_mask = _mm_setr_epi8(0x0, 0x4, 0x8, 0xc, 0x1, 0x5, 0x9, 0xd, 0x2, 0x6, 0xa, 0xe, 0x3, 0x7, 0xb, 0xf);
	return _mm_shuffle_epi8(input,transpose_4x4_8bit_mask);
}

//transpose 4x4, 32bit value matrix in input(4x128)
static inline void transpose_4x4_32bit(__m128i input[]){
	__m128i tmp0 = _mm_unpacklo_epi32(input[0], input[1]);
	__m128i tmp1 = _mm_unpacklo_epi32(input[2], input[3]);
	__m128i tmp2 = _mm_unpackhi_epi32(input[0], input[1]);
	__m128i tmp3 = _mm_unpackhi_epi32(input[2], input[3]);

	input[0] = _mm_unpacklo_epi64(tmp0, tmp1);
	input[1] = _mm_unpackhi_epi64(tmp0, tmp1);
	input[2] = _mm_unpacklo_epi64(tmp2, tmp3);
	input[3] = _mm_unpackhi_epi64(tmp2, tmp3);
}

static inline void into_byteslice_4x128(__m128i input[]){

	//transpose the 4x4, 8bit value matrices
	input[0] = transpose_4x4_8bit(input[0]);
	input[1] = transpose_4x4_8bit(input[1]);
	input[2] = transpose_4x4_8bit(input[2]);
	input[3] = transpose_4x4_8bit(input[3]);

	transpose_4x4_32bit(input);
}

static inline void from_byteslice_4x128(__m128i input[]){

	transpose_4x4_32bit(input);
	
	//transpose the 4x4, 8bit value matrices
	input[0] = transpose_4x4_8bit(input[0]);
	input[1] = transpose_4x4_8bit(input[1]);
	input[2] = transpose_4x4_8bit(input[2]);
	input[3] = transpose_4x4_8bit(input[3]);

}
#endif

