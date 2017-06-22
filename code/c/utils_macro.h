#ifndef _MULT_UTILS_MACRO
#define _MULT_UTILS_MACRO

__m128i rearrange_odds_evens_mask;
__m128i transpose_4x4_8bit_mask;

__m128i midword_zeroes;
__m128i shuffle_no_zeroing_mask;
__m128i shuffle_zeroing_odd_mask;
__m128i shuffle_zeroing_even_mask;

__m128i shuffle_x2_gf8;
__m128i shuffle_x4_gf8;
__m128i shuffle_x6_gf8;

__m128i mul_x2_low;
__m128i mul_x4_low;
__m128i mul_x6_low;

__m128i mul_x2_high;
__m128i mul_x4_high;
__m128i mul_x6_high;

__m128i shuffle_x2_gf16;
__m128i shuffle_x4_gf16;
__m128i shuffle_x6_gf16;

__m128i mul_x2_lower_cdef;
__m128i mul_x2_lower_89ab;
__m128i mul_x2_lower_0123;
__m128i mul_x2_upper_89ab;
__m128i mul_x2_upper_4567;
__m128i mul_x2_upper_0123;

__m128i mul_x4_lower_cdef;
__m128i mul_x4_lower_89ab;
__m128i mul_x4_lower_0123;
__m128i mul_x4_upper_89ab;
__m128i mul_x4_upper_4567;
__m128i mul_x4_upper_0123;

__m128i mul_x6_lower_cdef;
__m128i mul_x6_lower_89ab;
__m128i mul_x6_lower_0123;
__m128i mul_x6_upper_89ab;
__m128i mul_x6_upper_4567;
__m128i mul_x6_upper_0123;

#define init_macro_constants() do {\
	\
	rearrange_odds_evens_mask = _mm_setr_epi8(0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x1, 0x3, 0x5, 0x7, 0x9, 0xb, 0xd, 0xf);\
	transpose_4x4_8bit_mask = _mm_setr_epi8(0x0, 0x4, 0x8, 0xc, 0x1, 0x5, 0x9, 0xd, 0x2, 0x6, 0xa, 0xe, 0x3, 0x7, 0xb, 0xf);\
	midword_zeroes = _mm_setr_epi8(0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0);\
	shuffle_no_zeroing_mask = _mm_setr_epi8(0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f);\
	shuffle_zeroing_odd_mask = _mm_setr_epi8(0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80);\
	shuffle_zeroing_even_mask = _mm_setr_epi8(0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff,0x80,0xff);\
	\
	shuffle_x2_gf8 = _mm_setr_epi8(0x1,0x0,0x3,0x2,0x5,0x4,0x7,0x6,0x9,0x8,0xb,0xa,0xd,0xc,0xf,0xe);\
	shuffle_x4_gf8 = _mm_setr_epi8(0x2,0x3,0x0,0x1,0x6,0x7,0x4,0x5,0xa,0xb,0x8,0x9,0xe,0xf,0xc,0xd);\
	shuffle_x6_gf8 = _mm_setr_epi8(0x3,0x2,0x1,0x0,0x7,0x6,0x5,0x4,0xb,0xa,0x9,0x8,0xf,0xe,0xd,0xc);\
	\
	mul_x2_low = _mm_setr_epi8( 0x0,0x2,0x4,0x6,0x8,0xa,0xc,0xe,0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e );\
	mul_x4_low = _mm_setr_epi8( 0x0,0x4,0x8,0xc,0x10,0x14,0x18,0x1c,0x20,0x24,0x28,0x2c,0x30,0x34,0x38,0x3c );\
	mul_x6_low = _mm_setr_epi8( 0x0,0x6,0xc,0xa,0x18,0x1e,0x14,0x12,0x30,0x36,0x3c,0x3a,0x28,0x2e,0x24,0x22 );\
	\
	mul_x2_high = _mm_setr_epi8( 0x0,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0,0x1d,0x3d,0x5d,0x7d,0x9d,0xbd,0xdd,0xfd );\
	mul_x4_high = _mm_setr_epi8( 0x0,0x40,0x80,0xc0,0x1d,0x5d,0x9d,0xdd,0x3a,0x7a,0xba,0xfa,0x27,0x67,0xa7,0xe7 );\
	mul_x6_high = _mm_setr_epi8( 0x0,0x60,0xc0,0xa0,0x9d,0xfd,0x5d,0x3d,0x27,0x47,0xe7,0x87,0xba,0xda,0x7a,0x1a );\
	\
	shuffle_x2_gf16 = _mm_setr_epi8(0x2,0x3,0x0,0x1,0x6,0x7,0x4,0x5,0xa,0xb,0x8,0x9,0xe,0xf,0xc,0xd);\
	shuffle_x4_gf16 = _mm_setr_epi8(0x4,0x5,0x6,0x7,0x0,0x1,0x2,0x3,0xc,0xd,0xe,0xf,0x8,0x9,0xa,0xb);\
	shuffle_x6_gf16 = _mm_setr_epi8(0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1,0xe,0xf,0xc,0xd,0xa,0xb,0x8,0x9);\
	\
	mul_x2_lower_cdef = _mm_setr_epi8( 0x0,0x2,0x4,0x6,0x8,0xa,0xc,0xe,0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e );\
	mul_x2_lower_89ab = _mm_setr_epi8( 0x0,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0,0x0,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0 );\
	mul_x2_lower_0123 = _mm_setr_epi8( 0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x2d,0x2d,0x2d,0x2d,0x2d,0x2d,0x2d,0x2d );\
	mul_x2_upper_89ab = _mm_setr_epi8( 0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x1,0x1,0x1,0x1,0x1,0x1,0x1,0x1 );\
	mul_x2_upper_4567 = _mm_setr_epi8( 0x0,0x2,0x4,0x6,0x8,0xa,0xc,0xe,0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e );\
	mul_x2_upper_0123 = _mm_setr_epi8( 0x0,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0,0x0,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0 );\
	\
	mul_x4_lower_cdef = _mm_setr_epi8( 0x0,0x4,0x8,0xc,0x10,0x14,0x18,0x1c,0x20,0x24,0x28,0x2c,0x30,0x34,0x38,0x3c );\
	mul_x4_lower_89ab = _mm_setr_epi8( 0x0,0x40,0x80,0xc0,0x0,0x40,0x80,0xc0,0x0,0x40,0x80,0xc0,0x0,0x40,0x80,0xc0 );\
	mul_x4_lower_0123 = _mm_setr_epi8( 0x0,0x0,0x0,0x0,0x2d,0x2d,0x2d,0x2d,0x5a,0x5a,0x5a,0x5a,0x77,0x77,0x77,0x77 );\
	mul_x4_upper_89ab = _mm_setr_epi8( 0x0,0x0,0x0,0x0,0x1,0x1,0x1,0x1,0x2,0x2,0x2,0x2,0x3,0x3,0x3,0x3 );\
	mul_x4_upper_4567 = _mm_setr_epi8( 0x0,0x4,0x8,0xc,0x10,0x14,0x18,0x1c,0x20,0x24,0x28,0x2c,0x30,0x34,0x38,0x3c );\
	mul_x4_upper_0123 = _mm_setr_epi8( 0x0,0x40,0x80,0xc0,0x0,0x40,0x80,0xc0,0x0,0x40,0x80,0xc0,0x0,0x40,0x80,0xc0 );\
	\
	mul_x6_lower_cdef = _mm_setr_epi8( 0x0,0x6,0xc,0xa,0x18,0x1e,0x14,0x12,0x30,0x36,0x3c,0x3a,0x28,0x2e,0x24,0x22 );\
	mul_x6_lower_89ab = _mm_setr_epi8( 0x0,0x60,0xc0,0xa0,0x80,0xe0,0x40,0x20,0x0,0x60,0xc0,0xa0,0x80,0xe0,0x40,0x20 );\
	mul_x6_lower_0123 = _mm_setr_epi8( 0x0,0x0,0x0,0x0,0x2d,0x2d,0x2d,0x2d,0x77,0x77,0x77,0x77,0x5a,0x5a,0x5a,0x5a );\
	mul_x6_upper_89ab = _mm_setr_epi8( 0x0,0x0,0x0,0x0,0x1,0x1,0x1,0x1,0x3,0x3,0x3,0x3,0x2,0x2,0x2,0x2 );\
	mul_x6_upper_4567 = _mm_setr_epi8( 0x0,0x6,0xc,0xa,0x18,0x1e,0x14,0x12,0x30,0x36,0x3c,0x3a,0x28,0x2e,0x24,0x22 );\
	mul_x6_upper_0123 = _mm_setr_epi8( 0x0,0x60,0xc0,0xa0,0x80,0xe0,0x40,0x20,0x0,0x60,0xc0,0xa0,0x80,0xe0,0x40,0x20 );\
} while(0);

//set all 8 128 output registers to input
#define set_8x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7) do {\
	m128i_output0 = m128i_input0;\
	m128i_output1 = m128i_input1;\
	m128i_output2 = m128i_input2;\
	m128i_output3 = m128i_input3;\
	m128i_output4 = m128i_input4;\
	m128i_output5 = m128i_input5;\
	m128i_output6 = m128i_input6;\
	m128i_output7 = m128i_input7;\
} while(0);

//set all 32 128 output registers to input
#define set_32x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31) do {\
	set_8x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7);\
	set_8x128_macro(m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15);\
	set_8x128_macro(m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23);\
	set_8x128_macro(m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31);\
} while(0);

// //output = output XOR input
// //8x128 bit bitsliced vectors
#define xor_8x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7) do {\
	m128i_output0 = _mm_xor_si128(m128i_output0,m128i_input0);\
	m128i_output1 = _mm_xor_si128(m128i_output1,m128i_input1);\
	m128i_output2 = _mm_xor_si128(m128i_output2,m128i_input2);\
	m128i_output3 = _mm_xor_si128(m128i_output3,m128i_input3);\
	m128i_output4 = _mm_xor_si128(m128i_output4,m128i_input4);\
	m128i_output5 = _mm_xor_si128(m128i_output5,m128i_input5);\
	m128i_output6 = _mm_xor_si128(m128i_output6,m128i_input6);\
	m128i_output7 = _mm_xor_si128(m128i_output7,m128i_input7);\
} while(0);

//transpose 4x4, 32bit value matrix of (a,b,c,d)
#define transpose_4x4_32bit_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_tmp0,m128i_tmp1,m128i_tmp2,m128i_tmp3) do {\
	m128i_tmp0 = _mm_unpacklo_epi32(m128i_input0, m128i_input1);\
	m128i_tmp1 = _mm_unpacklo_epi32(m128i_input2, m128i_input3);\
	m128i_tmp2 = _mm_unpackhi_epi32(m128i_input0, m128i_input1);\
	m128i_tmp3 = _mm_unpackhi_epi32(m128i_input2, m128i_input3);\
	\
	m128i_input0 = _mm_unpacklo_epi64(m128i_tmp0, m128i_tmp1);\
	m128i_input1 = _mm_unpackhi_epi64(m128i_tmp0, m128i_tmp1);\
	m128i_input2 = _mm_unpacklo_epi64(m128i_tmp2, m128i_tmp3);\
	m128i_input3 = _mm_unpackhi_epi64(m128i_tmp2, m128i_tmp3);\
} while(0);

#define transpose_8x8_16bit_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_tmp0,m128i_tmp1,m128i_tmp2,m128i_tmp3,m128i_tmp4,m128i_tmp5,m128i_tmp6,m128i_tmp7) do {\
	\
	m128i_tmp0 = _mm_unpacklo_epi16(m128i_input0, m128i_input1);\
	m128i_tmp1 = _mm_unpacklo_epi16(m128i_input2, m128i_input3);\
	m128i_tmp2 = _mm_unpacklo_epi16(m128i_input4, m128i_input5);\
	m128i_tmp3 = _mm_unpacklo_epi16(m128i_input6, m128i_input7);\
	m128i_tmp4 = _mm_unpackhi_epi16(m128i_input0, m128i_input1);\
	m128i_tmp5 = _mm_unpackhi_epi16(m128i_input2, m128i_input3);\
	m128i_tmp6 = _mm_unpackhi_epi16(m128i_input4, m128i_input5);\
	m128i_tmp7 = _mm_unpackhi_epi16(m128i_input6, m128i_input7);\
	\
	transpose_4x4_32bit_macro(m128i_tmp0,m128i_tmp1,m128i_tmp2,m128i_tmp3,m128i_input0,m128i_input1,m128i_input2,m128i_input3);\
	transpose_4x4_32bit_macro(m128i_tmp4,m128i_tmp5,m128i_tmp6,m128i_tmp7,m128i_input0,m128i_input1,m128i_input2,m128i_input3);\
	\
	set_8x128_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_tmp0,m128i_tmp1,m128i_tmp2,m128i_tmp3,m128i_tmp4,m128i_tmp5,m128i_tmp6,m128i_tmp7);\
} while(0);


#endif