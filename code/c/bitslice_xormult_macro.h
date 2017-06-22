#ifndef _MULT_BITSLICE_XORMULT_MACRO
#define _MULT_BITSLICE_XORMULT_MACRO

// //calculates multiplication by two in the Anubis Field
// //input/output: 8x128 bit bitsliced vectors
#define x2_8x128_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_tmp0) do {\
	m128i_tmp0 = m128i_input0;\
	m128i_input0 = m128i_input1;\
	m128i_input1 = m128i_input2;\
	m128i_input2 = m128i_input3;\
	m128i_input3 = _mm_xor_si128(m128i_input4,m128i_tmp0);\
	m128i_input4 = _mm_xor_si128(m128i_input5,m128i_tmp0);\
	m128i_input5 = _mm_xor_si128(m128i_input6,m128i_tmp0);\
	m128i_input6 = m128i_input7;\
	m128i_input7 = m128i_tmp0;\
} while(0);

// //calculates multiplication by 6/4 in the Anubis Field
// //input/output: 8x128 bit bitsliced vectors
#define x6div4_8x128_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_tmp7,m128i_tmp4xor7) do {\
	m128i_tmp7 = m128i_input7;\
	m128i_input7 = _mm_xor_si128 (m128i_input6,m128i_input7);\
	m128i_input6 = _mm_xor_si128 (m128i_input5,m128i_input7);\
	m128i_tmp4xor7 = _mm_xor_si128 (m128i_input4,m128i_tmp7);\
	m128i_input5 = _mm_xor_si128 (m128i_input5,m128i_tmp4xor7);\
	m128i_input4 = _mm_xor_si128 (m128i_input3,m128i_tmp4xor7);\
	m128i_input3 = _mm_xor_si128 (m128i_input2,m128i_input3);\
	m128i_input2 = _mm_xor_si128 (m128i_input1,m128i_input2);\
	m128i_input1 = _mm_xor_si128 (m128i_input0,m128i_input1);\
	m128i_input0 = _mm_xor_si128 (m128i_input0,m128i_tmp7);\
} while(0);

// //calculates multiplication by four in the Anubis Field
// //input/output: 32x128 bit bitsliced vectors
#define anubis_mult_128_xor_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31,m128i_tmp0,m128i_tmp1) do {\
	\
	set_32x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31);\
	\
	x2_8x128_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_tmp0);\
	x2_8x128_macro(m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_tmp0);\
	x2_8x128_macro(m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_tmp0);\
	x2_8x128_macro(m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31,m128i_tmp0);\
	\
	xor_8x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15);\
	xor_8x128_macro(m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7);\
	xor_8x128_macro(m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31);\
	xor_8x128_macro(m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23);\
	\
	x2_8x128_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_tmp0);\
	x2_8x128_macro(m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_tmp0);\
	x2_8x128_macro(m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_tmp0);\
	x2_8x128_macro(m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31,m128i_tmp0);\
	\
	xor_8x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23);\
	xor_8x128_macro(m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31);\
	xor_8x128_macro(m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7);\
	xor_8x128_macro(m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15);\
	\
	x6div4_8x128_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_tmp0,m128i_tmp1);\
	x6div4_8x128_macro(m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_tmp0,m128i_tmp1);\
	x6div4_8x128_macro(m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_tmp0,m128i_tmp1);\
	x6div4_8x128_macro(m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31,m128i_tmp0,m128i_tmp1);\
	\
	xor_8x128_macro(m128i_output0,m128i_output1,m128i_output2,m128i_output3,m128i_output4,m128i_output5,m128i_output6,m128i_output7,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31);\
	xor_8x128_macro(m128i_output8,m128i_output9,m128i_output10,m128i_output11,m128i_output12,m128i_output13,m128i_output14,m128i_output15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23);\
	xor_8x128_macro(m128i_output16,m128i_output17,m128i_output18,m128i_output19,m128i_output20,m128i_output21,m128i_output22,m128i_output23,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15);\
	xor_8x128_macro(m128i_output24,m128i_output25,m128i_output26,m128i_output27,m128i_output28,m128i_output29,m128i_output30,m128i_output31,m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7);\
} while(0);


#define sse_trans_slice_x4_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3) do {\
	m128i_input0 = sse_trans_slice(m128i_input0);\
	m128i_input1 = sse_trans_slice(m128i_input1);\
	m128i_input2 = sse_trans_slice(m128i_input2);\
	m128i_input3 = sse_trans_slice(m128i_input3);\
} while(0);

#define sse_trans_slice_4x8_macro(m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31) do {\
	sse_trans_slice_x4_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3);\
	sse_trans_slice_x4_macro(m128i_input4,m128i_input5,m128i_input6,m128i_input7);\
	sse_trans_slice_x4_macro(m128i_input8,m128i_input9,m128i_input10,m128i_input11);\
	sse_trans_slice_x4_macro(m128i_input12,m128i_input13,m128i_input14,m128i_input15);\
	sse_trans_slice_x4_macro(m128i_input16,m128i_input17,m128i_input18,m128i_input19);\
	sse_trans_slice_x4_macro(m128i_input20,m128i_input21,m128i_input22,m128i_input23);\
	sse_trans_slice_x4_macro(m128i_input24,m128i_input25,m128i_input26,m128i_input27);\
	sse_trans_slice_x4_macro(m128i_input28,m128i_input29,m128i_input30,m128i_input31);\
} while(0);

#define sse_trans_slice_x4_rev_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3) do {\
	m128i_input0 = sse_trans_slice_rev(m128i_input0);\
	m128i_input1 = sse_trans_slice_rev(m128i_input1);\
	m128i_input2 = sse_trans_slice_rev(m128i_input2);\
	m128i_input3 = sse_trans_slice_rev(m128i_input3);\
} while(0);

#define sse_trans_slice_4x8_rev_macro(m128i_input1,m128i_input2,m128i_input3,m128i_input4,m128i_input5,m128i_input6,m128i_input7,m128i_input8,m128i_input9,m128i_input10,m128i_input11,m128i_input12,m128i_input13,m128i_input14,m128i_input15,m128i_input16,m128i_input17,m128i_input18,m128i_input19,m128i_input20,m128i_input21,m128i_input22,m128i_input23,m128i_input24,m128i_input25,m128i_input26,m128i_input27,m128i_input28,m128i_input29,m128i_input30,m128i_input31) do {\
	sse_trans_slice_x4_rev_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3);\
	sse_trans_slice_x4_rev_macro(m128i_input4,m128i_input5,m128i_input6,m128i_input7);\
	sse_trans_slice_x4_rev_macro(m128i_input8,m128i_input9,m128i_input10,m128i_input11);\
	sse_trans_slice_x4_rev_macro(m128i_input12,m128i_input13,m128i_input14,m128i_input15);\
	sse_trans_slice_x4_rev_macro(m128i_input16,m128i_input17,m128i_input18,m128i_input19);\
	sse_trans_slice_x4_rev_macro(m128i_input20,m128i_input21,m128i_input22,m128i_input23);\
	sse_trans_slice_x4_rev_macro(m128i_input24,m128i_input25,m128i_input26,m128i_input27);\
	sse_trans_slice_x4_rev_macro(m128i_input28,m128i_input29,m128i_input30,m128i_input31);\
} while(0);


#define bytesliceArrange_128x4_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3) do {\
	\
	m128i_input0 = _mm_shuffle_epi8(m128i_input0,transpose_4x4_8bit_mask);\
	m128i_input1 = _mm_shuffle_epi8(m128i_input1,transpose_4x4_8bit_mask);\
	m128i_input2 = _mm_shuffle_epi8(m128i_input2,transpose_4x4_8bit_mask);\
	m128i_input3 = _mm_shuffle_epi8(m128i_input3,transpose_4x4_8bit_mask);\
	\
	transpose_4x4_32bit_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_tmp0,m128i_tmp1,m128i_tmp2,m128i_tmp3);\
} while(0);

#define bytesliceArrange_128x4_rev_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3) do {\
	\
	transpose_4x4_32bit_macro(m128i_input0,m128i_input1,m128i_input2,m128i_input3,m128i_tmp0,m128i_tmp1,m128i_tmp2,m128i_tmp3);\
	\
	m128i_input0 = _mm_shuffle_epi8(m128i_input0,transpose_4x4_8bit_mask);\
	m128i_input1 = _mm_shuffle_epi8(m128i_input1,transpose_4x4_8bit_mask);\
	m128i_input2 = _mm_shuffle_epi8(m128i_input2,transpose_4x4_8bit_mask);\
	m128i_input3 = _mm_shuffle_epi8(m128i_input3,transpose_4x4_8bit_mask);\
} while(0);

#endif