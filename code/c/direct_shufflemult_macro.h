#ifndef _MULT_DIRECT_SHUFFLEMULT_MACRO
#define _MULT_DIRECT_SHUFFLEMULT_MACRO

#define anubis_mult_4_shuffle_macro(m128i_in,m128i_in_high,m128i_tmp,m128i_ret) do {\
	\
	m128i_ret = m128i_in;\
	\
	m128i_in_high = _mm_srli_epi16(_mm_and_si128(m128i_in,midword_zeroes),4);\
	\
	m128i_in = _mm_and_si128(m128i_in,shuffle_no_zeroing_mask);\
	\
	m128i_tmp = _mm_shuffle_epi8(mul_x2_low,m128i_in);\
	m128i_tmp = _mm_xor_si128(m128i_tmp,_mm_shuffle_epi8(mul_x2_high,m128i_in_high));\
	m128i_ret = _mm_xor_si128(_mm_shuffle_epi8(m128i_tmp,shuffle_x2_gf8),m128i_ret);\
	\
	m128i_tmp = _mm_shuffle_epi8(mul_x4_low,m128i_in);\
	m128i_tmp = _mm_xor_si128(m128i_tmp,_mm_shuffle_epi8(mul_x4_high,m128i_in_high));\
	m128i_ret = _mm_xor_si128(_mm_shuffle_epi8(m128i_tmp,shuffle_x4_gf8),m128i_ret);\
	\
	m128i_tmp = _mm_shuffle_epi8(mul_x6_low,m128i_in);\
	m128i_tmp = _mm_xor_si128(m128i_tmp,_mm_shuffle_epi8(mul_x6_high,m128i_in_high));\
	m128i_ret = _mm_xor_si128(_mm_shuffle_epi8(m128i_tmp,shuffle_x6_gf8),m128i_ret);\
} while(0);

#define shuffle_mult_gf16_macro(m128i_ret,m128i_in_0123,m128i_in_4567,m128i_in_89ab,m128i_in_cdef,m128i_mul_lower_0123,m128i_mul_upper_0123,m128i_mul_upper_4567,m128i_mul_lower_89ab,m128i_mul_upper_89ab,m128i_mul_lower_cdef) do {\
	\
	m128i_ret = _mm_shuffle_epi8(m128i_mul_lower_cdef,m128i_in_cdef);\
	m128i_ret = _mm_xor_si128(m128i_ret, _mm_shuffle_epi8(m128i_mul_lower_89ab,m128i_in_89ab));\
	m128i_ret = _mm_xor_si128(m128i_ret, _mm_shuffle_epi8(m128i_mul_upper_4567,m128i_in_4567));\
	m128i_ret = _mm_xor_si128(m128i_ret, _mm_shuffle_epi8(m128i_mul_upper_0123,m128i_in_0123));\
	\
	m128i_ret = _mm_xor_si128(m128i_ret, _mm_srli_epi16(_mm_shuffle_epi8(m128i_mul_upper_89ab,m128i_in_89ab), 8));\
	m128i_ret = _mm_xor_si128(m128i_ret, _mm_slli_epi16(_mm_shuffle_epi8(m128i_mul_lower_0123,m128i_in_0123), 8));\
} while(0);

#define gf16_anubis_mult_2_shuffle_macro(m128i_in,m128i_in_high,m128i_in_0123,m128i_in_4567,m128i_in_89ab,m128i_in_cdef,m128i_tmp,m128i_ret) do {\
	\
	m128i_ret = m128i_in;\
	\
	m128i_in_high = _mm_srli_epi16(_mm_and_si128(m128i_in,midword_zeroes),4);\
	m128i_in = _mm_and_si128(m128i_in,shuffle_no_zeroing_mask);\
	\
	m128i_in_0123 = _mm_and_si128(m128i_in_high,shuffle_zeroing_odd_mask);\
	m128i_in_4567 = _mm_and_si128(m128i_in,shuffle_zeroing_odd_mask);\
	m128i_in_89ab = _mm_and_si128(m128i_in_high,shuffle_zeroing_even_mask);\
	m128i_in_cdef = _mm_and_si128(m128i_in,shuffle_zeroing_even_mask);\
	\
	shuffle_mult_gf16_macro(m128i_tmp,m128i_in_0123,m128i_in_4567,m128i_in_89ab,m128i_in_cdef,mul_x2_lower_0123,mul_x2_upper_0123,mul_x2_upper_4567,mul_x2_lower_89ab,mul_x2_upper_89ab,mul_x2_lower_cdef);\
	m128i_ret = _mm_xor_si128(_mm_shuffle_epi8(m128i_tmp,shuffle_x2_gf16),m128i_ret);\
	\
	shuffle_mult_gf16_macro(m128i_tmp,m128i_in_0123,m128i_in_4567,m128i_in_89ab,m128i_in_cdef,mul_x4_lower_0123,mul_x4_upper_0123,mul_x4_upper_4567,mul_x4_lower_89ab,mul_x4_upper_89ab,mul_x4_lower_cdef);\
	m128i_ret = _mm_xor_si128(_mm_shuffle_epi8(m128i_tmp,shuffle_x4_gf16),m128i_ret);\
	\
	shuffle_mult_gf16_macro(m128i_tmp,m128i_in_0123,m128i_in_4567,m128i_in_89ab,m128i_in_cdef,mul_x6_lower_0123,mul_x6_upper_0123,mul_x6_upper_4567,mul_x6_lower_89ab,mul_x6_upper_89ab,mul_x6_lower_cdef);\
	m128i_ret = _mm_xor_si128(_mm_shuffle_epi8(m128i_tmp,shuffle_x6_gf16),m128i_ret);\
} while(0);

#endif

	



