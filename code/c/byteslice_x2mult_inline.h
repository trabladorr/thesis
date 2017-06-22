//Anubis matrix is:
//     Hadamard(0x01,0x02,0x04,0x06), in GF(2^8)/(x^8 + x^4 + x^3 + x^2 + 1)=0x11d
//Other Hadamard matrix is:
//     Hadamard(0x01,0x02,0xb0,0xb2), in GF(2^8)/(x^8 + x^6 + x^3 + x^2 + 1)=0x165

#ifndef _MULT_BYTESLICE_X2MULT_INLINE
#define _MULT_BYTESLICE_X2MULT_INLINE

//calculate times two multiplication in GF2^8, with irreducible polynomial irr_poly, for every byte in a
static inline __m128i x2(__m128i a, __m128i irr_poly){
	const __m128i midword_zeroes = _mm_setr_epi8(0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f,0x7f);
	const __m128i zeroes = _mm_setzero_si128();

	__m128i b, c;

	//store a copy, to act as mask for blendv
	b = a;
	// sll shifts 16 bits instead of 8, insert zeroes into the high bit of every second byte
	b = _mm_and_si128 (b, midword_zeroes);

	//shift all bytes by one left
	b = _mm_slli_epi16 (b, 1);

	//set b to the irr poly for every byte that overflowed in the shift
	c = _mm_blendv_epi8 (zeroes, irr_poly, a);

	//return the xor of the two
	return _mm_xor_si128 (b, c);
}


//calculate times 4 (100) multiplication in GF2^8, with irreducible polynomial irr_poly, for every byte in a
static inline __m128i x4(__m128i a, __m128i irr_poly){
	return x2(x2(a, irr_poly), irr_poly);
}

//calculate times 6 (110) multiplication in GF2^8, with irreducible polynomial irr_poly, for every byte in a
static inline __m128i x6(__m128i a, __m128i irr_poly){
	__m128i b = x2(a, irr_poly);
	a = x2(b, irr_poly);
	return _mm_xor_si128 (a, b);
}

//calculate times 176 (10110000) multiplication in GF2^8, with irreducible polynomial irr_poly, for every byte in a
static inline __m128i xb0(__m128i a, __m128i irr_poly){
	//x16 multiplication
	a = x4(x4(a, irr_poly), irr_poly);
	//x8 multiplication of a, equivalent to x128 of input
	__m128i b = x2(x4(a, irr_poly), irr_poly);
	//return  of x16 XOR x32 (x2(x16)) XOR x128
	return _mm_xor_si128 (a, _mm_xor_si128 (x2(a, irr_poly), b));
}

//calculate times 178 (10110010) multiplication in GF2^8, with irreducible polynomial irr_poly, for every byte in a
static inline __m128i xb2(__m128i a, __m128i irr_poly){
	return _mm_xor_si128 (x2(a, irr_poly), xb0(a, irr_poly));
}


//input v contains 16 vectors to mutliply, bytesliced, output ret contains 16 vectors mutliplied, bytesliced
//multiply using times two multiplication and XORs, anubis matrix/field
static inline void anubis_mult_16_x2(__m128i v[],__m128i ret[]){
	const __m128i irr_poly = _mm_setr_epi8(0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d); //ANUBIS irr poly

	ret[0] = v[0];
	ret[0] = _mm_xor_si128 (ret[0], x2(v[1], irr_poly));
	ret[0] = _mm_xor_si128 (ret[0], x4(v[2], irr_poly));
	ret[0] = _mm_xor_si128 (ret[0], x6(v[3], irr_poly));

	ret[1] = x2(v[0], irr_poly);
	ret[1] = _mm_xor_si128 (ret[1], v[1]);
	ret[1] = _mm_xor_si128 (ret[1], x6(v[2], irr_poly));
	ret[1] = _mm_xor_si128 (ret[1], x4(v[3], irr_poly));

	ret[2] = x4(v[0], irr_poly);
	ret[2] = _mm_xor_si128 (ret[2], x6(v[1], irr_poly));
	ret[2] = _mm_xor_si128 (ret[2], v[2]);
	ret[2] = _mm_xor_si128 (ret[2], x2(v[3], irr_poly));

	ret[3] = x6(v[0], irr_poly);
	ret[3] = _mm_xor_si128 (ret[3], x4(v[1], irr_poly));
	ret[3] = _mm_xor_si128 (ret[3], x2(v[2], irr_poly));
	ret[3] = _mm_xor_si128 (ret[3], v[3]);
}

//input v contains 16 vectors to mutliply, bytesliced, output ret contains 16 vectors mutliplied, bytesliced
//multiply using times two multiplication and XORs, other Hadamard matrix/field
static inline void h2_mult_16_x2(__m128i v[],__m128i ret[]){
	const __m128i irr_poly = _mm_setr_epi8(0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65,0x65); //other Hadamard matrix irr poly

	ret[0] = v[0], irr_poly;
	ret[0] = _mm_xor_si128 (ret[0], x2(v[1], irr_poly));
	ret[0] = _mm_xor_si128 (ret[0], xb0(v[2], irr_poly));
	ret[0] = _mm_xor_si128 (ret[0], xb2(v[3], irr_poly));

	ret[1] = x2(v[0], irr_poly);
	ret[1] = _mm_xor_si128 (ret[1], v[1]);
	ret[1] = _mm_xor_si128 (ret[1], xb2(v[2], irr_poly));
	ret[1] = _mm_xor_si128 (ret[1], xb0(v[3], irr_poly));

	ret[2] = xb0(v[0], irr_poly);
	ret[2] = _mm_xor_si128 (ret[2], xb2(v[1], irr_poly));
	ret[2] = _mm_xor_si128 (ret[2], v[2]);
	ret[2] = _mm_xor_si128 (ret[2], x2(v[3], irr_poly));

	ret[3] = xb2(v[0], irr_poly);
	ret[3] = _mm_xor_si128 (ret[3], xb0(v[1], irr_poly));
	ret[3] = _mm_xor_si128 (ret[3], x2(v[2], irr_poly));
	ret[3] = _mm_xor_si128 (ret[3], v[3]);
}

//input/output pointer contains 16 vectors to mutliply
//multiply using times two multiplication and XORs, anubis matrix/field
static inline void anubis_mult_16_x2_mem(__m128i input[],__m128i output[]){

	into_byteslice_4x128(input);
	
	anubis_mult_16_x2(input,output);

	from_byteslice_4x128(output);
}

//input/output pointer contains 16 vectors to multiply
//multiply using times two multiplication and XORs, other Hadamard matrix/field
static inline void h2_mult_16_x2_mem(__m128i input[],__m128i output[]){

	into_byteslice_4x128(input);
	
	h2_mult_16_x2(input,output);

	from_byteslice_4x128(output);
}

#endif