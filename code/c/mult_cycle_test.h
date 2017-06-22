#ifndef _MULT_CYCLE_TEST_H
#define _MULT_CYCLE_TEST_H

uint64_t test_mult_4_shuffle(int zero);
uint64_t test_mult_4_shuffle_macro(int zero);

uint64_t test_mult_2_shuffle_gf16(int zero);
uint64_t test_mult_2_shuffle_gf16_macro(int zero);

uint64_t test_mult_16_shuffle_anubis(int zero);
uint64_t test_mult_16_shuffle_other(int zero);

uint64_t test_mult_16_shuffle_load_anubis(int zero);
uint64_t test_mult_16_shuffle_load_other(int zero);

uint64_t test_mult_16_x2_anubis(int zero);
uint64_t test_mult_16_x2_other(int zero);

uint64_t test_mult_16_x2_load_anubis(int zero);
uint64_t test_mult_16_x2_load_other(int zero);

uint64_t test_mult_128_xor_anubis(int zero);
uint64_t test_mult_128_xor_anubis_opt_compiler(int zero);
uint64_t test_mult_128_xor_anubis_opt_tmp(int zero);
uint64_t test_mult_128_xor_anubis_macro(int zero);

uint64_t test_mult_128_xor_load_anubis(int zero);

#define DIFFTESTS 17

#endif