/**
 * @file tinymt64.h
 *
 * @brief Tiny Mersenne Twister only 127 bit internal state
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2011 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdint.h>
#include <inttypes.h>

#define TINYMT64_MEXP 127
#define TINYMT64_SH0 12
#define TINYMT64_SH1 11
#define TINYMT64_SH8 8
#define TINYMT64_MASK UINT64_C(0x7fffffffffffffff)
#define TINYMT64_MUL (1.0 / 18446744073709551616.0)

/*
 * tinymt64 default parameters
 */
#ifndef TINYMT64_MAT1
#  define TINYMT64_MAT1 0x7a840f50
#endif
#ifndef TINYMT64_MAT2
#  define TINYMT64_MAT2 0xf3d8fcf6
#endif
#ifndef TINYMT64_TMAT
#  define TINYMT64_TMAT 0x9746beffffbffffe
#endif

static const uint32_t mat1 = TINYMT64_MAT1;
static const uint32_t mat2 = TINYMT64_MAT2;
static const uint64_t tmat = TINYMT64_TMAT;

void tinymt64_init(uint64_t * random, uint64_t seed);
#if 0
void tinymt64_init_by_array(uint64_t * random, const uint64_t init_key[],
			    int key_length);
#endif
uint64_t tinymt64_generate_uint64(uint64_t * random);
double tinymt64_generate_double(uint64_t * random);

#if 0
#if defined(__GNUC__)
/**
 * This function always returns 127
 * @param random not used
 * @return always 127
 */
static int tinymt64_get_mexp(
    uint64_t * random  __attribute__((unused))) {
    return TINYMT64_MEXP;
}
#else
static int tinymt64_get_mexp(uint64_t * random) {
    return TINYMT64_MEXP;
}
#endif
#endif

/**
 * This function changes internal state of tinymt64.
 * Users should not call this function directly.
 * @param random tinymt internal status
 */
static void tinymt64_next_state(uint64_t * random) {
    uint64_t x;

    random[0] &= TINYMT64_MASK;
    x = random[0] ^ random[1];
    x ^= x << TINYMT64_SH0;
    x ^= x >> 32;
    x ^= x << 32;
    x ^= x << TINYMT64_SH1;
    random[0] = random[1];
    random[1] = x;
    random[0] ^= -((int64_t)(x & 1)) & mat1;
    random[1] ^= -((int64_t)(x & 1)) & (((uint64_t)mat2) << 32);
}

/**
 * This function outputs 64-bit unsigned integer from internal state.
 * Users should not call this function directly.
 * @param random tinymt internal status
 * @return 64-bit unsigned pseudorandom number
 */
static uint64_t uint64_temper(uint64_t * random) {
    uint64_t x;
#if defined(LINEARITY_CHECK)
    x = random[0] ^ random[1];
#else
    x = random[0] + random[1];
#endif
    x ^= random[0] >> TINYMT64_SH8;
    x ^= -((int64_t)(x & 1)) & tmat;
    return x;
}

#if 0
/**
 * This function outputs floating point number from internal state.
 * Users should not call this function directly.
 * @param random tinymt internal status
 * @return floating point number r (1.0 <= r < 2.0)
 */
static double uint64_temper_conv(uint64_t * random) {
    uint64_t x;
    union {
	uint64_t u;
	double d;
    } conv;
#if defined(LINEARITY_CHECK)
    x = random[0] ^ random[1];
#else
    x = random[0] + random[1];
#endif
    x ^= random[0] >> TINYMT64_SH8;
    conv.u = ((x ^ (-((int64_t)(x & 1)) & tmat)) >> 12)
	| UINT64_C(0x3ff0000000000000);
    return conv.d;
}

/**
 * This function outputs floating point number from internal state.
 * Users should not call this function directly.
 * @param random tinymt internal status
 * @return floating point number r (1.0 < r < 2.0)
 */
static double uint64_temper_conv_open(uint64_t * random) {
    uint64_t x;
    union {
	uint64_t u;
	double d;
    } conv;
#if defined(LINEARITY_CHECK)
    x = random[0] ^ random[1];
#else
    x = random[0] + random[1];
#endif
    x ^= random[0] >> TINYMT64_SH8;
    conv.u = ((x ^ (-((int64_t)(x & 1)) & tmat)) >> 12)
	| UINT64_C(0x3ff0000000000001);
    return conv.d;
}
#endif

/**
 * This function outputs 64-bit unsigned integer from internal state.
 * @param random tinymt internal status
 * @return 64-bit unsigned integer r (0 <= r < 2^64)
 */
uint64_t tinymt64_generate_uint64(uint64_t * random) {
    tinymt64_next_state(random);
    return uint64_temper(random);
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using multiplying by 1 / 2^64.
 * @param random tinymt internal status
 * @return floating point number r (0.0 <= r < 1.0)
 */
double tinymt64_generate_double(uint64_t * random) {
    tinymt64_next_state(random);
    return uint64_temper(random) * TINYMT64_MUL;
}

#if 0
/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param random tinymt internal status
 * @return floating point number r (0.0 <= r < 1.0)
 */
static double tinymt64_generate_double01(uint64_t * random) {
    tinymt64_next_state(random);
    return uint64_temper_conv(random) - 1.0;
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param random tinymt internal status
 * @return floating point number r (1.0 <= r < 2.0)
 */
static double tinymt64_generate_double12(uint64_t * random) {
    tinymt64_next_state(random);
    return uint64_temper_conv(random);
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param random tinymt internal status
 * @return floating point number r (0.0 < r <= 1.0)
 */
static double tinymt64_generate_doubleOC(uint64_t * random) {
    tinymt64_next_state(random);
    return 2.0 - uint64_temper_conv(random);
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param random tinymt internal status
 * @return floating point number r (0.0 < r < 1.0)
 */
static double tinymt64_generate_doubleOO(uint64_t * random) {
    tinymt64_next_state(random);
    return uint64_temper_conv_open(random) - 1.0;
}
#endif

/**
 * @file tinymt64.c
 *
 * @brief 64-bit Tiny Mersenne Twister only 127 bit internal state
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2011 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#define MIN_LOOP 8

#if 0
/**
 * This function represents a function used in the initialization
 * by init_by_array
 * @param[in] x 64-bit integer
 * @return 64-bit integer
 */
static uint64_t ini_func1(uint64_t x) {
    return (x ^ (x >> 59)) * UINT64_C(2173292883993);
}

/**
 * This function represents a function used in the initialization
 * by init_by_array
 * @param[in] x 64-bit integer
 * @return 64-bit integer
 */
static uint64_t ini_func2(uint64_t x) {
    return (x ^ (x >> 59)) * UINT64_C(58885565329898161);
}

/**
 * This function certificate the period of 2^127-1.
 * @param random tinymt state vector.
 */
static void period_certification(uint64_t * random) {
    if ((random[0] & TINYMT64_MASK) == 0 &&
	random[1] == 0) {
	random[0] = 'T';
	random[1] = 'M';
    }
}
#endif

/**
 * This function initializes the internal state array with a 64-bit
 * unsigned integer seed.
 * @param random tinymt state vector.
 * @param seed a 64-bit unsigned integer used as a seed.
 */
void tinymt64_init(uint64_t * random, uint64_t seed) {
    random[0] = seed ^ ((uint64_t)mat1 << 32);
    random[1] = mat2 ^ tmat;
    for (int i = 1; i < MIN_LOOP; i++) {
	random[i & 1] ^= i + UINT64_C(6364136223846793005)
	    * (random[(i - 1) & 1]
	       ^ (random[(i - 1) & 1] >> 62));
    }
#if 0
    period_certification(random);
#endif
}

#if 0
/**
 * This function initializes the internal state array,
 * with an array of 64-bit unsigned integers used as seeds
 * @param random tinymt state vector.
 * @param init_key the array of 64-bit integers, used as a seed.
 * @param key_length the length of init_key.
 */
void tinymt64_init_by_array(uint64_t * random, const uint64_t init_key[],
			    int key_length) {
    const int lag = 1;
    const int mid = 1;
    const int size = 4;
    int i, j;
    int count;
    uint64_t r;
    uint64_t st[4];

    st[0] = 0;
    st[1] = mat1;
    st[2] = mat2;
    st[3] = tmat;
    if (key_length + 1 > MIN_LOOP) {
	count = key_length + 1;
    } else {
	count = MIN_LOOP;
    }
    r = ini_func1(st[0] ^ st[mid % size]
		  ^ st[(size - 1) % size]);
    st[mid % size] += r;
    r += key_length;
    st[(mid + lag) % size] += r;
    st[0] = r;
    count--;
    for (i = 1, j = 0; (j < count) && (j < key_length); j++) {
	r = ini_func1(st[i] ^ st[(i + mid) % size] ^ st[(i + size - 1) % size]);
	st[(i + mid) % size] += r;
	r += init_key[j] + i;
	st[(i + mid + lag) % size] += r;
	st[i] = r;
	i = (i + 1) % size;
    }
    for (; j < count; j++) {
	r = ini_func1(st[i] ^ st[(i + mid) % size] ^ st[(i + size - 1) % size]);
	st[(i + mid) % size] += r;
	r += i;
	st[(i + mid + lag) % size] += r;
	st[i] = r;
	i = (i + 1) % size;
    }
    for (j = 0; j < size; j++) {
	r = ini_func2(st[i] + st[(i + mid) % size] + st[(i + size - 1) % size]);
	st[(i + mid) % size] ^= r;
	r -= i;
	st[(i + mid + lag) % size] ^= r;
	st[i] = r;
	i = (i + 1) % size;
    }
    random[0] = st[0] ^ st[1];
    random[1] = st[2] ^ st[3];
    period_certification(random);
}
#endif
