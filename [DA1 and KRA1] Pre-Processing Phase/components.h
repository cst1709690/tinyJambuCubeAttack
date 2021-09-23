/*
     TinyJAMBU-128: 128-bit key, 96-bit IV
	 Reference implementation for 32-bit CPU
     The state consists of four 32-bit registers
     state[3] || state[2] || state[1] || state[0]

	 Implemented by: Hongjun Wu
*/

#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cstdint>
#include <iostream>

#define FrameBitsIV  0x10
#define FrameBitsAD  0x30
#define FrameBitsPC  0x50  //Framebits for plaintext/ciphertext
#define FrameBitsFinalization 0x70

struct TinyJambu {
    int NROUND1;
    int NROUND2;
    uint8_t mac[8];
	uint32_t state[4];

    TinyJambu(int, int);
    uint32_t state_update(uint32_t *, const uint8_t *, uint32_t);
    uint32_t initialization(const uint8_t *, const uint8_t *, uint32_t *);
    void process_ad(const uint8_t *, const uint8_t *, uint64_t, uint32_t *);
    int crypto_aead_encrypt(uint8_t *, uint64_t *, const uint8_t *, uint64_t, const uint8_t *, uint64_t, const uint8_t *, const uint8_t *, const uint8_t *);
};

TinyJambu::TinyJambu(int NROUND1 = 384, int NROUND2 = 1024)
{
    this -> NROUND1 = NROUND1;
    this -> NROUND2 = NROUND2;
}

/*no-optimized date update function*/
uint32_t TinyJambu::state_update(uint32_t *state, const uint8_t *key, uint32_t number_of_steps)
{
	uint32_t i;
	uint32_t t1, t2, t3, t4, feedback;
	for (i = 0; i < (number_of_steps >> 5); i++)
	{
		t1 = (state[1] >> 15) | (state[2] << 17);  // 47 = 1*32+15
		t2 = (state[2] >> 6)  | (state[3] << 26);  // 47 + 23 = 70 = 2*32 + 6
		t3 = (state[2] >> 21) | (state[3] << 11);  // 47 + 23 + 15 = 85 = 2*32 + 21
		t4 = (state[2] >> 27) | (state[3] << 5);   // 47 + 23 + 15 + 6 = 91 = 2*32 + 27
		feedback = state[0] ^ t1 ^ (~(t2 & t3)) ^ t4 ^ ((uint32_t*)key)[i & 3];
 		// shift 32 bit positions
		state[0] = state[1]; state[1] = state[2]; state[2] = state[3];
		state[3] = feedback ;
	}

	return state[2];
}

// The initialization
/* The input to initialization is the 128-bit key; 96-bit IV;*/
uint32_t TinyJambu::initialization(const uint8_t *key, const uint8_t *iv, uint32_t *state)
{
    int i;
    uint32_t ks;

    //initialize the state as 0
    for (i = 0; i < 4; i++) state[i] = 0;

    //update the state with the key
    ks = state_update(state, key, NROUND2);

    //introduce IV into the state
    for (i = 0;  i < 3; i++)
    {
        state[1] ^= FrameBitsIV;
        ks = state_update(state, key, NROUND1);
        state[3] ^= ((uint32_t*)iv)[i];
    }

    return ks;
}

//process the associated data
void TinyJambu::process_ad(const uint8_t *k, const uint8_t *ad, uint64_t adlen, uint32_t *state)
{
	uint64_t i;
	uint32_t j;

	for (i = 0; i < (adlen >> 2); i++)
	{
		state[1] ^= FrameBitsAD;
		state_update(state, k, NROUND1);
		state[3] ^= ((uint32_t*)ad)[i];
	}

	// if adlen is not a multiple of 4, we process the remaining bytes
	if ((adlen & 3) > 0)
	{
		state[1] ^= FrameBitsAD;
		state_update(state, k, NROUND1);
		for (j = 0; j < (adlen & 3); j++)  ((uint8_t*)state)[12 + j] ^= ad[(i << 2) + j];
		state[1] ^= adlen & 3;
	}
}

//encrypt plaintext
int TinyJambu::crypto_aead_encrypt(
	uint8_t *c,uint64_t *clen,
	const uint8_t *m,uint64_t mlen,
	const uint8_t *ad,uint64_t adlen,
	const uint8_t *nsec,
	const uint8_t *npub,
	const uint8_t *k
	)
{
    uint64_t i;
	uint32_t j;

    //initialization stage
    printf("\nInitalisation:");
    initialization(k, npub, state);

    //process the associated data
    printf("\nAD Processing:");
	process_ad(k, ad, adlen, state);

	//process the plaintext
	printf("\nPlaintext Processing");
	for (i = 0; i < (mlen >> 2); i++)
	{
		state[1] ^= FrameBitsPC;
		state_update(state, k, NROUND2);
		state[3] ^= ((uint32_t*)m)[i];
		((uint32_t*)c)[i] = state[2] ^ ((uint32_t*)m)[i];
	}
	// if mlen is not a multiple of 4, we process the remaining bytes
	if ((mlen & 3) > 0)
	{
		state[1] ^= FrameBitsPC;
		state_update(state, k, NROUND2);
		for (j = 0; j < (mlen & 3); j++)
		{
			((uint8_t*)state)[12 + j] ^= m[(i << 2) + j];
			c[(i << 2) + j] = ((uint8_t*)state)[8 + j] ^ m[(i << 2) + j];
		}
		state[1] ^= mlen & 3;
	}

	//finalization stage, we assume that the tag length is 8 bytes
	printf("\nFinalisation:");
	state[1] ^= FrameBitsFinalization;
	state_update(state, k, NROUND2);
	((uint32_t*)mac)[0] = state[2];

	state[1] ^= FrameBitsFinalization;
	state_update(state, k, NROUND1);
	((uint32_t*)mac)[1] = state[2];

    *clen = mlen + 8;
    memcpy(c + mlen, mac, 8);

    return 0;
}
#endif
