// An implementation of the Discrete Cosine Algorithm for JPEG encoding
#include "dct.h"

//float c_fl[N] = {1,0.981,0.924,0.831,0.707,0.556,0.383,0.195};
//Fixed point: Scale Factor = 2^15

// coefficients
uint16_t c[8] = {32768, 32138, 30274, 27246, 23170, 18205, 12540, 6393};

uint32_t butterfly (uint32_t i , uint32_t c) {
	int16_t i1 = i & 0xffff;
	int16_t i0 = (i >> 16);
	int16_t c1 = c & 0xffff;
	int16_t c0 = (c >> 16);
	int16_t out0, out1;
	out0 = c0 * (i0 + i1);
	out1 = c1 * (i0 - i1);
	return (out0 << 16) | (out1);
}

inline uint32_t squash(int16_t a, int16_t b) {
    return ((uint32_t)a << 16) | ((uint32_t)b & 0xffff);
}

uint32_t quadra(int16_t *in, int16_t *co) {
    // refer to report for definition
    uint32_t out_a = butterfly(squash(in[0],in[1]), squash(co[0],co[1]));
    uint32_t out_b = butterfly(squash(in[2],in[3]), squash(co[2],co[3]));
    uint32_t out_c = butterfly(squash(in[4],in[5]), squash(co[4],co[5]));
    uint32_t out_d = butterfly(squash(in[6],in[7]), squash(co[6],co[7]));

    int16_t out_0 = (int16_t)(out_a>>16) + (int16_t)(out_b>>16) + 
        (int16_t)(out_c>>16) + (int16_t)(out_d>>16);

    int16_t out_1 = (int16_t)(out_a & 0xffff) + (int16_t)(out_b & 0xffff) + 
        (int16_t)(out_c & 0xffff) + (int16_t)(out_d & 0xffff);

    return squash(out_0, out_1);
}

void dct_1d(uint8_t *x, uint8_t *out) {
    uint16_t in[N] = {x[0],x[7],x[1],x[6],x[2],x[5],x[3],x[4]};

    int16_t coeffs_0_1[N] = {c[4],c[1],c[4],c[3],c[4],c[5],c[4],c[7]};
    int16_t coeffs_2_3[N] = {c[2],c[3],c[6],-c[7],-c[6],-c[1],-c[2],-c[5]};
    int16_t coeffs_4_5[N] = {c[4],c[5],-c[4],-c[1],-c[4],c[7],c[4],c[3]};
    int16_t coeffs_6_7[N] = {c[6],c[7],c[2],-c[5],c[2],c[3],-c[3],-c[1]};

    uint32_t out_0_1 = quadra(in, coeffs_0_1);
    uint32_t out_2_3 = quadra(in, coeffs_2_3);
    uint32_t out_4_5 = quadra(in, coeffs_4_5);
    uint32_t out_6_7 = quadra(in, coeffs_6_7);

    out[0] = (int16_t)(out_0_1 >> 16);
    out[1] = (int16_t)(out_0_1 & 0xffff);
    out[2] = (int16_t)(out_2_3 >> 16);
    out[3] = (int16_t)(out_2_3 & 0xffff);
    out[4] = (int16_t)(out_4_5 >> 16);
    out[5] = (int16_t)(out_4_5 & 0xffff);
    out[6] = (int16_t)(out_6_7 >> 16);
    out[7] = (int16_t)(out_6_7 & 0xffff);
}

void dct(uint8_t *in, uint8_t *out) {
    // `in` => 8x8 block of input
    // `out` => 8x8 block of the DCT of `in`
    uint8_t inter[8][8];     // intermediary matrix
    int r,i,j;
    uint8_t tmp;

    // 1D-DCT across rows
    for(r=0; r<8; r++) {
        dct_1d(d[r], inter[r]);
    }

    // transpose h
    for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
            tmp = inter[i][j];
            inter[i][j] = inter[j][i];
            inter[j][i] = tmp;
        }
    }

    // 1D-DCT across columns
    for(r=0; r<N; r++) {
        dct_1d(inter[r], outfix[r]);
    }

    // converting back to floating point
    for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
            out[j][i] = ((float)outfix[i][j] / (2^15));
        }
    }
}
