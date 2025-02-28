/**********************************************************/
/* Inverse two dimensional DCT, AAN algorithm             */
/* 5 multiplications, 29 additions per 8-point DCT        */
/* Adapted for JPEG decoder at github.com/yakeshi/jpeg_decoder */
/**********************************************************/

#ifndef JPEG_idct_aan_H
#define JPEG_idct_aan_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#ifndef DCTSIZE
#define DCTSIZE 8
#endif
#define COLADDR(col)    ((col) * DCTSIZE)

class idct_aan
{
public:
    idct_aan() { assert(DCTSIZE == 8); reset(); }
    void reset(void) { }

private:
void rowIDCT_aan(int* blk) {
    printf("rowIDCT_aan input: ");
    for (int i = 0; i < 8; i++) printf("%d ", blk[i]);
    printf("\n");

    const int FIX_0_707106781 = 5793;  // sqrt(2)/2 * 2^13
    const int FIX_0_382683433 = 3135;  // cos(6π/16) * 2^13
    const int FIX_0_541196100 = 4433;  // cos(4π/16 + π/4) * 2^13
    const int FIX_1_306562965 = 10703; // cos(4π/16 - π/4) * 2^13
    const int FIX_1_847759065 = 15137; // cos(π/8) * 2^13

    int tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    int z1, z2, z3, z4, z5;

    // Stricter DC-only check (tolerate small noise if needed)
    tmp1 = blk[4]; tmp2 = blk[2]; tmp3 = blk[6];
    tmp4 = blk[1]; tmp5 = blk[5]; tmp6 = blk[3]; tmp7 = blk[7];
    if (!(tmp1 | tmp2 | tmp3 | tmp4 | tmp5 | tmp6 | tmp7)) {
        tmp0 = blk[0] << 3;  // Match idct_fast
        for (int i = 0; i < DCTSIZE; ++i) {
            blk[i] = tmp0;
        }
        return;
    }

    // Full AAN IDCT computation (unchanged for brevity)
    tmp0 = blk[0] << 11;
    tmp1 = blk[4] << 11;
    tmp2 = blk[2] << 11;
    tmp3 = blk[6] << 11;
    tmp4 = blk[1] << 11;
    tmp5 = blk[5] << 11;
    tmp6 = blk[3] << 11;
    tmp7 = blk[7] << 11;

    // Even part
    z1 = tmp0 + tmp1; z2 = tmp0 - tmp1;
    z3 = tmp2 + tmp3; z4 = tmp2 - tmp3;
    tmp0 = z1 + z3; tmp2 = z1 - z3;
    tmp1 = z2 + (z4 * FIX_0_707106781 >> 13);
    tmp3 = z2 - (z4 * FIX_0_707106781 >> 13);

    // Odd part
    z1 = tmp4 + tmp7; z2 = tmp5 + tmp6;
    z3 = tmp4 + tmp6; z4 = tmp5 + tmp7;
    z5 = (z3 - z4) * FIX_1_847759065 >> 13;
    tmp7 = z1 + z2;
    tmp6 = (z1 * FIX_0_541196100 >> 13) + z5;
    tmp5 = (z2 * FIX_1_306562965 >> 13) - z5;
    tmp4 = (z3 + z4) * FIX_0_382683433 >> 13;

    blk[0] = (tmp0 + tmp7) >> 8;
    blk[7] = (tmp0 - tmp7) >> 8;
    blk[1] = (tmp1 + tmp6) >> 8;
    blk[6] = (tmp1 - tmp6) >> 8;
    blk[2] = (tmp2 + tmp5) >> 8;
    blk[5] = (tmp2 - tmp5) >> 8;
    blk[3] = (tmp3 + tmp4) >> 8;
    blk[4] = (tmp3 - tmp4) >> 8;
}

void colIDCT_aan(const int* blk, int* out, int stride) {
    const int FIX_0_707106781 = 5793;  // sqrt(2)/2 * 2^13
    const int FIX_0_382683433 = 3135;  // cos(6π/16) * 2^13
    const int FIX_0_541196100 = 4433;  // cos(4π/16 + π/4) * 2^13
    const int FIX_1_306562965 = 10703; // cos(4π/16 - π/4) * 2^13
    const int FIX_1_847759065 = 15137; // cos(π/8) * 2^13

    int tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    int z1, z2, z3, z4, z5;

    // DC-only check
    if (!(blk[COLADDR(4)] | blk[COLADDR(2)] | blk[COLADDR(6)] |
          blk[COLADDR(1)] | blk[COLADDR(5)] | blk[COLADDR(3)] | blk[COLADDR(7)])) {
        tmp0 = (blk[0] + 32) >> 6;  // Match idct_fast scaling
        for (int i = 0; i < DCTSIZE; ++i) {
            *out = tmp0;
            out += stride;
        }
        return;
    }

    // Full AAN IDCT computation
    tmp0 = blk[COLADDR(0)] << 8;
    tmp1 = blk[COLADDR(4)] << 8;
    tmp2 = blk[COLADDR(2)] << 8;
    tmp3 = blk[COLADDR(6)] << 8;
    tmp4 = blk[COLADDR(1)] << 8;
    tmp5 = blk[COLADDR(5)] << 8;
    tmp6 = blk[COLADDR(3)] << 8;
    tmp7 = blk[COLADDR(7)] << 8;

    // Even part
    z1 = tmp0 + tmp1; z2 = tmp0 - tmp1;
    z3 = tmp2 + tmp3; z4 = tmp2 - tmp3;
    tmp0 = z1 + z3; tmp2 = z1 - z3;
    tmp1 = z2 + (z4 * FIX_0_707106781 >> 13);
    tmp3 = z2 - (z4 * FIX_0_707106781 >> 13);

    // Odd part
    z1 = tmp4 + tmp7; z2 = tmp5 + tmp6;
    z3 = tmp4 + tmp6; z4 = tmp5 + tmp7;
    z5 = (z3 - z4) * FIX_1_847759065 >> 13;
    tmp7 = z1 + z2;
    tmp6 = (z1 * FIX_0_541196100 >> 13) + z5;
    tmp5 = (z2 * FIX_1_306562965 >> 13) - z5;
    tmp4 = (z3 + z4) * FIX_0_382683433 >> 13;

    *out = (tmp0 + tmp7 + 32) >> 6; out += stride;
    *out = (tmp1 + tmp6 + 32) >> 6; out += stride;
    *out = (tmp2 + tmp5 + 32) >> 6; out += stride;
    *out = (tmp3 + tmp4 + 32) >> 6; out += stride;
    *out = (tmp3 - tmp4 + 32) >> 6; out += stride;
    *out = (tmp2 - tmp5 + 32) >> 6; out += stride;
    *out = (tmp1 - tmp6 + 32) >> 6; out += stride;
    *out = (tmp0 - tmp7 + 32) >> 6;
}

public:
    void process(int* data_in, int* data_out) {
        int coef;

        // Apply 1D IDCT to rows
        for (coef = 0; coef < (DCTSIZE * DCTSIZE); coef += 8) {
            rowIDCT_aan(data_in + coef);
        }
        // Apply 1D IDCT to columns
        for (coef = 0; coef < DCTSIZE; ++coef) {
            colIDCT_aan(data_in + coef, data_out + coef, DCTSIZE);
        }
    }
};

#endif // JPEG_idct_aan_H