#ifndef JPEG_IDCT_WINOGRAD_H
#define JPEG_IDCT_WINOGRAD_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#define DCTSIZE 8
#define DCT_SCALE_BITS 7
#define DCT_SCALE (1U << DCT_SCALE_BITS)
#define DESCALE(x) (((int32_t)(x) + (1 << (DCT_SCALE_BITS - 1))) >> DCT_SCALE_BITS)
#define WINOGRAD_QUANT_SCALE_BITS 10

class jpeg_idct_winograd
{
public:
    jpeg_idct_winograd() { assert(DCTSIZE == 8); reset(); }
    void reset(void) { memcpy(m_quant, gWinogradQuant, sizeof(gWinogradQuant)); }

    void process(int* data_in, int* data_out) {
        // Remove redundant quantization; assume data_in is already dequantized
        int* pSrc = data_in;

        // Row-wise IDCT
        for (int i = 0; i < DCTSIZE; i++) {
            rowIDCT(pSrc);
            pSrc += DCTSIZE;
        }

        // Column-wise IDCT
        pSrc = data_in;
        for (int i = 0; i < DCTSIZE; i++) {
            colIDCT(pSrc, data_out);
            pSrc++;
            data_out++;
        }
    }

private:
    void rowIDCT(int* blk) {
        int src0 = blk[0], src1 = blk[4], src2 = blk[2], src3 = blk[6];
        int src4 = blk[5], src5 = blk[1], src6 = blk[3], src7 = blk[7];

        if (!(src1 | src2 | src3 | src4 | src5 | src6 | src7)) {
            int dc = src0;
            for (int i = 0; i < DCTSIZE; ++i) blk[i] = dc;
            return;
        }

        int x4 = src4 - src7, x7 = src4 + src7;
        int x5 = src5 + src6, x6 = src5 - src6;
        int tmp1 = imul_b5(x4 - x6);
        int stg26 = imul_b4(x6) - tmp1;
        int x24 = tmp1 - imul_b2(x4);
        int x15 = x5 - x7, x17 = x5 + x7;
        int tmp2 = stg26 - x17;
        int tmp3 = imul_b1_b3(x15) - tmp2;
        int x44 = tmp3 + x24;

        int x30 = src0 + src1, x31 = src0 - src1;
        int x12 = src2 - src3, x13 = src2 + src3;
        int x32 = imul_b1_b3(x12) - x13;

        int x40 = x30 + x13, x43 = x30 - x13;
        int x41 = x31 + x32, x42 = x31 - x32;

        blk[0] = x40 + x17; blk[1] = x41 + tmp2;
        blk[2] = x42 + tmp3; blk[3] = x43 - x44;
        blk[4] = x43 + x44; blk[5] = x42 - tmp3;
        blk[6] = x41 - tmp2; blk[7] = x40 - x17;

        printf("rowIDCT: ");
        for (int i = 0; i < DCTSIZE; i++) printf("%d ", blk[i]);
        printf("\n");
    }

    void colIDCT(const int* blk, int* out) {
        int src0 = blk[0], src1 = blk[4*8], src2 = blk[2*8], src3 = blk[6*8];
        int src4 = blk[5*8], src5 = blk[1*8], src6 = blk[3*8], src7 = blk[7*8];

        if (!(src1 | src2 | src3 | src4 | src5 | src6 | src7)) {
            int dc = DESCALE(src0);
            for (int i = 0; i < DCTSIZE; ++i) out[i * DCTSIZE] = dc;
            return;
        }

        int x4 = src4 - src7, x7 = src4 + src7;
        int x5 = src5 + src6, x6 = src5 - src6;
        int tmp1 = imul_b5(x4 - x6);
        int stg26 = imul_b4(x6) - tmp1;
        int x24 = tmp1 - imul_b2(x4);
        int x15 = x5 - x7, x17 = x5 + x7;
        int tmp2 = stg26 - x17;
        int tmp3 = imul_b1_b3(x15) - tmp2;
        int x44 = tmp3 + x24;

        int x30 = src0 + src1, x31 = src0 - src1;
        int x12 = src2 - src3, x13 = src2 + src3;
        int x32 = imul_b1_b3(x12) - x13;

        int x40 = x30 + x13, x43 = x30 - x13;
        int x41 = x31 + x32, x42 = x31 - x32;

        out[0 * DCTSIZE] = DESCALE(x40 + x17);
        out[1 * DCTSIZE] = DESCALE(x41 + tmp2);
        out[2 * DCTSIZE] = DESCALE(x42 + tmp3);
        out[3 * DCTSIZE] = DESCALE(x43 - x44);
        out[4 * DCTSIZE] = DESCALE(x43 + x44);
        out[5 * DCTSIZE] = DESCALE(x42 - tmp3);
        out[6 * DCTSIZE] = DESCALE(x41 - tmp2);
        out[7 * DCTSIZE] = DESCALE(x40 - x17);

        printf("colIDCT: ");
        for (int i = 0; i < DCTSIZE; i++) printf("%d ", out[i * DCTSIZE]);
        printf("\n");
    }

    static const int W_B1_B3 = 362, W_B2 = 669, W_B4 = 277, W_B5 = 196;
    static const uint8_t gWinogradQuant[DCTSIZE * DCTSIZE];
    uint8_t m_quant[DCTSIZE * DCTSIZE];

    inline int imul_b1_b3(int w) {
        long x = (w * W_B1_B3) + 128;
        return (int)(x >> 8);
    }

    inline int imul_b2(int w) {
        long x = (w * W_B2) + 128;
        return (int)(x >> 8);
    }

    inline int imul_b4(int w) {
        long x = (w * W_B4) + 128;
        return (int)(x >> 8);
    }

    inline int imul_b5(int w) {
        long x = (w * W_B5) + 128;
        return (int)(x >> 8);
    }
};

const uint8_t jpeg_idct_winograd::gWinogradQuant[DCTSIZE * DCTSIZE] = {
    128, 178, 178, 167, 246, 167, 151, 232,
    232, 151, 128, 209, 219, 209, 128, 101,
    178, 197, 197, 178, 101,  69, 139, 167,
    177, 167, 139,  69,  35,  96, 131, 151,
    151, 131,  96,  35,  49,  91, 118, 128,
    118,  91,  49,  46,  81, 101, 101,  81,
     46,  42,  69,  79,  69,  42,  35,  54,
     54,  35,  28,  37,  28,  19,  19,  10
};

#endif