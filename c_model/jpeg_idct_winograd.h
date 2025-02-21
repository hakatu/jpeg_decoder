/**********************************************************/
/* inverse two dimensional DCT, Winograd algorithm        */
/* 5 multiplies per row/col, up to 80 muls for 2D IDCT   */
/**********************************************************/

#ifndef JPEG_IDCT_WINOGRAD_H
#define JPEG_IDCT_WINOGRAD_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#ifndef DCTSIZE
#define DCTSIZE 8
#endif
#define COLADDR(col) ((col) * DCTSIZE)

// Scaling and helper macros
#define DCT_SCALE_BITS 7
#define DCT_SCALE (1U << DCT_SCALE_BITS)
#define DESCALE(x) (((int32_t)(x) + (1 << (DCT_SCALE_BITS - 1))) >> DCT_SCALE_BITS) // Updated for arithmetic shift
#define WINOGRAD_QUANT_SCALE_BITS 10

//#define FINAL_SCALE_BITS 13
//#define DESCALE(x) (((int32_t)(x) + (1 << (FINAL_SCALE_BITS - 1))) >> FINAL_SCALE_BITS)

class jpeg_idct_winograd
{
public:
    jpeg_idct_winograd() { assert(DCTSIZE == 8); reset(); }
    void reset(void) { memcpy(m_quant, gWinogradQuant, sizeof(gWinogradQuant)); }

private:
    void rowIDCT(int* blk) {
        int src0, src1, src2, src3, src4, src5, src6, src7;
        int x4, x5, x6, x7, x12, x13, x15, x17, x24, x30, x31, x32;
        int x40, x41, x42, x43, x44, tmp1, tmp2, tmp3, stg26;

        // Load inputs
        src0 = blk[0];
        src1 = blk[4];
        src2 = blk[2];
        src3 = blk[6];
        src4 = blk[5];
        src5 = blk[1];
        src6 = blk[3];
        src7 = blk[7];

        // Short-circuit if only DC component is non-zero
        if (!(src1 | src2 | src3 | src4 | src5 | src6 | src7)) {
            int dc = src0 << 3; // Similar scaling to Chen-Wang
            for (int i = 0; i < DCTSIZE; ++i) {
                blk[i] = dc;
            }
            //printf("rowIDCT (DC only): ");
            //for (int i = 0; i < DCTSIZE; i++) {
            //    printf("%d ", blk[i]);
            //}
            //printf("\n");
            return;
        }

        // Stage 1: Compute intermediate values
        x4 = src4 - src7;
        x7 = src4 + src7;
        x5 = src5 + src6;
        x6 = src5 - src6;

        // Stage 2: Multiplications
        tmp1 = imul_b5(x4 - x6);
        stg26 = imul_b4(x6) - tmp1;
        x24 = tmp1 - imul_b2(x4);
        x15 = x5 - x7;
        x17 = x5 + x7;

        tmp2 = stg26 - x17;
        tmp3 = imul_b1_b3(x15) - tmp2;
        x44 = tmp3 + x24;

        // Stage 3: Even coefficients
        x30 = src0 + src1;
        x31 = src0 - src1;
        x12 = src2 - src3;
        x13 = src2 + src3;

        x32 = imul_b1_b3(x12) - x13;

        // Stage 4: Final combinations
        x40 = x30 + x13;
        x43 = x30 - x13;
        x41 = x31 + x32;
        x42 = x31 - x32;

        // Output (no descaling yet, done after columns)
        blk[0] = x40 + x17;
        blk[1] = x41 + tmp2;
        blk[2] = x42 + tmp3;
        blk[3] = x43 - x44;
        blk[4] = x43 + x44;
        blk[5] = x42 - tmp3;
        blk[6] = x41 - tmp2;
        blk[7] = x40 - x17;

        // Debug print
        printf("rowIDCT: ");
        for (int i = 0; i < DCTSIZE; i++) {
            printf("%d ", blk[i]);
        }
        printf("\n");
    }

    void colIDCT(const int* blk, int* out, int stride) {
        int src0, src1, src2, src3, src4, src5, src6, src7;
        int x4, x5, x6, x7, x12, x13, x15, x17, x24, x30, x31, x32;
        int x40, x41, x42, x43, x44, tmp1, tmp2, tmp3, stg26;

        // Load column inputs
        src0 = blk[COLADDR(0)];
        src1 = blk[COLADDR(4)];
        src2 = blk[COLADDR(2)];
        src3 = blk[COLADDR(6)];
        src4 = blk[COLADDR(5)];
        src5 = blk[COLADDR(1)];
        src6 = blk[COLADDR(3)];
        src7 = blk[COLADDR(7)];

        // Short-circuit if only DC component is non-zero
        if (!(src1 | src2 | src3 | src4 | src5 | src6 | src7)) {
            int dc = DESCALE(src0 << 3); // Scale DC appropriately
            for (int i = 0; i < DCTSIZE; ++i) {
                out[i * stride] = dc;
            }
            //printf("colIDCT (DC only): ");
            //for (int i = 0; i < DCTSIZE; i++) {
            //    printf("%d ", out[i * stride]);
            //}
            //printf("\n");
            return;
        }

        // Stage 1: Compute intermediate values
        x4 = src4 - src7;
        x7 = src4 + src7;
        x5 = src5 + src6;
        x6 = src5 - src6;

        // Stage 2: Multiplications
        tmp1 = imul_b5(x4 - x6);
        stg26 = imul_b4(x6) - tmp1;
        x24 = tmp1 - imul_b2(x4);
        x15 = x5 - x7;
        x17 = x5 + x7;

        tmp2 = stg26 - x17;
        tmp3 = imul_b1_b3(x15) - tmp2;
        x44 = tmp3 + x24;

        // Stage 3: Even coefficients
        x30 = src0 + src1;
        x31 = src0 - src1;
        x12 = src2 - src3;
        x13 = src2 + src3;

        x32 = imul_b1_b3(x12) - x13;

        // Stage 4: Final combinations
        x40 = x30 + x13;
        x43 = x30 - x13;
        x41 = x31 + x32;
        x42 = x31 - x32;

        // Output with adjusted scaling to match JPEG range (-128 to 127)
        out[0 * stride] = DESCALE(x40 + x17);  // Additional shift to normalize
        out[1 * stride] = DESCALE(x41 + tmp2);
        out[2 * stride] = DESCALE(x42 + tmp3);
        out[3 * stride] = DESCALE(x43 - x44);
        out[4 * stride] = DESCALE(x43 + x44);
        out[5 * stride] = DESCALE(x42 - tmp3);
        out[6 * stride] = DESCALE(x41 - tmp2);
        out[7 * stride] = DESCALE(x40 - x17);

        // Debug print
        printf("colIDCT: ");
        for (int i = 0; i < DCTSIZE; i++) {
            printf("%d ", out[i * stride]);
        }
        printf("\n");
    }

public:
    void process(int* data_in, int* data_out) {
        int coef;

        // Apply quantization scaling
        for (coef = 0; coef < DCTSIZE * DCTSIZE; ++coef) {
            long x = data_in[coef];
            x *= m_quant[coef];
            data_in[coef] = (int)((x + (1 << (WINOGRAD_QUANT_SCALE_BITS - DCT_SCALE_BITS - 1))) >> 
                                 (WINOGRAD_QUANT_SCALE_BITS - DCT_SCALE_BITS));
        }

        // Row-wise IDCT
        for (coef = 0; coef < (DCTSIZE * DCTSIZE); coef += 8) {
            rowIDCT(data_in + coef);
        }

        // Column-wise IDCT
        for (coef = 0; coef < DCTSIZE; ++coef) {
            colIDCT(data_in + coef, data_out + coef, DCTSIZE);
        }

    }

private:
    // Winograd quantization table
    static const uint8_t gWinogradQuant[DCTSIZE * DCTSIZE];
    uint8_t m_quant[DCTSIZE * DCTSIZE]; // Local copy for potential modification

    // Multiplication coefficients (scaled by 2^8)
    static const int W_B1_B3 = 362; // 1/cos(4*pi/16)
    static const int W_B2 = 669;    // 1/cos(6*pi/16)
    static const int W_B4 = 277;    // 1/cos(2*pi/16)
    static const int W_B5 = 196;    // 1/(cos(2*pi/16) + cos(6*pi/16))

    // Multiplication helper functions
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

    // Clamping to 8-bit range
    inline int clamp(int s) {
        if ((unsigned)s > 255U) {
            if (s < 0) return 0;
            else if (s > 255) return 255;
        }
        return s;
    }
};

// Static initialization of Winograd quantization table
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

#endif // JPEG_IDCT_WINOGRAD_H