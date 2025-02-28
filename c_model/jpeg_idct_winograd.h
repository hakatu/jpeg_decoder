#ifndef JPEG_IDCT_WINOGRAD_H
#define JPEG_IDCT_WINOGRAD_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

//-----------------------------------------------------------------------------
// jpeg_idct: Class to perform IDCT using the Winograd algorithm
//-----------------------------------------------------------------------------
class jpeg_idct_winograd
{
public:
    jpeg_idct_winograd() { reset(); }
    void reset(void) { }

    // Winograd multiplication macros, adapted from picojpeg
    // These perform fixed-point multiplications with a right shift of 8
    #define IMUL_B1_B3(w) (((w) * 362) >> 8)  // 1/cos(4pi/16)
    #define IMUL_B2(w)    (((w) * 669) >> 8)  // 1/cos(6pi/16)
    #define IMUL_B4(w)    (((w) * 277) >> 8)  // 1/cos(2pi/16)
    #define IMUL_B5(w)    (((w) * 196) >> 8)  // 1/(cos(2pi/16) + cos(6pi/16))

    #define SHIFT 7  // Descaling shift, same as PJPG_DCT_SCALE_BITS in picojpeg

    //-----------------------------------------------------------------------------
    // process: Perform 2D Winograd IDCT on an 8x8 block of dequantized data
    // Input:  data_in  - 64-element array of dequantized DCT coefficients
    // Output: data_out - 64-element array of spatial domain values
    //-----------------------------------------------------------------------------
    void process(int *data_in, int *data_out)
    {
        int temp_buf[64];  // Temporary buffer for row-wise IDCT results

        // Step 1: IDCT on rows
        for (int i = 0; i < 8; i++) {
            int *pSrc = data_in + i * 8;
            int *pDst = temp_buf + i * 8;

            // Check if all AC coefficients are zero (DC-only case)
            if (pSrc[1] == 0 && pSrc[2] == 0 && pSrc[3] == 0 && pSrc[4] == 0 &&
                pSrc[5] == 0 && pSrc[6] == 0 && pSrc[7] == 0) {
                int src0 = pSrc[0];
                for (int j = 0; j < 8; j++) {
                    pDst[j] = src0;  // Propagate DC value across the row
                }
            } else {
                // Full Winograd IDCT for the row
                int src4 = pSrc[5];
                int src7 = pSrc[3];
                int x4  = src4 - src7;
                int x7  = src4 + src7;

                int src5 = pSrc[1];
                int src6 = pSrc[7];
                int x5  = src5 + src6;
                int x6  = src5 - src6;

                int tmp1 = IMUL_B5(x4 - x6);
                int stg26 = IMUL_B4(x6) - tmp1;
                int x24 = tmp1 - IMUL_B2(x4);

                int x15 = x5 - x7;
                int x17 = x5 + x7;

                int tmp2 = stg26 - x17;
                int tmp3 = IMUL_B1_B3(x15) - tmp2;
                int x44 = tmp3 + x24;

                int src0 = pSrc[0];
                int src1 = pSrc[4];
                int x30 = src0 + src1;
                int x31 = src0 - src1;

                int src2 = pSrc[2];
                int src3 = pSrc[6];
                int x12 = src2 - src3;
                int x13 = src2 + src3;

                int x32 = IMUL_B1_B3(x12) - x13;

                int x40 = x30 + x13;
                int x43 = x30 - x13;
                int x41 = x31 + x32;
                int x42 = x31 - x32;

                pDst[0] = x40 + x17;
                pDst[1] = x41 + tmp2;
                pDst[2] = x42 + tmp3;
                pDst[3] = x43 - x44;
                pDst[4] = x43 + x44;
                pDst[5] = x42 - tmp3;
                pDst[6] = x41 - tmp2;
                pDst[7] = x40 - x17;
            }
        }

        // Step 2: IDCT on columns
        for (int i = 0; i < 8; i++) {
            int *pSrc = temp_buf + i;
            int *pDst = data_out + i;

            // Check if all AC coefficients are zero (DC-only case)
            if (pSrc[8] == 0 && pSrc[16] == 0 && pSrc[24] == 0 && pSrc[32] == 0 &&
                pSrc[40] == 0 && pSrc[48] == 0 && pSrc[56] == 0) {
                int src0 = pSrc[0];
                int val = (src0 + (1 << (SHIFT - 1))) >> SHIFT;  // Descaling with rounding
                for (int j = 0; j < 8; j++) {
                    pDst[j * 8] = val;  // Propagate DC value across the column
                }
            } else {
                // Full Winograd IDCT for the column
                int src4 = pSrc[40];  // pSrc[5*8]
                int src7 = pSrc[24];  // pSrc[3*8]
                int x4  = src4 - src7;
                int x7  = src4 + src7;

                int src5 = pSrc[8];   // pSrc[1*8]
                int src6 = pSrc[56];  // pSrc[7*8]
                int x5  = src5 + src6;
                int x6  = src5 - src6;

                int tmp1 = IMUL_B5(x4 - x6);
                int stg26 = IMUL_B4(x6) - tmp1;
                int x24 = tmp1 - IMUL_B2(x4);

                int x15 = x5 - x7;
                int x17 = x5 + x7;

                int tmp2 = stg26 - x17;
                int tmp3 = IMUL_B1_B3(x15) - tmp2;
                int x44 = tmp3 + x24;

                int src0 = pSrc[0];
                int src1 = pSrc[32];  // pSrc[4*8]
                int x30 = src0 + src1;
                int x31 = src0 - src1;

                int src2 = pSrc[16];  // pSrc[2*8]
                int src3 = pSrc[48];  // pSrc[6*8]
                int x12 = src2 - src3;
                int x13 = src2 + src3;

                int x32 = IMUL_B1_B3(x12) - x13;

                int x40 = x30 + x13;
                int x43 = x30 - x13;
                int x41 = x31 + x32;
                int x42 = x31 - x32;

                // Descaling with rounding, no level shift or clamping
                pDst[0 * 8] = (x40 + x17 + (1 << (SHIFT - 1))) >> SHIFT;
                pDst[1 * 8] = (x41 + tmp2 + (1 << (SHIFT - 1))) >> SHIFT;
                pDst[2 * 8] = (x42 + tmp3 + (1 << (SHIFT - 1))) >> SHIFT;
                pDst[3 * 8] = (x43 - x44 + (1 << (SHIFT - 1))) >> SHIFT;
                pDst[4 * 8] = (x43 + x44 + (1 << (SHIFT - 1))) >> SHIFT;
                pDst[5 * 8] = (x42 - tmp3 + (1 << (SHIFT - 1))) >> SHIFT;
                pDst[6 * 8] = (x41 - tmp2 + (1 << (SHIFT - 1))) >> SHIFT;
                pDst[7 * 8] = (x40 - x17 + (1 << (SHIFT - 1))) >> SHIFT;
            }
        }
    }

private:
    // No constants retained from the original, as Winograd uses different coefficients
};

#endif // JPEG_IDCT_H