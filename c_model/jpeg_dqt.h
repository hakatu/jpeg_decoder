#ifndef JPEG_DQT_H
#define JPEG_DQT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

// Zigzag table
static const int m_zigzag_table[] = {
     0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63,
    0
};

// Winograd-specific quantization scale factors
#define DCTSIZE 8
#define DCT_SCALE_BITS 7
#define DCT_SCALE (1U << DCT_SCALE_BITS)
#define WINOGRAD_QUANT_SCALE_BITS 10

static const uint8_t gWinogradQuant[DCTSIZE * DCTSIZE] = {
    128, 178, 178, 167, 246, 167, 151, 232,
    232, 151, 128, 209, 219, 209, 128, 101,
    178, 197, 197, 178, 101,  69, 139, 167,
    177, 167, 139,  69,  35,  96, 131, 151,
    151, 131,  96,  35,  49,  91, 118, 128,
    118,  91,  49,  46,  81, 101, 101,  81,
     46,  42,  69,  79,  69,  42,  35,  54,
     54,  35,  28,  37,  28,  19,  19,  10
};

#define dprintf

//-----------------------------------------------------------------------------
// jpeg_dqt:
//-----------------------------------------------------------------------------
class jpeg_dqt
{
public:
    jpeg_dqt() { reset(); }

    //-------------------------------------------------------------------------
    // reset: Reset DQT tables and initialize Winograd quantization
    //-------------------------------------------------------------------------
    void reset(void)
    {
        memset(&m_table_dqt[0], 0, 64 * 4);
        createWinogradQuant(); // Precompute adjusted quantization tables
    }

    //-------------------------------------------------------------------------
    // process: Store DQT table from input stream
    //-------------------------------------------------------------------------
    int process(uint8_t *data, int len)
    {
        uint8_t *buf = data;

        // Table number
        uint8_t table_num = (*buf++) & 0x3;
        dprintf(" DQT: Table %d\n", table_num);

        for (int x = 0; x < 64; x++)
        {
            // 8-bit
            uint8_t qv = *buf++;
            dprintf(" %d: %x\n", x, qv);
            m_table_dqt[table_num][x] = qv;
        }

        // Update Winograd-adjusted table after loading new DQT
        createWinogradQuant();

        return buf - data;
    }

    //-------------------------------------------------------------------------
    // lookup: DQT table entry lookup (original table)
    //-------------------------------------------------------------------------
    uint8_t lookup(int table_num, int position)
    {
        return m_table_dqt[table_num][position];
    }

    //-------------------------------------------------------------------------
    // process_samples: Multiply out samples with Winograd-adjusted quantization
    // and de-zigzag ready for Winograd IDCT
    // samples: (idx, value)
    //-------------------------------------------------------------------------
    void process_samples(int quant_table, int *sample_in, int *block_out, int count)
    {
        // Apply quantization and zigzag with Winograd scaling
        memset(block_out, 0, sizeof(block_out[0]) * 64);
        for (int i = 0; i < count; i++)
        {
            int16_t smpl = (int16_t)(sample_in[i] & 0xFFFF);
            int block_idx = (sample_in[i] >> 16);
            int qv = m_table_dqt_winograd[quant_table][block_idx]; // Use Winograd-adjusted table
            int value = smpl * qv;
            dprintf("DEQ: %d: %d * %d -> %d @ %d\n", block_idx, smpl, qv, value, m_zigzag_table[block_idx]);
            block_out[m_zigzag_table[block_idx]] = value;
        }
    }

private:
    uint8_t m_table_dqt[4][64];           // Original JPEG quantization tables
    int16_t m_table_dqt_winograd[4][64];  // Winograd-adjusted quantization tables

    //-------------------------------------------------------------------------
    // createWinogradQuant: Adjust quantization tables for Winograd IDCT
    //-------------------------------------------------------------------------
    void createWinogradQuant(void)
    {
        for (int table = 0; table < 4; table++) {
            for (int i = 0; i < 64; i++) {
                long x = m_table_dqt[table][i];
                x *= gWinogradQuant[i];
                m_table_dqt_winograd[table][i] = (int16_t)((x + (1 << (WINOGRAD_QUANT_SCALE_BITS - DCT_SCALE_BITS - 1))) >> 
                                                           (WINOGRAD_QUANT_SCALE_BITS - DCT_SCALE_BITS));
            }
        }
    }
};

#endif