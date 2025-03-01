// DESCRIPTION: col inverse DCT, ANN algorithm
// LICENSE: GNU Lesser General Public License Version 3 (LGPLv3)
// AUTHOR: Haka & grok 3
// YEAR: 2025

module jpeg_idct_aan_y
//-----------------------------------------------------------------
// Params
//-----------------------------------------------------------------
#(
     parameter OUT_SHIFT        = 14
)
//-----------------------------------------------------------------
// Ports
//-----------------------------------------------------------------
(
    // Inputs
     input           clk_i
    ,input           rst_i
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 31:0]  inport_data0_i
    ,input  [ 31:0]  inport_data1_i
    ,input  [ 31:0]  inport_data2_i
    ,input  [ 31:0]  inport_data3_i
    ,input  [  2:0]  inport_idx_i

    // Outputs
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

// Fixed-point constants scaled by 2^14 (16384)
localparam FIX_0_707106781 = 16'd11585;  // 0.707106781 * 16384
localparam FIX_0_382683433 = 16'd6269;   // 0.382683433 * 16384
localparam FIX_0_541196100 = 16'd8867;   // 0.541196100 * 16384
localparam FIX_1_306562965 = 16'd21407;  // 1.306562965 * 16384
localparam FIX_1_847759065 = 16'd30274;  // 1.847759065 * 16384

// Internal registers
reg signed [31:0] coef_buffer [0:7][0:7];  // 8x8 coefficient buffer
reg [2:0] col_idx;
reg computing;
reg outport_valid_o;
reg [31:0] outport_data_o;
reg [5:0] outport_idx_o;

// AAN IDCT computation task for a single column
task compute_idct;
    input [2:0] col;
    output [31:0] result;
    reg signed [39:0] tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    reg signed [39:0] z1, z2, z3, z4, z5;
    reg signed [39:0] coef [0:7];
    integer i;
begin
    // Load coefficients for the column
    for (i = 0; i < 8; i = i + 1) begin
        coef[i] = coef_buffer[i][col];
    end

    // Even part
    z1 = coef[0] + coef[4];
    z2 = coef[0] - coef[4];
    z3 = coef[2] + coef[6];
    z4 = coef[2] - coef[6];

    tmp0 = z1 + z3;
    tmp2 = z1 - z3;
    tmp1 = z2 + ((z4 * FIX_0_707106781) >>> 14);
    tmp3 = z2 - ((z4 * FIX_0_707106781) >>> 14);

    // Odd part
    z1 = coef[1] + coef[7];
    z2 = coef[5] + coef[3];
    z3 = coef[1] + coef[3];
    z4 = coef[5] + coef[7];

    z5 = ((z3 - z4) * FIX_1_847759065) >>> 14;
    tmp7 = z1 + z2;
    tmp6 = ((z1 * FIX_0_541196100) >>> 14) + z5;
    tmp5 = ((z2 * FIX_1_306562965) >>> 14) - z5;
    tmp4 = ((z3 + z4) * FIX_0_382683433) >>> 14;

    // Combine and scale output (only two outputs per cycle here)
    result[15:0]  = (tmp0 + tmp7) >>> OUT_SHIFT;
    result[31:16] = (tmp0 - tmp7) >>> OUT_SHIFT;
end
endtask

// Main sequential logic
always @(posedge clk_i or posedge rst_i) begin
    if (rst_i) begin
        outport_valid_o <= 1'b0;
        outport_data_o  <= 32'd0;
        outport_idx_o   <= 6'd0;
        col_idx         <= 3'd0;
        computing       <= 1'b0;
        // Clear coefficient buffer
        for (integer i = 0; i < 8; i = i + 1) begin
            for (integer j = 0; j < 8; j = j + 1) begin
                coef_buffer[i][j] <= 32'd0;
            end
        end
    end else begin
        if (inport_valid_i) begin
            // Load four coefficients per cycle into the buffer
            coef_buffer[inport_idx_i][0] <= inport_data0_i;
            coef_buffer[inport_idx_i][1] <= inport_data1_i;
            coef_buffer[inport_idx_i][2] <= inport_data2_i;
            coef_buffer[inport_idx_i][3] <= inport_data3_i;
        end

        if (computing) begin
            compute_idct(col_idx, outport_data_o);
            outport_valid_o <= 1'b1;
            outport_idx_o   <= {col_idx, 3'd0};  // Simplified indexing
            col_idx         <= col_idx + 1;
            if (col_idx == 7) begin
                computing <= 1'b0;
            end
        end else if (img_start_i) begin
            computing <= 1'b1;
            col_idx   <= 3'd0;
            outport_valid_o <= 1'b0;
        end else begin
            outport_valid_o <= 1'b0;
        end
    end
end

endmodule