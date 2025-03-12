module jpeg_idct_x
(
    // Inputs
    input           clk_i,          // Clock input
    input           rst_i,          // Reset input (active high)
    input           img_start_i,    // Start of image signal (unused in this implementation)
    input           img_end_i,      // End of image signal (unused in this implementation)
    input           inport_valid_i, // Input data valid signal
    input  [15:0]   inport_data0_i, // First input coefficient
    input  [15:0]   inport_data1_i, // Second input coefficient
    input  [15:0]   inport_data2_i, // Third input coefficient
    input  [15:0]   inport_data3_i, // Fourth input coefficient
    input  [2:0]    inport_idx_i,   // Row index for input data

    // Outputs
    output reg      outport_valid_o, // Output data valid signal
    output reg [31:0] outport_data_o, // Output transformed value
    output reg [5:0]  outport_idx_o   // Output index (row + column)
);

// ### Constants for AAN Algorithm
// These are fixed-point coefficients scaled by 2^13
localparam signed [13:0] FIX_0_707106781 = 5793;  // sqrt(2)/2
localparam signed [13:0] FIX_0_382683433 = 3135;  // cos(6π/16)
localparam signed [13:0] FIX_0_541196100 = 4433;  // cos(4π/16 + π/4)
localparam signed [13:0] FIX_1_306562965 = 10703; // cos(4π/16 - π/4)
localparam signed [13:0] FIX_1_847759065 = 15137; // cos(π/8)

// ### State Machine States
localparam STATE_IDLE     = 3'd0, // Waiting for input
           STATE_COLLECT1 = 3'd1, // Collecting first half of coefficients
           STATE_COLLECT2 = 3'd2, // Collecting second half of coefficients
           STATE_COMPUTE  = 3'd3, // Computing IDCT
           STATE_OUTPUT   = 3'd4; // Outputting results

// ### Registers
reg [2:0]  state;                // Current state
reg [15:0] coef [0:7];           // Array to store 8 input coefficients (16-bit each)
reg [31:0] output_values [0:7];  // Array to store 8 output values (32-bit each)
reg [2:0]  row_idx;              // Current row index
reg [2:0]  output_count;         // Counter for outputting values

// ### Input Scaling
// Scale 16-bit coefficients to 32-bit by sign-extending and shifting left by 11
wire signed [31:0] tmp_in [0:7];
genvar i;
generate
    for (i = 0; i < 8; i = i + 1) begin : gen_tmp_in
        assign tmp_in[i] = {{5{coef[i][15]}}, coef[i], 11'd0}; // Sign-extend and shift
    end
endgenerate

// ### Even Part Computations
wire signed [31:0] z1_e = tmp_in[0] + tmp_in[4];
wire signed [31:0] z2_e = tmp_in[0] - tmp_in[4];
wire signed [31:0] z3_e = tmp_in[2] + tmp_in[6];
wire signed [31:0] z4_e = tmp_in[2] - tmp_in[6];

wire signed [31:0] tmp0_e = z1_e + z3_e;
wire signed [31:0] tmp2_e = z1_e - z3_e;

// Multiplication and shifting for even part
wire signed [45:0] mult0_e = z4_e * FIX_0_707106781;         // 32-bit * 14-bit = 46-bit
wire signed [45:0] shifted_mult0_e = mult0_e >>> 13;         // Shift right by 13
wire signed [31:0] shifted_mult0_e_32 = shifted_mult0_e[31:0]; // Extract lower 32 bits

wire signed [31:0] tmp1_e = z2_e + shifted_mult0_e_32;
wire signed [31:0] tmp3_e = z2_e - shifted_mult0_e_32;

// ### Odd Part Computations
wire signed [31:0] z1_o = tmp_in[1] + tmp_in[7];
wire signed [31:0] z2_o = tmp_in[3] + tmp_in[5];
wire signed [31:0] z3_o = tmp_in[1] + tmp_in[5];
wire signed [31:0] z4_o = tmp_in[3] + tmp_in[7];

// Sign-extend z3_o and z4_o to 46 bits for consistent multiplication
wire signed [45:0] z3_o_ext = {{14{z3_o[31]}}, z3_o};
wire signed [45:0] z4_o_ext = {{14{z4_o[31]}}, z4_o};

// Multiplications for odd part
wire signed [45:0] mult1_o = (z3_o_ext - z4_o_ext) * FIX_1_847759065;
wire signed [45:0] mult2_o = z1_o * FIX_0_541196100;
wire signed [45:0] mult3_o = z2_o * FIX_1_306562965;
wire signed [45:0] mult4_o = (z3_o_ext + z4_o_ext) * FIX_0_382683433;

// Shifted multiplication results
wire signed [45:0] shifted_mult1_o = mult1_o >>> 13;
wire signed [45:0] shifted_mult2_o = mult2_o >>> 13;
wire signed [45:0] shifted_mult3_o = mult3_o >>> 13;
wire signed [45:0] shifted_mult4_o = mult4_o >>> 13;

// Extract lower 32 bits
wire signed [31:0] shifted_mult1_o_32 = shifted_mult1_o[31:0];
wire signed [31:0] shifted_mult2_o_32 = shifted_mult2_o[31:0];
wire signed [31:0] shifted_mult3_o_32 = shifted_mult3_o[31:0];
wire signed [31:0] shifted_mult4_o_32 = shifted_mult4_o[31:0];

// Compute z5_o for odd part
wire signed [31:0] z5_o = shifted_mult1_o_32;

// Compute temporary values for odd part
wire signed [31:0] tmp7_o = z1_o + z2_o;
wire signed [31:0] tmp6_o = shifted_mult2_o_32 + z5_o;
wire signed [31:0] tmp5_o = shifted_mult3_o_32 - z5_o;
wire signed [31:0] tmp4_o = shifted_mult4_o_32;

// ### Final Output Computations
// Combine even and odd parts, shift right by 8 for final scaling
wire signed [31:0] out0 = (tmp0_e + tmp7_o) >>> 8;
wire signed [31:0] out1 = (tmp1_e + tmp6_o) >>> 8;
wire signed [31:0] out2 = (tmp2_e + tmp5_o) >>> 8;
wire signed [31:0] out3 = (tmp3_e + tmp4_o) >>> 8;
wire signed [31:0] out4 = (tmp3_e - tmp4_o) >>> 8;
wire signed [31:0] out5 = (tmp2_e - tmp5_o) >>> 8;
wire signed [31:0] out6 = (tmp1_e - tmp6_o) >>> 8;
wire signed [31:0] out7 = (tmp0_e - tmp7_o) >>> 8;

// ### Sequential Logic (State Machine)
always @(posedge clk_i) begin
    if (rst_i) begin
        // Reset state
        state <= STATE_IDLE;
        outport_valid_o <= 1'b0;
        outport_data_o <= 32'd0;
        outport_idx_o <= 6'd0;
        output_count <= 3'd0;
    end else begin
        case (state)
            STATE_IDLE: begin
                outport_valid_o <= 1'b0;
                if (inport_valid_i) begin
                    // Collect first 4 coefficients
                    coef[0] <= inport_data0_i;
                    coef[1] <= inport_data1_i;
                    coef[2] <= inport_data2_i;
                    coef[3] <= inport_data3_i;
                    row_idx <= inport_idx_i;
                    state <= STATE_COLLECT2;
                end
            end

            STATE_COLLECT2: begin
                if (inport_valid_i && inport_idx_i == row_idx) begin
                    // Collect last 4 coefficients
                    coef[4] <= inport_data0_i;
                    coef[5] <= inport_data1_i;
                    coef[6] <= inport_data2_i;
                    coef[7] <= inport_data3_i;
                    state <= STATE_COMPUTE;
                end
            end

            STATE_COMPUTE: begin
                // Store computed outputs
                output_values[0] <= out0;
                output_values[1] <= out1;
                output_values[2] <= out2;
                output_values[3] <= out3;
                output_values[4] <= out4;
                output_values[5] <= out5;
                output_values[6] <= out6;
                output_values[7] <= out7;
                output_count <= 3'd0;
                state <= STATE_OUTPUT;
            end

            STATE_OUTPUT: begin
                // Output one value per cycle
                outport_valid_o <= 1'b1;
                outport_data_o <= output_values[output_count];
                outport_idx_o <= {row_idx, output_count}; // 6-bit index: [row_idx, column_idx]
                if (output_count < 7) begin
                    output_count <= output_count + 1;
                end else begin
                    outport_valid_o <= 1'b0;
                    state <= STATE_IDLE;
                end
            end

            default: state <= STATE_IDLE;
        endcase
    end
end

//State machine debug

//always @(posedge clk_i) begin
//    $display("Time: %t, State_y: %d, rst_i: %b, inport_valid_i: %b, inport_idx_i: %d, col_idx: %d",
//             $time, state, rst_i, inport_valid_i, inport_idx_i, row_idx);
//end

endmodule