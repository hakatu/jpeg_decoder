module jpeg_idct_y
//-----------------------------------------------------------------
// Parameters
//-----------------------------------------------------------------
#(
     parameter OUT_SHIFT        = 15
    ,parameter INPUT_WIDTH      = 32
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
    ,input  [31:0]   inport_data0_i
    ,input  [31:0]   inport_data1_i
    ,input  [31:0]   inport_data2_i
    ,input  [31:0]   inport_data3_i
    ,input  [2:0]    inport_idx_i
    // Outputs
    ,output          outport_valid_o
    ,output [31:0]   outport_data_o
    ,output [5:0]    outport_idx_o
);

//-----------------------------------------------------------------
// Constants for AAN Algorithm (Scaled by 2^13)
//-----------------------------------------------------------------
localparam [13:0] FIX_0_707106781 = 5793;  // sqrt(2)/2
localparam [13:0] FIX_0_382683433 = 3135;  // cos(6π/16)
localparam [13:0] FIX_0_541196100 = 4433;  // cos(4π/16 + π/4)
localparam [13:0] FIX_1_306562965 = 10703; // cos(4π/16 - π/4)
localparam [13:0] FIX_1_847759065 = 15137; // cos(π/8)

//-----------------------------------------------------------------
// State Machine States
//-----------------------------------------------------------------
localparam STATE_IDLE     = 3'd0,
           STATE_COLLECT1 = 3'd1,
           STATE_COLLECT2 = 3'd2,
           STATE_COMPUTE  = 3'd3,
           STATE_OUTPUT   = 3'd4;

//-----------------------------------------------------------------
// Registers
//-----------------------------------------------------------------
reg [2:0]  state;              // Current state
reg [31:0] coef [0:7];         // Input coefficients
reg [31:0] output_values [0:7]; // Transformed outputs
reg [2:0]  col_idx;            // Column index (0-7)
reg [2:0]  output_count;       // Output cycle counter
reg [7:0]  valid_q;            // Output valid shift register
reg [5:0]  ptr_q;              // Output index pointer

//-----------------------------------------------------------------
// Input Scaling
//-----------------------------------------------------------------
wire signed [31:0] tmp_in [0:7];
genvar i;
generate
    for (i = 0; i < 8; i = i + 1) begin : gen_tmp_in
        assign tmp_in[i] = coef[i] << 8; // Scale inputs by 2^8
    end
endgenerate

//-----------------------------------------------------------------
// Even Part Computations
//-----------------------------------------------------------------
wire signed [31:0] z1_e = tmp_in[0] + tmp_in[4];
wire signed [31:0] z2_e = tmp_in[0] - tmp_in[4];
wire signed [31:0] z3_e = tmp_in[2] + tmp_in[6];
wire signed [31:0] z4_e = tmp_in[2] - tmp_in[6];

wire signed [31:0] tmp0_e = z1_e + z3_e;
wire signed [31:0] tmp2_e = z1_e - z3_e;

wire signed [45:0] mult0_e = z4_e * FIX_0_707106781;
// wire signed [31:0] shifted_mult0_e = (mult0_e >>> 13); // Truncate to 32 bits

wire signed [31:0] tmp1_e = z2_e + shifted_mult0_e;
wire signed [31:0] tmp3_e = z2_e - shifted_mult0_e;

//-----------------------------------------------------------------
// Odd Part Computations
//-----------------------------------------------------------------
wire signed [31:0] z1_o = tmp_in[1] + tmp_in[7];
wire signed [31:0] z2_o = tmp_in[3] + tmp_in[5];
wire signed [31:0] z3_o = tmp_in[1] + tmp_in[5];
wire signed [31:0] z4_o = tmp_in[3] + tmp_in[7];

// Sign-extend to 46 bits for multiplication
wire signed [45:0] z3_o_ext = {{14{z3_o[31]}}, z3_o};
wire signed [45:0] z4_o_ext = {{14{z4_o[31]}}, z4_o};

wire signed [45:0] mult1_o = (z3_o_ext - z4_o_ext) * FIX_1_847759065;
wire signed [45:0] mult2_o = z1_o * FIX_0_541196100;
wire signed [45:0] mult3_o = z2_o * FIX_1_306562965;
wire signed [45:0] mult4_o = (z3_o_ext + z4_o_ext) * FIX_0_382683433;

wire signed [45:0] shifted_mult0_e_temp = mult0_e >>> 13;
wire signed [45:0] shifted_mult1_o_temp = mult1_o >>> 13;
wire signed [45:0] shifted_mult2_o_temp = mult2_o >>> 13;
wire signed [45:0] shifted_mult3_o_temp = mult3_o >>> 13;
wire signed [45:0] shifted_mult4_o_temp = mult4_o >>> 13;

wire signed [31:0] shifted_mult0_e = shifted_mult0_e_temp[31:0];
wire signed [31:0] shifted_mult1_o = shifted_mult1_o_temp[31:0];
wire signed [31:0] shifted_mult2_o = shifted_mult2_o_temp[31:0];
wire signed [31:0] shifted_mult3_o = shifted_mult3_o_temp[31:0];
wire signed [31:0] shifted_mult4_o = shifted_mult4_o_temp[31:0];

wire signed [31:0] z5_o = shifted_mult1_o;

wire signed [31:0] tmp7_o = z1_o + z2_o;
wire signed [31:0] tmp6_o = shifted_mult2_o + z5_o;
wire signed [31:0] tmp5_o = shifted_mult3_o - z5_o;
wire signed [31:0] tmp4_o = shifted_mult4_o;

//-----------------------------------------------------------------
// Final Outputs
//-----------------------------------------------------------------
wire signed [31:0] out0 = (tmp0_e + tmp7_o + 32) >>> 6; // Rounding and scaling
wire signed [31:0] out1 = (tmp1_e + tmp6_o + 32) >>> 6;
wire signed [31:0] out2 = (tmp2_e + tmp5_o + 32) >>> 6;
wire signed [31:0] out3 = (tmp3_e + tmp4_o + 32) >>> 6;
wire signed [31:0] out4 = (tmp3_e - tmp4_o + 32) >>> 6;
wire signed [31:0] out5 = (tmp2_e - tmp5_o + 32) >>> 6;
wire signed [31:0] out6 = (tmp1_e - tmp6_o + 32) >>> 6;
wire signed [31:0] out7 = (tmp0_e - tmp7_o + 32) >>> 6;

//-----------------------------------------------------------------
// State Machine
//-----------------------------------------------------------------
always @(posedge clk_i) begin
    if (rst_i) begin
        state <= STATE_IDLE;
        valid_q <= 8'b0;
        ptr_q <= 6'd0;
        output_count <= 3'd0;
    end else begin
        case (state)
            STATE_IDLE: begin
                if (inport_valid_i) begin
                    coef[0] <= inport_data0_i;
                    coef[1] <= inport_data1_i;
                    coef[2] <= inport_data2_i;
                    coef[3] <= inport_data3_i;
                    col_idx <= inport_idx_i;
                    state <= STATE_COLLECT2;
                end
            end

            STATE_COLLECT2: begin
                if (inport_valid_i && inport_idx_i == col_idx) begin
                    coef[4] <= inport_data0_i;
                    coef[5] <= inport_data1_i;
                    coef[6] <= inport_data2_i;
                    coef[7] <= inport_data3_i;
                    state <= STATE_COMPUTE;
                end
            end

            STATE_COMPUTE: begin
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
                if (output_count < 7) begin
                    output_count <= output_count + 1;
                end else begin
                    state <= STATE_IDLE;
                end
                valid_q <= {valid_q[6:0], 1'b1};
                ptr_q <= ptr_q + 6'd1;
            end

            default: state <= STATE_IDLE;
        endcase

        if (img_start_i) begin
            ptr_q <= 6'd0;
            valid_q <= 8'b0;
        end
    end
end

//-----------------------------------------------------------------
// Output Assignments
//-----------------------------------------------------------------
assign outport_valid_o = valid_q[6];
assign outport_data_o  = output_values[ptr_q[2:0]];

//-----------------------------------------------------------------
// Pointer Conversion
//-----------------------------------------------------------------
function [5:0] ptr_conv;
    input [5:0] idx;
    reg [5:0] out_idx;
begin
    case (idx)
    6'd0:  out_idx = 6'd0;   6'd1:  out_idx = 6'd8;   6'd2:  out_idx = 6'd16;  6'd3:  out_idx = 6'd24;
    6'd4:  out_idx = 6'd32;  6'd5:  out_idx = 6'd40;  6'd6:  out_idx = 6'd48;  6'd7:  out_idx = 6'd56;
    6'd8:  out_idx = 6'd1;   6'd9:  out_idx = 6'd9;   6'd10: out_idx = 6'd17;  6'd11: out_idx = 6'd25;
    6'd12: out_idx = 6'd33;  6'd13: out_idx = 6'd41;  6'd14: out_idx = 6'd49;  6'd15: out_idx = 6'd57;
    6'd16: out_idx = 6'd2;   6'd17: out_idx = 6'd10;  6'd18: out_idx = 6'd18;  6'd19: out_idx = 6'd26;
    6'd20: out_idx = 6'd34;  6'd21: out_idx = 6'd42;  6'd22: out_idx = 6'd50;  6'd23: out_idx = 6'd58;
    6'd24: out_idx = 6'd3;   6'd25: out_idx = 6'd11;  6'd26: out_idx = 6'd19;  6'd27: out_idx = 6'd27;
    6'd28: out_idx = 6'd35;  6'd29: out_idx = 6'd43;  6'd30: out_idx = 6'd51;  6'd31: out_idx = 6'd59;
    6'd32: out_idx = 6'd4;   6'd33: out_idx = 6'd12;  6'd34: out_idx = 6'd20;  6'd35: out_idx = 6'd28;
    6'd36: out_idx = 6'd36;  6'd37: out_idx = 6'd44;  6'd38: out_idx = 6'd52;  6'd39: out_idx = 6'd60;
    6'd40: out_idx = 6'd5;   6'd41: out_idx = 6'd13;  6'd42: out_idx = 6'd21;  6'd43: out_idx = 6'd29;
    6'd44: out_idx = 6'd37;  6'd45: out_idx = 6'd45;  6'd46: out_idx = 6'd53;  6'd47: out_idx = 6'd61;
    6'd48: out_idx = 6'd6;   6'd49: out_idx = 6'd14;  6'd50: out_idx = 6'd22;  6'd51: out_idx = 6'd30;
    6'd52: out_idx = 6'd38;  6'd53: out_idx = 6'd46;  6'd54: out_idx = 6'd54;  6'd55: out_idx = 6'd62;
    6'd56: out_idx = 6'd7;   6'd57: out_idx = 6'd15;  6'd58: out_idx = 6'd23;  6'd59: out_idx = 6'd31;
    6'd60: out_idx = 6'd39;  6'd61: out_idx = 6'd47;  6'd62: out_idx = 6'd55;  default: out_idx = 6'd63;
    endcase
    ptr_conv = out_idx;
end
endfunction

assign outport_idx_o = ptr_conv(ptr_q);

//State machine debug
always @(posedge clk_i) begin
    $display("Time: %t, State_y: %d, rst_i: %b, inport_valid_i: %b, inport_idx_i: %d, col_idx: %d",
             $time, state, rst_i, inport_valid_i, inport_idx_i, col_idx);
end

endmodule