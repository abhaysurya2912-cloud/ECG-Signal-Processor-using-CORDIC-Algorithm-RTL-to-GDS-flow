# ECG Signal Processor using CORDIC Algorithm — RTL to GDS Flow

> A complete digital VLSI implementation of an ECG signal processing core using the CORDIC (Coordinate Rotation Digital Computer) algorithm, taken through a full RTL-to-GDS flow using open-source and/or industry-standard EDA tools.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Block Diagram](#block-diagram)
3. [CORDIC Algorithm — Theory](#cordic-algorithm--theory)
4. [ECG Signal Processing Pipeline](#ecg-signal-processing-pipeline)
5. [ECG Signal Extraction from MIT-BIH Database](#ECG-Signal-Extraction-from-MIT-BIH-Database)
6. [RTL Design](#rtl-design)
   - [Top-Level: qrs_cordic_detector](#top-level-qrs_cordic_detector)
   - [Bandpass Filter (5-15 Hz)](#bandpass-filter-515-hz)
   - [Pipelined CORDIC Magnitude](#pipelined-cordic-magnitude)
   - [QRS R-Peak Detector](#qrs-r-peak-detector)
   - [RR Interval and BPM Calculator](#rr-interval-and-bpm-calculator)
   - [Integer Square Root](#integer-square-root)
   - [HRV Metrics Module](#hrv-metrics-module)
7. [Testbench](#testbench)
8. [Functional Simulation](#functional-simulation)
9. [Synthesis](#synthesis)
10. [Placement](#Placement)
11. [CTS](#CTS)
12. [Routing](#Routing)
13. [DRC-Connectivity-Verification](#DRC-Connectivity-verification)
14. [Sign-off and GDS Generation](#sign-off-and-gds-generation)
15. [Results and Reports](#results-and-reports)
16. [References](#references)

---

## Project Overview

This project implements a **hardware ECG signal processor** that performs real-time computation of key cardiac parameters — including QRS complex detection, R-peak identification, heart rate estimation, and Heart Rate Variability (HRV) metrics — using the **CORDIC algorithm** as the computational backbone for hardware-efficient vector magnitude operations (no multipliers required).

The design is taken through a complete **RTL-to-GDS (tape-out-ready) physical design flow**, covering RTL coding, functional simulation, logic synthesis, DFT scan insertion, static timing analysis, formal verification, place and route, physical verification, and GDS-II generation.

### Key Design Parameters

| Parameter | Value |
|-----------|-------|
| Technology node | [e.g., 45 nm / 28 nm / Sky130 PDK] |
| Sampling frequency | 360 Hz (MIT-BIH standard) |
| Clock frequency | [e.g., 100 MHz] |
| CORDIC iterations | 16 |
| Pipeline stages | 16 (fully pipelined CORDIC) |
| Word length | 16-bit signed fixed-point |
| Bandpass filter | 5–15 Hz, 5-tap symmetric FIR |
| HRV metrics | RMSSD, SDNN, pNN50, Pearson Correlation |
| Scan coverage target | ≥ 95% |
| Core area | [fill after PnR] |
| Power (typical) | [fill after PnR] |

---

## Block Diagram

<img src="WhatsApp Image 2026-04-19 at 11.27.33 PM.jpeg" width="300">

---

## CORDIC Algorithm — Theory

<img src="Screenshot 2026-04-20 183536.png" width="300">

The **CORDIC (Coordinate Rotation Digital Computer)** algorithm computes trigonometric, hyperbolic, and vector magnitude functions using only **additions, subtractions, and bit shifts** — making it ideal for resource-constrained hardware with no dedicated multiplier units.

### Vector Mode — used in this design

At each iteration `i`, the micro-rotation rule is:

```
if y[i] >= 0:
    x[i+1] = x[i] + (y[i] >> i)
    y[i+1] = y[i] - (x[i] >> i)
else:
    x[i+1] = x[i] - (y[i] >> i)
    y[i+1] = y[i] + (x[i] >> i)
```

After `N` iterations, `x[N]` converges to `K * sqrt(x^2 + y^2)` where `K ≈ 1.6468` is the CORDIC gain constant. This design takes:
- `x_in` = current filtered ECG sample
- `y_in` = differentiated ECG sample (captures slope/edge energy)

The magnitude output captures the combined signal energy envelope, enabling robust R-peak detection that is resistant to amplitude drift.

### Implementation Parameters

| Parameter | Value |
|-----------|-------|
| Iterations | 16 |
| Data format | 16-bit signed fixed-point |
| Mode | Vector (magnitude computation) |
| Architecture | Fully pipelined (one result per clock after 16-cycle latency) |

> **CORDIC Gain Note:** Output is scaled by K ≈ 1.6468. Thresholds in `qrs_peak_detect` are calibrated against raw CORDIC output. For absolute magnitude, multiply output by `1/K ≈ 0.6073` (fixed-point: multiply by 39243, right-shift by 16).

---

## ECG Signal Processing Pipeline

### Stage 0 — Bandpass Filter (5–15 Hz)
5-tap symmetric FIR filter targeting the QRS energy band. Removes baseline wander (< 5 Hz) and high-frequency EMG noise (> 15 Hz). Coefficients: `[119, 206, 226, 206, 119]`, normalized by 876.

### Stage 1 — Differentiator
Computes `ecg_diff = ecg_filtered[n] − ecg_filtered[n−1]` to highlight the steep slope characteristic of a QRS complex.

### Stage 2 — Pipelined CORDIC Magnitude
Computes `mag ≈ K * sqrt(ecg_filtered^2 + ecg_diff^2)` using the 16-stage pipelined CORDIC core.

### Stage 3 — QRS R-Peak Detector
Dual-threshold hysteresis detector (HIGH = 6000, LOW = 4000) with a 140-sample refractory period (~389 ms @ 360 Hz) to suppress double-detection within a single beat.

### Stage 4 — RR Interval and BPM
Counts sample ticks between consecutive R-peaks. `RR_ms = (count * 1000) / FS_HZ`, `BPM = (FS_HZ * 60) / count`.

### Stage 5 — HRV Metrics
Computes RMSSD, SDNN, pNN50, and Pearson Correlation Coefficient between consecutive RR pairs using running accumulators. Valid RR range: 300–2000 ms.

### Disease Classification (Testbench)

| Classification | Criterion |
|----------------|-----------|
| BRADYCARDIA | Average BPM < 60 |
| TACHYCARDIA | Average BPM > 100 |
| ARRHYTHMIA_SUSPECTED | SDNN > 100 ms OR RMSSD > 80 ms OR pNN50 > 20% |
| NORMAL | None of the above |

---


---

## RTL Design

### Top-Level: qrs_cordic_detector

The top-level module instantiates and wires all sub-modules in the signal processing pipeline.

```verilog
//============================================================
// QRS complex detection using CORDIC magnitude (pipelined)
// + HRV metrics: RMSSD, SDNN, pNN50, Correlation Coefficient
// + Bandpass filter (5-15 Hz) for noise reduction
//============================================================
module qrs_cordic_detector #(
    parameter FS_HZ          = 360,
    parameter DATA_WIDTH     = 16,
    parameter CORDIC_ITER    = 16,
    parameter THRESHOLD_MAG  = 8000
)(
    input  wire                         clk,
    input  wire                         reset_n,
    input  wire signed [DATA_WIDTH-1:0] ecg_sample,
    input  wire                         sample_valid,
    output wire                         r_peak_pulse,
    output wire [31:0]                  rr_interval_ms,
    output wire [15:0]                  bpm,
    output wire                         rr_valid,
    output wire [31:0]                  rmssd_out,
    output wire [31:0]                  sdnn_out,
    output wire [15:0]                  pnn50_out,
    output wire signed [15:0]           correlation_coeff,
    output wire                         hrv_valid
);

    // -----------------------------
    // 0) Bandpass filter (5-15 Hz)
    // -----------------------------
    wire signed [DATA_WIDTH-1:0] ecg_filtered;
    wire                         filter_valid;

    bandpass_filter_5_15 u_bpf (
        .clk       (clk),
        .reset_n   (reset_n),
        .ecg_in    (ecg_sample),
        .in_valid  (sample_valid),
        .ecg_out   (ecg_filtered),
        .out_valid (filter_valid)
    );

    // -----------------------------
    // 1) Simple differentiator
    // -----------------------------
    reg signed [DATA_WIDTH-1:0] ecg_prev;
    wire signed [DATA_WIDTH:0]   ecg_diff;

    always @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            ecg_prev <= 'd0;
        end else if (filter_valid) begin
            ecg_prev <= ecg_filtered;
        end
    end

    assign ecg_diff = ecg_filtered - ecg_prev;

    // -----------------------------
    // 2) Pipelined CORDIC magnitude
    // -----------------------------
    localparam CORDIC_W = CORDIC_ITER;

    wire signed [CORDIC_W-1:0] cordic_x_in;
    wire signed [CORDIC_W-1:0] cordic_y_in;

    assign cordic_x_in = ecg_filtered[DATA_WIDTH-1 -: CORDIC_W];
    assign cordic_y_in = ecg_diff    [DATA_WIDTH   -: CORDIC_W];

    wire [CORDIC_W-1:0] mag_out;
    wire                mag_valid;

    cordic_mag_pipelined #(
        .WIDTH (CORDIC_W),
        .ITER  (CORDIC_ITER)
    ) u_cordic (
        .clk       (clk),
        .reset_n   (reset_n),
        .x_in      (cordic_x_in),
        .y_in      (cordic_y_in),
        .in_valid  (filter_valid),
        .mag_out   (mag_out),
        .mag_valid (mag_valid)
    );

    // -----------------------------
    // 3) QRS / R-peak detector
    // -----------------------------
    qrs_peak_detect #(
        .WIDTH          (CORDIC_W),
        .THRESHOLD_HIGH (6000),
        .THRESHOLD_LOW  (4000),
        .REFRAC_SAMPLES (140)
    ) u_peak (
        .clk           (clk),
        .reset_n       (reset_n),
        .mag_sample    (mag_out),
        .mag_valid     (mag_valid),
        .r_peak_pulse  (r_peak_pulse)
    );

    // -----------------------------
    // 4) RR and BPM
    // -----------------------------
    rr_bpm_calc #(
        .FS_HZ (FS_HZ)
    ) u_rr (
        .clk            (clk),
        .reset_n        (reset_n),
        .sample_tick    (sample_valid),
        .r_peak_pulse   (r_peak_pulse),
        .rr_interval_ms (rr_interval_ms),
        .bpm            (bpm),
        .rr_valid       (rr_valid)
    );

    // -----------------------------
    // 5) HRV metrics (with correlation)
    // -----------------------------
    hrv_metrics u_hrv (
        .clk                (clk),
        .reset_n            (reset_n),
        .rr_interval_ms     (rr_interval_ms),
        .rr_valid           (rr_valid),
        .rmssd_out          (rmssd_out),
        .sdnn_out           (sdnn_out),
        .pnn50_out          (pnn50_out),
        .correlation_coeff  (correlation_coeff),
        .hrv_valid          (hrv_valid)
    );

endmodule
```

---

### Bandpass Filter (5–15 Hz)

5-tap symmetric FIR filter. Suppresses baseline wander and high-frequency noise outside the QRS band.

```verilog
//============================================================
// Bandpass filter (5-15 Hz @ 360 Hz)
//============================================================
module bandpass_filter_5_15 (
    input  wire        clk,
    input  wire        reset_n,
    input  wire signed [15:0] ecg_in,
    input  wire        in_valid,
    output reg  signed [15:0] ecg_out,
    output reg         out_valid
);
    localparam COEFF0 = 119;
    localparam COEFF1 = 206;
    localparam COEFF2 = 226;

    reg signed [15:0] delay [0:4];
    reg signed [31:0] sum;
    integer i;

    always @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            for (i = 0; i < 5; i = i + 1) delay[i] <= 16'sd0;
            out_valid <= 1'b0;
            ecg_out   <= 16'sd0;
        end else begin
            out_valid <= in_valid;

            if (in_valid) begin
                for (i = 4; i > 0; i = i - 1)
                    delay[i] <= delay[i-1];
                delay[0] <= ecg_in;

                sum = COEFF0 * (delay[0] + delay[4])
                    + COEFF1 * (delay[1] + delay[3])
                    + COEFF2 * (delay[2]);

                ecg_out <= sum / 876;
            end
        end
    end
endmodule
```

---

### Pipelined CORDIC Magnitude

16-stage fully pipelined CORDIC core in vector mode. Each stage implements one micro-rotation. Latency = ITER clock cycles; throughput = 1 result per clock after pipeline fill.

```verilog
//============================================================
// Pipelined CORDIC magnitude
//============================================================
module cordic_mag_pipelined #(
    parameter WIDTH = 16,
    parameter ITER  = 16
)(
    input  wire                    clk,
    input  wire                    reset_n,
    input  wire signed [WIDTH-1:0] x_in,
    input  wire signed [WIDTH-1:0] y_in,
    input  wire                    in_valid,
    output wire [WIDTH-1:0]        mag_out,
    output wire                    mag_valid
);

    reg signed [WIDTH-1:0] x [0:ITER];
    reg signed [WIDTH-1:0] y [0:ITER];
    reg                    v [0:ITER];
    integer i;

    always @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            for (i = 0; i <= ITER; i = i + 1) begin
                x[i] <= 'd0;
                y[i] <= 'd0;
                v[i] <= 1'b0;
            end
        end else begin
            x[0] <= x_in;
            y[0] <= y_in;
            v[0] <= in_valid;

            for (i = 0; i < ITER; i = i + 1) begin
                if (y[i] >= 0) begin
                    x[i+1] <= x[i] + (y[i] >>> i);
                    y[i+1] <= y[i] - (x[i] >>> i);
                end else begin
                    x[i+1] <= x[i] - (y[i] >>> i);
                    y[i+1] <= y[i] + (x[i] >>> i);
                end
                v[i+1] <= v[i];
            end
        end
    end

    assign mag_out   = x[ITER];
    assign mag_valid = v[ITER];

endmodule
```

---

### QRS R-Peak Detector

Dual-threshold hysteresis detector with post-detection refractory period.

| Parameter | Default | Description |
|-----------|---------|-------------|
| THRESHOLD_HIGH | 6000 | Trigger level — must exceed to fire |
| THRESHOLD_LOW | 4000 | Reset level — must fall below to re-arm |
| REFRAC_SAMPLES | 140 | ~389 ms blackout after each detection |

```verilog
//============================================================
// QRS R-peak detector
//============================================================
module qrs_peak_detect #(
    parameter WIDTH           = 16,
    parameter THRESHOLD_HIGH  = 6000,
    parameter THRESHOLD_LOW   = 4000,
    parameter REFRAC_SAMPLES  = 180
)(
    input  wire             clk,
    input  wire             reset_n,
    input  wire [WIDTH-1:0] mag_sample,
    input  wire             mag_valid,
    output reg              r_peak_pulse
);

    reg armed;
    reg [$clog2(REFRAC_SAMPLES+1)-1:0] refrac_cnt;

    wire above_high = (mag_sample >= THRESHOLD_HIGH);
    wire below_low  = (mag_sample <= THRESHOLD_LOW);

    always @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            armed        <= 1'b1;
            refrac_cnt   <= 'd0;
            r_peak_pulse <= 1'b0;
        end else begin
            r_peak_pulse <= 1'b0;

            if (mag_valid) begin
                if (refrac_cnt != 0) begin
                    refrac_cnt <= refrac_cnt - 1;
                    armed      <= 1'b0;
                end else begin
                    if (armed && above_high) begin
                        r_peak_pulse <= 1'b1;
                        armed        <= 1'b0;
                        refrac_cnt   <= REFRAC_SAMPLES - 1;
                    end else if (!armed && below_low) begin
                        armed <= 1'b1;
                    end
                end
            end
        end
    end

endmodule
```

---

### RR Interval and BPM Calculator

Counts sample ticks between consecutive R-peaks and computes the RR interval (ms) and instantaneous heart rate (BPM).

```verilog
//============================================================
// RR interval and BPM calculator
//============================================================
module rr_bpm_calc #(
    parameter FS_HZ = 360
)(
    input  wire        clk,
    input  wire        reset_n,
    input  wire        sample_tick,
    input  wire        r_peak_pulse,
    output reg  [31:0] rr_interval_ms,
    output reg  [15:0] bpm,
    output reg         rr_valid
);

    reg [31:0] sample_count;
    reg        first_peak_seen;

    always @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            sample_count    <= 32'd0;
            first_peak_seen <= 1'b0;
            rr_interval_ms  <= 32'd0;
            bpm             <= 16'd0;
            rr_valid        <= 1'b0;
        end else begin
            rr_valid <= 1'b0;

            if (sample_tick)
                sample_count <= sample_count + 1;

            if (r_peak_pulse) begin
                if (!first_peak_seen) begin
                    first_peak_seen <= 1'b1;
                    sample_count    <= 32'd0;
                end else begin
                    rr_interval_ms <= (sample_count * 32'd1000) / FS_HZ;

                    if (sample_count != 0)
                        bpm <= (FS_HZ * 32'd60) / sample_count;
                    else
                        bpm <= 16'd0;

                    rr_valid     <= 1'b1;
                    sample_count <= 32'd0;
                end
            end
        end
    end

endmodule
```

---

### Integer Square Root

Combinational 32-bit integer square root, producing a 16-bit result. Used by `hrv_metrics` for RMSSD, SDNN, and correlation denominator computation.

```verilog
//============================================================
// Integer square root (32-bit input -> 16-bit output)
//============================================================
module isqrt32 (
    input  wire [31:0] x,
    output wire [15:0] y
);
    reg [31:0] rem;
    reg [15:0] root;
    reg [31:0] test;
    integer i;

    always @(*) begin
        rem  = 32'd0;
        root = 16'd0;
        for (i = 15; i >= 0; i = i - 1) begin
            rem  = {rem[29:0], x[2*i+1], x[2*i]};
            test = {root, 2'b01};
            if (rem >= {root, 1'b0, 1'b0} + 1) begin
                rem  = rem - ({root, 1'b0} + 1);
                root = {root[14:0], 1'b1};
            end else begin
                root = {root[14:0], 1'b0};
            end
        end
    end

    assign y = root;
endmodule
```

> **Synthesis note:** This is a 16-iteration combinational loop — it will unroll fully and may become a timing-critical path. Registering its inputs/outputs or pipelining may be needed after STA.

---

### HRV Metrics Module

Computes four HRV metrics using running accumulators updated on each valid RR interval. The Pearson Correlation Coefficient is computed between consecutive RR pairs, scaled by 1000 and represented as a signed 16-bit integer (range −1000 to +1000). Outputs are registered and updated after each new beat.

```verilog
//============================================================
// HRV Metrics Module (with CORRECTED Pearson Correlation)
//============================================================
module hrv_metrics (
    input  wire        clk,
    input  wire        reset_n,
    input  wire [31:0] rr_interval_ms,
    input  wire        rr_valid,
    output reg  [31:0] rmssd_out,
    output reg  [31:0] sdnn_out,
    output reg  [15:0] pnn50_out,
    output reg  signed [15:0] correlation_coeff,
    output reg         hrv_valid
);

    localparam RR_MIN_MS = 300;
    localparam RR_MAX_MS = 2000;

    // Basic statistics
    reg [31:0]  N;
    reg [63:0]  sum_rr;
    reg [63:0]  sum_rr2;
    reg [31:0]  prev_rr;
    reg         prev_valid;
    reg [63:0]  sum_drr2;
    reg [31:0]  N_diff;
    reg [31:0]  cnt_nn50;

    // Correlation coefficient accumulators
    reg [63:0]  sum_xy;
    reg [63:0]  sum_x;
    reg [63:0]  sum_y;
    reg [63:0]  sum_x2;
    reg [63:0]  sum_y2;
    reg [31:0]  N_corr;

    wire in_range = (rr_interval_ms >= RR_MIN_MS) &&
                    (rr_interval_ms <= RR_MAX_MS);

    wire signed [32:0] d_rr;
    assign d_rr = $signed({1'b0, rr_interval_ms}) -
                  $signed({1'b0, prev_rr});

    wire [32:0] abs_drr;
    assign abs_drr = d_rr[32] ? (~d_rr + 1'b1) : d_rr;

    wire [63:0] drr2;
    assign drr2 = abs_drr * abs_drr;

    wire [63:0] rr2;
    assign rr2 = rr_interval_ms * rr_interval_ms;

    always @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            N          <= 32'd0;
            sum_rr     <= 64'd0;
            sum_rr2    <= 64'd0;
            prev_rr    <= 32'd0;
            prev_valid <= 1'b0;
            sum_drr2   <= 64'd0;
            N_diff     <= 32'd0;
            cnt_nn50   <= 32'd0;

            sum_xy     <= 64'd0;
            sum_x      <= 64'd0;
            sum_y      <= 64'd0;
            sum_x2     <= 64'd0;
            sum_y2     <= 64'd0;
            N_corr     <= 32'd0;

            rmssd_out  <= 32'd0;
            sdnn_out   <= 32'd0;
            pnn50_out  <= 16'd0;
            correlation_coeff <= 16'sd0;
            hrv_valid  <= 1'b0;
        end else begin
            hrv_valid <= 1'b0;

            if (rr_valid && in_range) begin
                N       <= N + 1;
                sum_rr  <= sum_rr  + rr_interval_ms;
                sum_rr2 <= sum_rr2 + rr2;

                if (prev_valid) begin
                    sum_drr2 <= sum_drr2 + drr2;
                    N_diff   <= N_diff + 1;
                    if (abs_drr > 33'd50)
                        cnt_nn50 <= cnt_nn50 + 1;
                end

                if (prev_valid) begin
                    sum_xy <= sum_xy + (rr_interval_ms * prev_rr);
                    sum_x  <= sum_x  + rr_interval_ms;
                    sum_y  <= sum_y  + prev_rr;
                    sum_x2 <= sum_x2 + (rr_interval_ms * rr_interval_ms);
                    sum_y2 <= sum_y2 + (prev_rr * prev_rr);
                    N_corr <= N_corr + 1;
                end

                prev_rr    <= rr_interval_ms;
                prev_valid <= 1'b1;
                hrv_valid  <= 1'b1;
            end
        end
    end

    // --- RMSSD ---
    wire [63:0] mean_drr2 = (N_diff > 0) ? (sum_drr2 / N_diff) : 64'd0;
    wire [15:0] rmssd_tmp;
    isqrt32 u_sqrt_rmssd (.x(mean_drr2[31:0]), .y(rmssd_tmp));
    wire [31:0] rmssd_wire = {16'd0, rmssd_tmp};

    // --- SDNN ---
    wire [63:0] mean_rr    = (N > 0) ? (sum_rr  / N) : 64'd0;
    wire [63:0] mean_rr2   = (N > 0) ? (sum_rr2 / N) : 64'd0;
    wire [63:0] mean_rr_sq = mean_rr * mean_rr;
    wire [63:0] var_rr     = (mean_rr2 >= mean_rr_sq) ?
                             (mean_rr2 - mean_rr_sq) : 64'd0;
    wire [15:0] sdnn_tmp;
    isqrt32 u_sqrt_sdnn (.x(var_rr[31:0]), .y(sdnn_tmp));
    wire [31:0] sdnn_wire = {16'd0, sdnn_tmp};

    // --- pNN50 ---
    wire [31:0] pnn50_wire = (N_diff > 0) ?
                             ((cnt_nn50 * 32'd100) / N_diff) : 32'd0;

    // --- Pearson Correlation Coefficient ---
    wire [63:0] numerator_signed    = (N_corr * sum_xy) - (sum_x * sum_y);
    wire        numerator_negative  = numerator_signed[63];
    wire [63:0] numerator_abs       = numerator_negative ?
                                      (~numerator_signed + 1'b1) : numerator_signed;
    wire [63:0] denominator_x       = (N_corr * sum_x2) - (sum_x * sum_x);
    wire [63:0] denominator_y       = (N_corr * sum_y2) - (sum_y * sum_y);
    wire [63:0] denominator_product = denominator_x * denominator_y;
    wire [15:0] denominator_sqrt16;
    isqrt32 u_sqrt_corr (.x(denominator_product[31:0]), .y(denominator_sqrt16));
    wire [31:0] denominator_sqrt_full    = {16'd0, denominator_sqrt16};
    wire [31:0] correlation_scaled_abs   =
                (denominator_sqrt_full > 0) ?
                ((numerator_abs * 32'd1000) / denominator_sqrt_full) : 32'd0;
    wire [31:0] correlation_scaled       =
                (correlation_scaled_abs > 32'd1000) ? 32'd1000 : correlation_scaled_abs;
    wire signed [15:0] correlation_final =
                numerator_negative ?
                -$signed(correlation_scaled[15:0]) :
                 $signed(correlation_scaled[15:0]);

    // --- Register outputs ---
    always @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            rmssd_out         <= 32'd0;
            sdnn_out          <= 32'd0;
            pnn50_out         <= 16'd0;
            correlation_coeff <= 16'sd0;
        end else if (hrv_valid) begin
            rmssd_out <= rmssd_wire;
            sdnn_out  <= sdnn_wire;
            pnn50_out <= pnn50_wire[15:0];
            // Only output correlation if we have enough samples (>= 3 pairs)
            if (N_corr >= 3)
                correlation_coeff <= correlation_final;
            else
                correlation_coeff <= 16'sd0;
        end
    end

endmodule
```

---

## Testbench

### Generating the Input File

```python
# tb/ecg_samples/generate_ecg_txt.py
import wfdb

record = wfdb.rdrecord('100', sampto=3600, pn_dir='mitdb')
samples = record.p_signal[:, 0]          # Channel 0 (MLII lead)
adc = (samples * 200).astype(int)        # Scale to integer range
with open('100_ecg.txt', 'w') as f:
    for s in adc:
        f.write(f'{s}\n')
```

Run once before simulation:
```bash
pip install wfdb numpy
python3 tb/ecg_samples/generate_ecg_txt.py
cp tb/ecg_samples/100_ecg.txt sim/
```

### Output Files

| File | Content |
|------|---------|
| `detection_samples.txt` | Sample index of each detected R-peak (one per line) |
| `final_summary.txt` | CSV-style key-value metrics summary |

### Testbench Code

```systemverilog
`timescale 1ns/1ps

module tb_qrs_final_display;

    // -------------------------
    // Parameters
    // -------------------------
    localparam int FS_HZ      = 360;
    localparam int DATA_WIDTH = 16;
    localparam time CLK_PERIOD = 10ns;

    // -------------------------
    // DUT I/O signals
    // -------------------------
    logic clk;
    logic reset_n;

    logic signed [DATA_WIDTH-1:0] ecg_sample;
    logic                         sample_valid;

    logic                         r_peak_pulse;
    logic [31:0]                  rr_interval_ms;
    logic [15:0]                  bpm;
    logic                         rr_valid;

    // ---- HRV outputs ----
    logic [31:0]                  rmssd_out;
    logic [31:0]                  sdnn_out;
    logic [15:0]                  pnn50_out;
    logic signed [15:0]           correlation_coeff;
    logic                         hrv_valid;

    // -------------------------
    // TB accumulators
    // -------------------------
    logic [63:0] sum_rr_ms;
    int          rr_count;
    int          rr_min_ms;
    int          rr_max_ms;
    int          r_peak_count;
    int          sample_idx_save;

    // final results
    int    final_mean_rr_ms;
    int    final_avg_bpm;
    real   final_correlation;
    string final_disease_str;

    // -------------------------
    // Files & indices
    // -------------------------
    integer fd_in;
    integer fd_det;
    integer fd_summary;
    integer r;
    integer sample_int;
    integer sample_idx;

    // -------------------------
    // Clock
    // -------------------------
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk;
    end

    // -------------------------
    // Instantiate DUT
    // -------------------------
    qrs_cordic_detector #(
        .FS_HZ      (FS_HZ),
        .DATA_WIDTH (DATA_WIDTH),
        .CORDIC_ITER(16)
    ) dut (
        .clk                (clk),
        .reset_n            (reset_n),
        .ecg_sample         (ecg_sample),
        .sample_valid       (sample_valid),
        .r_peak_pulse       (r_peak_pulse),
        .rr_interval_ms     (rr_interval_ms),
        .bpm                (bpm),
        .rr_valid           (rr_valid),
        .rmssd_out          (rmssd_out),
        .sdnn_out           (sdnn_out),
        .pnn50_out          (pnn50_out),
        .correlation_coeff  (correlation_coeff),
        .hrv_valid          (hrv_valid)
    );

    // -------------------------
    // Stimulus: feed ECG file
    // -------------------------
    initial begin
        sample_valid  = 1'b0;
        sample_idx    = -1;

        reset_n = 1'b0;
        repeat (10) @(posedge clk);
        reset_n = 1'b1;

        fd_in = $fopen("100_ecg.txt", "r");
        if (fd_in == 0) begin
            $display("ERROR: Could not open 100_ecg.txt");
            $finish;
        end else begin
            $display("Opened 100_ecg.txt for reading.");
        end

        fd_det     = $fopen("detection_samples.txt", "w");
        fd_summary = $fopen("final_summary.txt",     "w");

        while (!$feof(fd_in)) begin
            r = $fscanf(fd_in, "%d\n", sample_int);
            if (r == 1) begin
                sample_idx = sample_idx + 1;

                @(posedge clk);
                ecg_sample   <= sample_int[DATA_WIDTH-1:0];
                sample_valid <= 1'b1;

                @(posedge clk);
                sample_valid <= 1'b0;
            end else begin
                @(posedge clk);
            end
        end

        $display("EOF reached. Waiting for pipeline to flush...");
        repeat (1000) @(posedge clk);

        // --------------------------------------------------
        // Compute final metrics
        // --------------------------------------------------
        if (rr_count > 0) begin
            final_mean_rr_ms = sum_rr_ms / rr_count;
            if (final_mean_rr_ms > 0)
                final_avg_bpm = 60000 / final_mean_rr_ms;
            else
                final_avg_bpm = 0;

            final_correlation = correlation_coeff / 1000.0;

            if (final_avg_bpm < 60)
                final_disease_str = "BRADYCARDIA";
            else if (final_avg_bpm > 100)
                final_disease_str = "TACHYCARDIA";
            else if (sdnn_out > 100 || rmssd_out > 80 || pnn50_out > 20)
                final_disease_str = "ARRHYTHMIA_SUSPECTED";
            else
                final_disease_str = "NORMAL";
        end else begin
            final_mean_rr_ms  = 0;
            final_avg_bpm     = 0;
            final_correlation = 0.0;
            final_disease_str = "NO_RR_DETECTIONS";
        end

        $display("\n===== FINAL SUMMARY =====");
        $display("R-peaks detected       : %0d",  r_peak_count);
        $display("RR intervals collected : %0d",  rr_count);
        $display("Mean RR (ms)           : %0d",  final_mean_rr_ms);
        $display("Average HR (bpm)       : %0d",  final_avg_bpm);
        $display("RR min..max (ms)       : %0d .. %0d", rr_min_ms, rr_max_ms);
        $display("--- HRV Metrics ---");
        $display("RMSSD (ms)             : %0d ms",  rmssd_out);
        $display("SDNN  (ms)             : %0d ms",  sdnn_out);
        $display("pNN50 (%%)              : %0d %%",  pnn50_out);
        $display("Correlation Coeff (r)  : %0.3f",   final_correlation);
        $display("Disease classification : %s",      final_disease_str);
        $display("=========================\n");

        $fdisplay(fd_summary, "R_PEAK_COUNT,%0d",   r_peak_count);
        $fdisplay(fd_summary, "RR_COUNT,%0d",        rr_count);
        $fdisplay(fd_summary, "MEAN_RR_MS,%0d",      final_mean_rr_ms);
        $fdisplay(fd_summary, "AVG_BPM,%0d",         final_avg_bpm);
        $fdisplay(fd_summary, "RR_MIN_MS,%0d",       rr_min_ms);
        $fdisplay(fd_summary, "RR_MAX_MS,%0d",       rr_max_ms);
        $fdisplay(fd_summary, "RMSSD_MS,%0d",        rmssd_out);
        $fdisplay(fd_summary, "SDNN_MS,%0d",         sdnn_out);
        $fdisplay(fd_summary, "PNN50_PCT,%0d",       pnn50_out);
        $fdisplay(fd_summary, "CORRELATION_R,%0.3f", final_correlation);
        $fdisplay(fd_summary, "DISEASE,%s",          final_disease_str);

        $fclose(fd_in);
        $fclose(fd_det);
        $fclose(fd_summary);
        $finish;
    end

    // -------------------------
    // Monitor & accumulators
    // -------------------------
    always_ff @(posedge clk or negedge reset_n) begin
        if (!reset_n) begin
            sum_rr_ms       <= 64'd0;
            rr_count        <= 0;
            rr_min_ms       <= 32'h7fffffff;
            rr_max_ms       <= 0;
            r_peak_count    <= 0;
            sample_idx_save <= 0;
        end else begin
            if (r_peak_pulse) begin
                r_peak_count    <= r_peak_count + 1;
                sample_idx_save <= sample_idx;
                $fdisplay(fd_det, "%0d", sample_idx);
                $display("[%0t] R-peak #%0d detected at sample %0d",
                         $time, r_peak_count + 1, sample_idx);
            end

            if (rr_valid) begin
                sum_rr_ms <= sum_rr_ms + rr_interval_ms;
                rr_count  <= rr_count  + 1;
                if (rr_interval_ms < rr_min_ms) rr_min_ms <= rr_interval_ms;
                if (rr_interval_ms > rr_max_ms) rr_max_ms <= rr_interval_ms;
                $display("[%0t] RR=%0d ms  BPM=%0d  count=%0d",
                         $time, rr_interval_ms, bpm, rr_count);
            end

            if (hrv_valid) begin
                $display("[%0t] HRV update -> RMSSD=%0d ms  SDNN=%0d ms  pNN50=%0d%%  r=%0.3f",
                         $time, rmssd_out, sdnn_out, pnn50_out,
                         correlation_coeff / 1000.0);
            end
        end
    end

endmodule
```

---

## Functional Simulation

### Running the Simulation

```bash
cd sim/

# Cadence Xcelium (recommended — full SystemVerilog support)
xrun -sv \
  ../rtl/qrs_cordic_detector.v \
  ../rtl/bandpass_filter_5_15.v \
  ../rtl/cordic_mag_pipelined.v \
  ../rtl/qrs_peak_detect.v \
  ../rtl/rr_bpm_calc.v \
  ../rtl/isqrt32.v \
  ../rtl/hrv_metrics.v \
  ../tb/tb_qrs_final_display.sv

# Synopsys VCS
vcs -full64 -sverilog ../rtl/*.v ../tb/tb_qrs_final_display.sv -o simv && ./simv

# ModelSim / Questa
vlog -sv ../rtl/*.v ../tb/tb_qrs_final_display.sv
vsim -c tb_qrs_final_display -do "run -all; quit"
```

> **Icarus Verilog note:** The testbench uses `logic`, `always_ff`, `string`, `real`, and `int` — compile with `iverilog -g2012`. `string` type support is limited; use VCS or Xcelium for full compatibility.

### Simulation Checklist

- [ ] CORDIC output verified against Python reference: `K * np.sqrt(x**2 + y**2)`
- [ ] R-peak count matches MIT-BIH Record 100 annotations (expect ~75 BPM, ~32 beats per 10 s)
- [ ] RMSSD and SDNN are physiologically plausible
- [ ] Correlation coefficient is within [−1.0, +1.0]
- [ ] No X-propagation on any output
- [ ] `detection_samples.txt` and `final_summary.txt` generated correctly

---

## Synthesis

### Tool: Synopsys Design Compiler

```bash
cd synthesis/
dc_shell -f scripts/dc_synthesis.tcl | tee logs/dc_synthesis.log
```

### Key SDC Constraints

```tcl
# constraints/ecg_top.sdc
create_clock -name clk -period 10.0 [get_ports clk]
set_input_delay  2.0 -clock clk [all_inputs]
set_output_delay 2.0 -clock clk [all_outputs]
set_max_fanout   16  [current_design]
set_max_transition 0.5 [current_design]
```

### Synthesis Results

| Metric | Value |
|--------|-------|
| Total cell area | [fill after synthesis] |
| Combinational area | [fill] |
| Sequential area | [fill] |
| WNS | [fill] ns |
| TNS | [fill] ns |

---

## Design for Testability (DFT)

### Tool: Synopsys DFT Compiler

```bash
cd dft/
dc_shell -f scripts/dft_compiler.tcl | tee logs/dft.log
```

### Scan Configuration Script

```tcl
set_dft_signal -view existing_dft -type ScanClock  -port clk -timing {45 55}
set_dft_signal -view spec         -type Reset       -port reset_n -active_state 0
set_dft_signal -view spec         -type ScanEnable  -port scan_en -active_state 1
set_dft_signal -view spec         -type ScanDataIn  -port scan_in
set_dft_signal -view spec         -type ScanDataOut -port scan_out

set_scan_configuration -chain_count 1
create_test_protocol
dft_drc
preview_dft
insert_dft

write_scan_def -output dft/def/ecg_top_scan.def
write_verilog  dft/netlists/ecg_top_scan.v
```

### DFT Results

| Metric | Value |
|--------|-------|
| Scan flip-flops | [fill] |
| Scan chains | 1 |
| Scan coverage | [fill] % |
| DFT DRC violations | 0 (target) |

---

## Static Timing Analysis (STA)

### Tool: Synopsys PrimeTime

```bash
cd sta/
pt_shell -f scripts/primetime_sta.tcl | tee logs/sta.log
```

| Check | Corner | WNS | TNS | Status |
|-------|--------|-----|-----|--------|
| Setup | WCS (Slow) | [fill] ns | [fill] ns | [PASS/FAIL] |
| Hold | BCF (Fast) | [fill] ns | [fill] ns | [PASS/FAIL] |

---

## Formal Verification

### Tool: Synopsys Formality

```bash
cd formal/
fm_shell -f scripts/formality.tcl | tee logs/formality.log
```

| Check | Result |
|-------|--------|
| RTL vs Post-synthesis netlist | [EQUIVALENT] |
| Post-synthesis vs Post-DFT netlist | [EQUIVALENT] |

---

## Physical Design (PnR)

### Tool: Cadence Innovus

```bash
cd pnr/
innovus -batch -src scripts/innovus_pnr.tcl | tee logs/innovus.log
```

Flow steps: Floorplanning → Power Planning → Placement → CTS → Post-CTS Optimization → Routing → Post-Route ECO → Filler/Decap Insertion.

### PnR Results

| Metric | Value |
|--------|-------|
| Die area | [fill] µm² |
| Core utilization | [fill] % |
| Clock skew | [fill] ps |
| Routing DRC violations | 0 (target) |

---

## Sign-off and GDS Generation

```bash
cd signoff/
calibre -drc -hier scripts/calibre_drc.tcl    # Physical DRC
calibre -lvs -hier scripts/calibre_lvs.tcl    # Layout vs Schematic
calibre -xrc       scripts/calibre_rcx.tcl    # RC Parasitic Extraction
```

### Sign-off Checklist

- [ ] DRC clean — 0 violations
- [ ] LVS clean — netlist matches layout
- [ ] Post-layout STA with extracted parasitics — timing closure confirmed
- [ ] IR drop analysis — passed
- [ ] GDS-II exported: `signoff/gds/ecg_top_final.gds`

---

## Tool Flow Summary

```
ECG RTL (Verilog / SystemVerilog)
         │
         ▼
Functional Simulation ─────────── (VCS / Xcelium / ModelSim)
         │
         ▼
Logic Synthesis ────────────────── (Synopsys Design Compiler)
         │
         ▼
DFT Scan Insertion ─────────────── (Synopsys DFT Compiler)
         │
         ▼
Formal Equivalence Check ──────── (Synopsys Formality)
         │
         ▼
Pre-Layout STA ─────────────────── (Synopsys PrimeTime)
         │
         ▼
Place and Route ────────────────── (Cadence Innovus)
         │
         ▼
Post-Route STA + Sign-off ──────── (PrimeTime + Calibre)
         │
         ▼
GDS-II Output ✓
```

---

## Results and Reports

### Area Breakdown (Post-PnR)

| Module | Area (µm²) |
|--------|------------|
| cordic_mag_pipelined | [fill] |
| bandpass_filter_5_15 | [fill] |
| qrs_peak_detect | [fill] |
| rr_bpm_calc | [fill] |
| hrv_metrics (incl. 3x isqrt32) | [fill] |
| **Total** | [fill] |

### Power and Coverage

| Metric | Value |
|--------|-------|
| Dynamic power | [fill] mW |
| Leakage power | [fill] mW |
| Scan coverage | [fill] % |

---

## How to Run

```bash
git clone https://github.com/<your-username>/ECG-Signal-Processor-using-CORDIC-Algorithm-RTL-to-GDS-flow.git
cd ECG-Signal-Processor-using-CORDIC-Algorithm-RTL-to-GDS-flow

make sim       # RTL simulation only
make synth     # Logic synthesis
make dft       # DFT scan insertion
make sta       # Static timing analysis
make pnr       # Place and route
make signoff   # DRC / LVS / GDS export
make all       # Full flow end to end
make clean     # Remove all generated outputs
```

---

## Dependencies and Setup

| Tool / Package | Purpose |
|----------------|---------|
| Synopsys VCS or Cadence Xcelium | RTL simulation |
| Synopsys Design Compiler | Logic synthesis |
| Synopsys DFT Compiler | Scan insertion |
| Synopsys PrimeTime | STA |
| Synopsys Formality | Formal verification |
| Cadence Innovus | Place and Route |
| Mentor Calibre | DRC / LVS |
| Python 3 + wfdb + numpy | ECG data preparation |

```bash
# Install Python dependencies
pip install wfdb numpy

# Set PDK environment variable
export PDK_ROOT=/path/to/your/pdk
```

---

## Future Work

- [ ] Add CORDIC gain compensation (fixed-point: multiply by 39243, right-shift by 16)
- [ ] Extend to full Pan-Tompkins algorithm (squaring + moving-window integration stages)
- [ ] Implement adaptive thresholding for THRESHOLD_HIGH / THRESHOLD_LOW
- [ ] Add AXI4-Lite slave interface for SoC integration
- [ ] Develop UVM testbench with functional coverage and constrained-random stimulus
- [ ] Extend DFT to include MBIST for any embedded SRAM
- [ ] Tape-out using Sky130 open-source PDK via OpenLane

---

## References

1. J. Volder, "The CORDIC Trigonometric Computing Technique," *IRE Transactions on Electronic Computers*, 1959.
2. J. Pan and W. J. Tompkins, "A Real-Time QRS Detection Algorithm," *IEEE Transactions on Biomedical Engineering*, vol. 32, no. 3, 1985.
3. Task Force of ESC and NASPE, "Heart Rate Variability: Standards of Measurement, Physiological Interpretation, and Clinical Use," *Circulation*, 1996.
4. MIT-BIH Arrhythmia Database — [PhysioNet](https://physionet.org/content/mitdb/1.0.0/)
5. Synopsys DFT Compiler User Guide.
6. Cadence Innovus Implementation System User Guide.

---

## License

This project is licensed under the [MIT License](LICENSE).

---

