clc;
clear all;
close all;

%prepare_roi_data_looper_stage1_deconv;
compute_roi_correlations_stage2_deconv([-50:1:20]);
compute_roi_correlations_stats_stage3_deconv(10000,7);