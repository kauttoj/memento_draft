clc;
clear all;
close all;

prepare_roi_data_looper_stage1;
DELAYS = -70:30;
compute_roi_correlations_stage2(DELAYS,3);

