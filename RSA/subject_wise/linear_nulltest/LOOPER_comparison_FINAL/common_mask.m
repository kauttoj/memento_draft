clear all;close all;clc;

data_path1{1}='/m/nbe/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/post_processed_normal/bramila/';
data_path1{2}='/m/nbe/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/post_processed_normal/bramila/';
data_path1{3}='/m/nbe/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/post_processed_normal/bramila/';

data_path2{1}='/m/nbe/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session1/post_processed_normal/bramila/';
data_path2{2}='/m/nbe/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session2/post_processed_normal/bramila/';
data_path2{3}='/m/nbe/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session3/post_processed_normal/bramila/';

S1 = {'S7','S8', 'S9', 'S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'};
S2 = {'KRON_2','KRON_3','KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13','KRON_15','KRON_16'};

EPI_mask1 = compute_analysis_mask(data_path1,S1);

EPI_mask2 = compute_analysis_mask(data_path2,S2);

EPI_mask = EPI_mask1.*EPI_mask2;

save_nii_oma(EPI_mask,'ORIG_KRON_mask.nii');