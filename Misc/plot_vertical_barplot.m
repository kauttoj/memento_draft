clc;
clear variables;
close all;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');

% initial & keyframe onset
cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/loop_keyframes');
TIMING_file = 'memento_loop_keyframe_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,2);
cd ..

% bw onset
cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/bw_scenes');
TIMING_file = 'memento_bw_scene_timings_FINAL.txt';
TIMING_data_bw = parse_timing_file(TIMING_file,3);
cd ..
