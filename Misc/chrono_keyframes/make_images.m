clear VIDEO_FILES
close all;
addpath('/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');
VIDEO_FILES.FIRST={...
    '/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M1_chrono_checkAudio_Mjpeg_April2015.avi',...
    '/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M2_chronological_stimuls_april2015.avi',...
    '/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M3_chronological_stimuls_april2015.avi'};
frame_plotter('memento_chrono_keyframes_FINAL.txt',VIDEO_FILES,'memento_chrono_keyframes',3);