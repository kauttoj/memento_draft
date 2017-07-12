function play_keyframes(TIMES,DELAY,DURATION)
%PLAY_VIDEO Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    DURATION=-1;
end

VIDEOS{1} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run1.avi';
VIDEOS{2} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run2.avi';
VIDEOS{3} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run3.avi';

VIDEOS{1} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M1_chrono_checkAudio_Mjpeg_April2015.avi';
VIDEOS{2} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M2_chronological_stimuls_april2015.avi';
VIDEOS{3} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M3_chronological_stimuls_april2015.avi';


for i=1:length(TIMES)
    
    tim = TIMES(i).SECOND_time(1)+DELAY;
    if DURATION<0
        dur = diff(TIMES(i).SECOND_time)-DELAY;
    else
        dur = DURATION;
    end
    vid = TIMES(i).SECOND_ses;
    
    close all;
    play_video_simple(VIDEOS(vid),tim,dur,1);
    
end
