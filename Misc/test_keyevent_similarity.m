clc;
clear variables;
close all;

DELAYS = [20,6];
STEPPING = 0.5;

if ~exist('TIMING_keys_tempdata.mat')
    TIMING_file = 'loop_keyframes/memento_loop_keyframe_timings_FINAL.txt';
    TIMING_keys_orig = parse_timing_file(TIMING_file,2);
    
    TIMING_file = 'chrono_keyframes/memento_chrono_keyframes_FINAL.txt';
    TIMING_keys_kron = parse_timing_file(TIMING_file,3);
else
    load('TIMING_keys_tempdata.mat')
end

VIDEOS_orig{1} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run1.avi';
VIDEOS_orig{2} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run2.avi';
VIDEOS_orig{3} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run3.avi';

VIDEOS_kron{1} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M1_chrono_checkAudio_Mjpeg_April2015.avi';
VIDEOS_kron{2} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M2_chronological_stimuls_april2015.avi';
VIDEOS_kron{3} = '/m/nbe/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M3_chronological_stimuls_april2015.avi';

FILEPRE = 'Memento_ORIG_keyframe_ID_';
TITLEPRE = 'ORIG ID ';

currtime = -DELAYS(1);
endtime = currtime+sum(DELAYS);

kk=0;
while currtime<=endtime
    
    for i = 1:length(TIMING_keys_orig)
        
        SES = TIMING_keys_orig(i).SECOND_ses;
        ID = TIMING_keys_orig(i).ID;
        targettime = TIMING_keys_orig(i).SECOND_time(1);
        videoname = [FILEPRE,num2str(ID),'.avi'];
        videoFReader = VideoReader(VIDEOS_orig{SES});
        videoFReader.CurrentTime = targettime - currtime;
        
        a = videoFReader.CurrentTime;
        endtime = a+STEPPING;
        temp = 0;
        k=0;
        while a<=endtime
            frame = readFrame(videoFReader);
            temp = temp + frame;
            k=k+1;
            a = videoFReader.CurrentTime;
        end
        temp = temp/k;
        allframes{i,1}=squeeze(temp(:,:,1));
        allframes{i,2}=squeeze(temp(:,:,2));
        allframes{i,3}=squeeze(temp(:,:,3));
        
    end
    
    kk=kk+1;
    corvals_orig(kk)=0;
    for k=1:3
        corvals_orig(kk)=corvals_orig(kk) + mean(corr(allframes{i,k}))/3;
    end
    currtime = currtime +STEPPING;
    
end

FILEPRE = 'Memento_KRON_keyevent_ID_';
TITLEPRE = 'KRON ID ';

currtime = -DELAYS(1);
endtime = currtime+sum(DELAYS);

kk=0;
while currtime<=endtime
    
    for i = 1:length(TIMING_keys_orig)
        
        SES = TIMING_keys_kron(i).FIRST_ses;
        ID = TIMING_keys_kron(i).ID;
        targettime = TIMING_keys_kron(i).FIRST_time(1);
        videoname = [FILEPRE,num2str(ID),'.avi'];
        videoFReader = VideoReader(VIDEOS_kron{SES});
        videoFReader.CurrentTime = targettime - currtime;
        
        a = videoFReader.CurrentTime;
        endtime = a+STEPPING;
        temp = 0;
        k=0;
        while a<=endtime
            frame = readFrame(videoFReader);
            temp = temp + frame;
            k=k+1;
            a = videoFReader.CurrentTime;
        end
        temp = temp/k;
        allframes{i,1}=squeeze(temp(:,:,1));
        allframes{i,2}=squeeze(temp(:,:,2));
        allframes{i,3}=squeeze(temp(:,:,3));
        
    end
    
    kk=kk+1;
    corvals_kron(kk)=0;
    for k=1:3
        corvals_kron(kk)=corvals_kron(kk) + mean(corr(allframes{i,k}))/3;
    end
    currtime = currtime +STEPPING;
    
end
