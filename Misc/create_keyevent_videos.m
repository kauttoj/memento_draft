clc;
clear variables;
close all;

DELAYS = [-20,20];

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

for i = 1:length(TIMING_keys_orig)
    
    close all;
    SES = TIMING_keys_orig(i).SECOND_ses;
    ID = TIMING_keys_orig(i).ID;
    targettime = TIMING_keys_orig(i).SECOND_time(1);
    videoname = [FILEPRE,num2str(ID),'.avi'];
    videoFReader = VideoReader(VIDEOS_orig{SES});
    videoFReader.CurrentTime = targettime - DELAYS(1);
    vidWidth = videoFReader.Width;
    vidHeight = videoFReader.Height;
    frameRate=videoFReader.FrameRate;
    
    videoout = VideoWriter(videoname,'Motion JPEG AVI');
    videoout.FrameRate = frameRate;%video.FrameRate;
    
    scrsz = get(groot,'ScreenSize');
    
    pos = [scrsz(3)/2,scrsz(4)/2,scrsz(3)/2,scrsz(4)/2];
    handle = figure('Position',pos);
    handle_ax=gca;
    set(handle_ax,'XTick',[],'YTick',[],'box','on');
    tithand = title([TITLEPRE,num2str(ID)]);
    
    fprintf('\nORIG ID %i starting\n',ID);
    currtime = videoFReader.CurrentTime;
    endtime = targettime + DELAYS(2);
    k=0;
    open(videoout);
    while currtime<=endtime
        try
            frame = readFrame(videoFReader);
            if exist(hand,'var')
                delete(hand);
            end
            hand = image(frame,'Parent',handle_ax);
            set(tithand,'String',['Onset: ',(currtime - targettime),'s']);
            currtime = currtime + (1/frameRate(1));
            writeVideo(videoout,frame);
            k = k+1;
            if mod(k,frameRate)==0                
                fprintf('... time %f\n',currtime - targettime);
            end            
        catch err
            warning('Stopped at frame %i: %s',k,err.message);
            break;
        end
        currtime = currtime + (1/frameRate);
    end
    
    close(videoout);
    
end

FILEPRE = 'Memento_KRON_keyevent_ID_';
TITLEPRE = 'KRON ID ';

for i = 1:length(TIMING_keys_kron)
    
    close all;
    SES = TIMING_keys_kron(i).SECOND_ses;
    ID = TIMING_keys_kron(i).ID;
    targettime = TIMING_keys_kron(i).FIRST_time(1);
    videoname = [FILEPRE,num2str(ID),'.avi'];
    videoFReader = VideoReader(VIDEOS_kron{SES});
    videoFReader.CurrentTime = targettime - DELAYS(1);
    vidWidth = videoFReader.Width;
    vidHeight = videoFReader.Height;
    frameRate=videoFReader.FrameRate;
    
    videoout = VideoWriter(videoname,'Motion JPEG AVI');
    videoout.FrameRate = frameRate;%video.FrameRate;    
    
    scrsz = get(groot,'ScreenSize');
    
    pos = [scrsz(3)/2,scrsz(4)/2,scrsz(3)/2,scrsz(4)/2];
    handle = figure('Position',pos);
    handle_ax=gca;
    set(handle_ax,'XTick',[],'YTick',[],'box','on');
    tithand = title([TITLEPRE,num2str(ID)]);
    
    fprintf('\nKRON ID %i starting\n',ID);
    currtime = videoFReader.CurrentTime;
    endtime = targettime + DELAYS(2);
    k=0;
    open(videoout);
    while currtime<=endtime
        try
            frame = readFrame(videoFReader);
            if exist(hand,'var')
                delete(hand);
            end
            hand = image(frame,'Parent',handle_ax);
            set(tithand,'String',['Onset: ',(currtime - targettime),'s']);
            currtime = currtime + (1/frameRate(1));
            writeVideo(videoout,frame);
            k = k+1;
            if mod(k,frameRate)==0                
                fprintf('... time %f\n',currtime - targettime);
            end
        catch err
            warning('Stopped at frame %i: %s',k,err.message);
            break;
        end
        currtime = currtime + (1/frameRate);
    end
    
    close(videoout);
    
end
