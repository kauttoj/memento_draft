clc;
clear all;
close all;

addpath('/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila');

ORIGINAL_videos{1} = '/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run1.avi';
ORIGINAL_videos{2} ='/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run2.avi';
ORIGINAL_videos{3} ='/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run3.avi';

KRON_videos{1} = '/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M1_chrono_checkAudio_Mjpeg_April2015.avi';
KRON_videos{2} ='/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M2_chronological_stimuls_april2015.avi';
KRON_videos{3} ='/scratch/braindata/kauttoj2/Memento/2015/Stimulus/M3_chronological_stimuls_april2015.avi';

TIMING_file = 'jelena_segment_timings_corrected_stage1.txt';
TIMING_data = parse_timing_file(TIMING_file);
GOOD_segment = true(1,length(TIMING_data));
for i=1:length(TIMING_data)
    if TIMING_data(i).time_difference>2
        GOOD_segment(i)=false;
    end
end
GOOD_segment_ind = find(GOOD_segment);
fprintf('\n\n>>>>> %i bad segments were discarded\n\n',length(GOOD_segment)-nnz(GOOD_segment));

TIMING_data = TIMING_data(GOOD_segment);

DIARY_file = 'memento_video_segment_analysis_stage1.txt';

INITIAL_OFFSET = 0.8;

CUT_WIDTH_PIXELS = 1;

MAXVAL_TOL=0.999;

delete(DIARY_file);

diary(DIARY_file)

%funAlign_createpool( 2 );
%a(1)=TIMING_data(15);
%a(2)=TIMING_data(16);
%a(3)=TIMING_data(17);
%TIMING_data=a;

L=length(TIMING_data);
full_corrdata(L) = struct();


for k = 1:L
    
    try
        ORIGINAL_vidObj = VideoReader(ORIGINAL_videos{TIMING_data(k).ORIGINAL_ses});
        KRON_vidObj = VideoReader(KRON_videos{TIMING_data(k).KRON_ses});
        
        full_corrdata(k).ORIGINAL_file = ORIGINAL_videos{TIMING_data(k).ORIGINAL_ses};
        full_corrdata(k).KRON_file = KRON_videos{TIMING_data(k).KRON_ses};
        
        full_corrdata(k).ORIGINAL_ses = TIMING_data(k).ORIGINAL_ses;
        full_corrdata(k).KRON_ses = TIMING_data(k).KRON_ses;
        
        vidHeight = ORIGINAL_vidObj.Height;
        vidWidth = ORIGINAL_vidObj.Width-CUT_WIDTH_PIXELS;
        N = vidHeight*vidWidth;
        
        vidHeight1 = KRON_vidObj.Height;
        vidWidth1 = KRON_vidObj.Width;
        N1 = vidHeight1*vidWidth1;
        
        if vidHeight~=vidHeight1 || vidWidth~=vidWidth1
            warning('Video size mismatch!')
        end
        
        full_corrdata(k).ID = TIMING_data(k).ID;
        
        fprintf(' Analyzing segment ID %i\n',TIMING_data(k).ID);
        fprintf(' ORIGINAL: %s (file %s)\n',TIMING_data(k).ORIGINAL_time_str,ORIGINAL_videos{TIMING_data(k).ORIGINAL_ses});
        fprintf(' KRON: %s (file %s)\n',TIMING_data(k).KRON_time_str,KRON_videos{TIMING_data(k).KRON_ses});
        
        ORIGINAL_vidObj.CurrentTime = max(eps,TIMING_data(k).ORIGINAL_time(1)-INITIAL_OFFSET);
        ORIGINAL_data = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),'colormap',[],'timepoint',-1);
        ORIGINAL_data1=ORIGINAL_data;
        
        i=0;
        while ORIGINAL_vidObj.CurrentTime < TIMING_data(k).ORIGINAL_time(1)+INITIAL_OFFSET
            i=i+1;
            try
                a = single(readFrame(ORIGINAL_vidObj));
                a = a(:,1:end-CUT_WIDTH_PIXELS,:);
                ORIGINAL_data(i).cdata = a;
                ORIGINAL_data1(i).cdata = a;
                ORIGINAL_data(i).cdata(:,:,1)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,1)));
                ORIGINAL_data(i).cdata(:,:,2)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,2)));
                ORIGINAL_data(i).cdata(:,:,3)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,3)));
                ORIGINAL_data(i).timepoint = ORIGINAL_vidObj.CurrentTime;
                ORIGINAL_data1(i).timepoint=ORIGINAL_data(i).timepoint;
            catch err
                fprintf('cannot read frames: %s\n',err.message);
                break;
            end
        end
        
        KRON_data = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),'colormap',[],'timepoint',-1);
        KRON_data1=KRON_data;
        KRON_vidObj.CurrentTime = max(eps,TIMING_data(k).KRON_time(1)-INITIAL_OFFSET);
        
        i=0;
        while KRON_vidObj.CurrentTime < TIMING_data(k).KRON_time(1)+INITIAL_OFFSET
            i=i+1;
            try
                KRON_data(i).cdata = single(readFrame(KRON_vidObj));
                KRON_data1(i).cdata=KRON_data(i).cdata;
                KRON_data(i).cdata(:,:,1)=image_zscore(squeeze(KRON_data(i).cdata(:,:,1)));
                KRON_data(i).cdata(:,:,2)=image_zscore(squeeze(KRON_data(i).cdata(:,:,2)));
                KRON_data(i).cdata(:,:,3)=image_zscore(squeeze(KRON_data(i).cdata(:,:,3)));
                KRON_data(i).timepoint = KRON_vidObj.CurrentTime;
                KRON_data1(i).timepoint=KRON_data(i).timepoint;
            catch err
                fprintf('cannot read frames: %s\n',err.message);
                break;
            end
        end
        
        print_int = round(linspace(length(ORIGINAL_data)/10,length(ORIGINAL_data),10));
        kk=1;
        fprintf(' ...running correlation analysis (beginning)\n');
        corrdata = nan(length(ORIGINAL_data),length(KRON_data));
        for i=1:length(ORIGINAL_data)
            if kk<11 && print_int(kk)==i
                fprintf(' ... %i/%i\n',kk,10);
                kk=kk+1;
            end
            a1 = squeeze(ORIGINAL_data(i).cdata(:,:,1));
            a2 = squeeze(ORIGINAL_data(i).cdata(:,:,2));
            a3 = squeeze(ORIGINAL_data(i).cdata(:,:,3));
            for j=1:length(KRON_data)
                b1 = squeeze(KRON_data(j).cdata(:,:,1));
                b2 = squeeze(KRON_data(j).cdata(:,:,2));
                b3 = squeeze(KRON_data(j).cdata(:,:,3));
                val = sum(sum(a1.*b1 + a2.*b2 + a3.*b3))/(3*(N-1));
                corrdata(i,j)=val;
            end
        end
        fprintf(' ...done\n');
        
        
        full_corrdata(k).beginning=corrdata;
        %full_corrdata(k).first = first;
        
        maxval = MAXVAL_TOL*max(corrdata(:));
        % maxval = 0.98*maxval;
        [i,j]=ind2sub(size(corrdata),find(maxval <= corrdata));
        z = i+j;
        [~,a]=min(z);
        i=i(a);
        j=j(a);
        
        full_corrdata(k).beginning_best_index = [i,j,corrdata(i,j)];
        full_corrdata(k).beginning_best_time = [ORIGINAL_data(i).timepoint,KRON_data(j).timepoint,corrdata(i,j)];
        
        print_frames(corrdata,ORIGINAL_data1,KRON_data1,i,j,['Memento_segment_',num2str(full_corrdata(k).ID),'_beginning'],vidWidth,vidHeight,TIMING_data(k).ORIGINAL_ses,TIMING_data(k).KRON_ses);
        
        %clear ORIGINAL_data KRON_data ORIGINAL_data1 KRON_data1 corrdata;
        
        %%-----
        
        ORIGINAL_vidObj.CurrentTime = max(eps,TIMING_data(k).ORIGINAL_time(2)-INITIAL_OFFSET);
        ORIGINAL_data = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),'colormap',[],'timepoint',-1);
        ORIGINAL_data1=ORIGINAL_data;
        
        i=0;
        while ORIGINAL_vidObj.CurrentTime < TIMING_data(k).ORIGINAL_time(2)+INITIAL_OFFSET
            i=i+1;
            try
                a = single(readFrame(ORIGINAL_vidObj));
                a = a(:,1:end-CUT_WIDTH_PIXELS,:);
                ORIGINAL_data(i).cdata = a;
                ORIGINAL_data1(i).cdata = a;
                ORIGINAL_data(i).cdata(:,:,1)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,1)));
                ORIGINAL_data(i).cdata(:,:,2)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,2)));
                ORIGINAL_data(i).cdata(:,:,3)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,3)));
                ORIGINAL_data(i).timepoint = ORIGINAL_vidObj.CurrentTime;
                ORIGINAL_data1(i).timepoint=ORIGINAL_data(i).timepoint;
            catch err
                fprintf('cannot read frames: %s\n',err.message);
                break;
            end
        end
        
        KRON_data = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),'colormap',[],'timepoint',-1);
        KRON_data1=KRON_data;
        KRON_vidObj.CurrentTime = max(0,TIMING_data(k).KRON_time(2)-INITIAL_OFFSET);
        
        i=0;
        while KRON_vidObj.CurrentTime < TIMING_data(k).KRON_time(2)+INITIAL_OFFSET
            i=i+1;
            try
                KRON_data(i).cdata = single(readFrame(KRON_vidObj));
                KRON_data1(i).cdata=KRON_data(i).cdata;
                KRON_data(i).cdata(:,:,1)=image_zscore(squeeze(KRON_data(i).cdata(:,:,1)));
                KRON_data(i).cdata(:,:,2)=image_zscore(squeeze(KRON_data(i).cdata(:,:,2)));
                KRON_data(i).cdata(:,:,3)=image_zscore(squeeze(KRON_data(i).cdata(:,:,3)));
                KRON_data(i).timepoint = KRON_vidObj.CurrentTime;
                KRON_data1(i).timepoint=KRON_data(i).timepoint;
            catch err
                fprintf('cannot read frames: %s\n',err.message);
                break;
            end
        end
        
        print_int = round(linspace(length(ORIGINAL_data)/10,length(ORIGINAL_data),10));
        kk=1;
        fprintf('\n ...running correlation analysis (ending)\n');
        corrdata = nan(length(ORIGINAL_data),length(KRON_data));
        for i=1:length(ORIGINAL_data)
            if kk<11 && print_int(kk)==i
                fprintf(' ... %i/%i\n',kk,10);
                kk=kk+1;
            end
            a1 = squeeze(ORIGINAL_data(i).cdata(:,:,1));
            a2 = squeeze(ORIGINAL_data(i).cdata(:,:,2));
            a3 = squeeze(ORIGINAL_data(i).cdata(:,:,3));
            for j=1:length(KRON_data)
                b1 = squeeze(KRON_data(j).cdata(:,:,1));
                b2 = squeeze(KRON_data(j).cdata(:,:,2));
                b3 = squeeze(KRON_data(j).cdata(:,:,3));
                val = sum(sum(a1.*b1 + a2.*b2 + a3.*b3))/(3*(N-1));
                corrdata(i,j)=val;
            end
        end
        fprintf(' ...done\n');
        
        
        full_corrdata(k).ending=corrdata;
        %full_corrdata(k).first = first;
        
        maxval = MAXVAL_TOL*max(corrdata(:));
        %maxval = 0.98*maxval;
        [i,j]=ind2sub(size(corrdata),find(maxval <= corrdata));
        z = i+j;
        [~,a]=max(z);
        i=i(a);
        j=j(a);
        
        full_corrdata(k).ending_best_index = [i,j,corrdata(i,j)];
        full_corrdata(k).ending_best_time = [ORIGINAL_data(i).timepoint,KRON_data(j).timepoint,corrdata(i,j)];
        
        print_frames(corrdata,ORIGINAL_data1,KRON_data1,i,j,['Memento_segment_',num2str(full_corrdata(k).ID),'_ending'],vidWidth,vidHeight,TIMING_data(k).ORIGINAL_ses,TIMING_data(k).KRON_ses);
        
        %clear ORIGINAL_data KRON_data ORIGINAL_data1 KRON_data1 corrdata;
        
        full_corrdata(k).durations = [full_corrdata(k).ending_best_time(1)-full_corrdata(k).beginning_best_time(1),...
            full_corrdata(k).ending_best_time(2)-full_corrdata(k).beginning_best_time(2)];
        full_corrdata(k).duration_diff = abs(diff(full_corrdata(k).durations));
        
        fprintf('>>>>>>>>> Segment ID %i final difference is %.2s\n',TIMING_data(k).ID,full_corrdata(k).duration_diff);

        %fprintf('Best correlation is %.4f is at %.2fs and %.2fs (inds %i and %i)\n\n',corrdata(i,j),times_first{iter}(i),times_second{iter}(j),i,j);
        
    catch err
        warning('Segment ID %i FAILED: %s',TIMING_data(k).ID,err.message);
    end
    %top = all_corrdata{iter};
    %top(top<0.95)=0;
    %handle = plot_matrix_simple(top,sprintf('Memento key-frame %i timing (>0.95)',iter),round(times_second{iter}),round(times_first{iter}));
    %set(handle,'Position',[75   275   921   819]);
    %saveas(handle,sprintf('Memento_keyframe_%i_matrix_over95.png',iter));
    %close(handle);
    
end

diary(DIARY_file)

diary off

for i=1:length(full_corrdata),fprintf('index %i, ID %i: %f\n',i,full_corrdata(i).ID,full_corrdata(i).duration_diff);end

fprintf('\n');
for i=1:length(full_corrdata),
    fprintf('Segment ID = %i:\n  ORIG duration %0.2fs, CHRONO duration %.2fs (difference %.2fs)\n  ORIG ses = %i, timeslot (%s -> %s)\n  CHRONO ses = %i, timeslot (%s -> %s)\n',...
        full_corrdata(i).ID,...        
        full_corrdata(i).durations(1),...
        full_corrdata(i).durations(2),...
        full_corrdata(i).duration_diff,...
        TIMING_data(i).ORIGINAL_ses,...
        sec2min(full_corrdata(i).beginning_best_time(1)),...
        sec2min(full_corrdata(i).ending_best_time(1)),...
        TIMING_data(i).KRON_ses,...
        sec2min(full_corrdata(i).beginning_best_time(2)),...
        sec2min(full_corrdata(i).ending_best_time(2))...
        );
end
fprintf('\n');

save('repeating_segment_data.mat','-v7.3');

%diary off
% image(s(5).cdata)
% set(gcf,'position',[150 150 vidObj.Width vidObj.Height]);
% set(gca,'units','pixels');
% set(gca,'position',[0 0 vidObj.Width vidObj.Height]);
% movie(s,1,vidObj.FrameRate);