clc;
clear all;
close all;

addpath('/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila');

ORIGINAL_videos{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run1.avi';
ORIGINAL_videos{2} ='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run2.avi';
ORIGINAL_videos{3} ='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run3.avi';

TIMING_file = 'memento_final_loop_timings.txt';
TIMING_data = parse_timing_file(TIMING_file,2);
DIARY_file = 'memento_loop_timings_stage1_diary.txt';

for i=1:length(TIMING_data)
    TIMING_data(i).time_difference = abs(diff(TIMING_data(i).FIRST_time)-diff(TIMING_data(i).SECOND_time));
    if TIMING_data(i).time_difference>1e-8
        TIMING_data(i).ID
        error('time difference!!!')
    end
end

MAX_OFFSET = 0.40;

CHECK_START = -0.1;
CHECK_END = 4;

CUT_WIDTH_PIXELS = 0;

delete(DIARY_file);
diary(DIARY_file)

L=length(TIMING_data);
full_corrdata(L) = struct();

for k = 1:L
    
    try
        ORIGINAL_vidObj = VideoReader(ORIGINAL_videos{TIMING_data(k).FIRST_ses});
        KRON_vidObj = VideoReader(ORIGINAL_videos{TIMING_data(k).SECOND_ses});
        
        framerate = [ORIGINAL_vidObj.FrameRate,KRON_vidObj.FrameRate];
        if abs(diff(framerate))>1e-6
            error('incorrect framerate!')
        end
        framerate=framerate(1);
        
        full_corrdata(k).FIRST_file = ORIGINAL_videos{TIMING_data(k).FIRST_ses};
        full_corrdata(k).SECOND_file = ORIGINAL_videos{TIMING_data(k).SECOND_ses};
        
        full_corrdata(k).FIRST_ses = TIMING_data(k).FIRST_ses;
        full_corrdata(k).SECOND_ses = TIMING_data(k).SECOND_ses;
        
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
        fprintf(' FIRST: %s (file %s)\n',TIMING_data(k).FIRST_time_str,ORIGINAL_videos{TIMING_data(k).FIRST_ses});
        fprintf(' SECOND (keyframe): %s (file %s)\n',TIMING_data(k).SECOND_time_str,ORIGINAL_videos{TIMING_data(k).SECOND_ses});
        
        ORIGINAL_vidObj.CurrentTime = max(eps,TIMING_data(k).FIRST_time(1)+CHECK_START);
        ORIGINAL_data = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),'colormap',[],'timepoint',-1);
        
        i=0;
        while ORIGINAL_vidObj.CurrentTime < TIMING_data(k).FIRST_time(1)+CHECK_END
            i=i+1;
            try
                a = single(readFrame(ORIGINAL_vidObj));
                a = a(:,1:end-CUT_WIDTH_PIXELS,:);
                ORIGINAL_data(i).cdata = a;
                ORIGINAL_data(i).cdata(:,:,1)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,1)));
                ORIGINAL_data(i).cdata(:,:,2)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,2)));
                ORIGINAL_data(i).cdata(:,:,3)=image_zscore(squeeze(ORIGINAL_data(i).cdata(:,:,3)));
                ORIGINAL_data(i).cdata(isnan(ORIGINAL_data(i).cdata))=0;
                ORIGINAL_data(i).timepoint = ORIGINAL_vidObj.CurrentTime;
            catch err
                fprintf('cannot read frames: %s\n',err.message);
                break;
            end
        end
        
        KRON_vidObj.CurrentTime = max(eps,TIMING_data(k).SECOND_time(1)+CHECK_START);
        KRON_data = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),'colormap',[],'timepoint',-1);
        
        i=0;
        while KRON_vidObj.CurrentTime < TIMING_data(k).SECOND_time(1)+CHECK_END
            i=i+1;
            try
                KRON_data(i).cdata = single(readFrame(KRON_vidObj));
                KRON_data(i).cdata(:,:,1)=image_zscore(squeeze(KRON_data(i).cdata(:,:,1)));
                KRON_data(i).cdata(:,:,2)=image_zscore(squeeze(KRON_data(i).cdata(:,:,2)));
                KRON_data(i).cdata(:,:,3)=image_zscore(squeeze(KRON_data(i).cdata(:,:,3)));                                      
                KRON_data(i).cdata(isnan(KRON_data(i).cdata))=0;
                KRON_data(i).timepoint = KRON_vidObj.CurrentTime;
            catch err
                fprintf('cannot read frames: %s\n',err.message);
                break;
            end
        end
        
        MAX_SHIFT = round(MAX_OFFSET/(1/framerate));
        
        print_int = round(linspace(MAX_SHIFT/10,MAX_SHIFT,10));
        kk=1;
        
        fprintf(' ...running jitter analysis\n');
        
        ind = 1:min(length(ORIGINAL_data),length(KRON_data));
        
        if MAX_SHIFT>length(ind)*0.25
            error('increase duration time!')
        end
                                
        shift_ind = -MAX_SHIFT:MAX_SHIFT;        
        total_corr=nan(1,length(shift_ind));
        
        z=0;
        for i=shift_ind
            
            z=z+1;
            
            ind1 = ind;
            ind2 = ind;
            
            if i<0
                ind1(1:(-i))=[];
                ind2 = ind2(1:length(ind1));
            elseif i>0
                ind1((end-i+1):end)=[];
                ind2(1:i)=[];
            else                
            end
            
            if length(ind1)~=length(ind2)
                error('!!!')
            end
            
            data1 = ORIGINAL_data(ind1);
            data2 = KRON_data(ind2);
            
            if kk<11 && print_int(kk)==i
                fprintf(' ... %i/%i\n',kk,10);
                kk=kk+1;
            end
            
            cor = 0;
            for j=1:length(data1)
                a1 = squeeze(data1(j).cdata(:,:,1));
                a2 = squeeze(data1(j).cdata(:,:,2));
                a3 = squeeze(data1(j).cdata(:,:,3));
                
                b1 = squeeze(data2(j).cdata(:,:,1));
                b2 = squeeze(data2(j).cdata(:,:,2));
                b3 = squeeze(data2(j).cdata(:,:,3));
                
                val = sum(sum(a1.*b1 + a2.*b2 + a3.*b3))/(3*(N-1));
                
                if isnan(val)
                    error('NaN found, needs revision')
                end
                
                cor = cor + val;
            end
            cor = cor/length(data1);
            
            total_corr(z)=cor;
        end                
        fprintf(' ...done\n');
        
        clear data1 data2 KRON_data ORIGINAL_data;
        
        [m,ind]=max(total_corr);
        
        full_corrdata(k).shift_step = shift_ind;
        full_corrdata(k).total_correlation=total_corr;
        full_corrdata(k).best_mean_corr = m;
        full_corrdata(k).best_shift = shift_ind(ind);
        full_corrdata(k).best_shift_time = shift_ind(ind)*(1/framerate);
               
        fprintf('>>>>>>>>> Segment ID %i recommended shift is %.2fs (%i frames)\n',TIMING_data(k).ID,full_corrdata(k).best_shift_time,abs(full_corrdata(k).best_shift));
        
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
    
    close all;
    plot(full_corrdata(k).shift_step,full_corrdata(k).total_correlation,'bo-')
    title(full_corrdata(k).ID);
    axis tight;
    grid on;
    pause(1);
    
end

diary(DIARY_file)

diary off

for i=1:length(full_corrdata),
    fprintf('index %i, ID %i: %.3fs (%i frames), best correlation %f\n',i,full_corrdata(i).ID,full_corrdata(i).best_shift_time,full_corrdata(i).best_shift,full_corrdata(i).best_mean_corr);
end

% fprintf('\n');
% for i=1:length(full_corrdata),
%     fprintf('Segment ID = %i:\n  ORIG duration %0.2fs, CHRONO duration %.2fs (difference %.2fs)\n  ORIG ses = %i, timeslot (%s -> %s)\n  CHRONO ses = %i, timeslot (%s -> %s)\n',...
%         full_corrdata(i).ID,...
%         full_corrdata(i).durations(1),...
%         full_corrdata(i).durations(2),...
%         full_corrdata(i).duration_diff,...
%         TIMING_data(i).ORIGINAL_ses,...
%         sec2min(full_corrdata(i).beginning_best_time(1)),...
%         sec2min(full_corrdata(i).ending_best_time(1)),...
%         TIMING_data(i).KRON_ses,...
%         sec2min(full_corrdata(i).beginning_best_time(2)),...
%         sec2min(full_corrdata(i).ending_best_time(2))...
%         );
% end
% fprintf('\n');

save('memento_loop_final_timings_jitter.mat','-v7.3');
