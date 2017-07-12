clc;
clear all;
close all;

addpath('/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila');

ORIGINAL_videos{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run1.avi';
ORIGINAL_videos{2} ='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run2.avi';
ORIGINAL_videos{3} ='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run3.avi';

TIMING_file =  'memento_final_loop_timings.txt';
TIMING_data = parse_timing_file(TIMING_file,2);

DIARY_file = 'memento_final_loop_timing_correlation_timecourses_diary.txt';

START_OFFSET = -0.05;
END_OFFSET = -0.05;
%DURATION = 10;

CUT_WIDTH_PIXELS = 0;

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
         
        ORIGINAL_vidObj.CurrentTime = max(eps,TIMING_data(k).FIRST_time(1)-START_OFFSET);       
        KRON_vidObj.CurrentTime = max(eps,TIMING_data(k).SECOND_time(1)-START_OFFSET);
        
        DURATION = TIMING_data(k).FIRST_time(2)-TIMING_data(k).FIRST_time(1) + END_OFFSET;
        
        total_corr=[];
        ORIGINAL_time=[];
        KRON_time = [];
        
        i=0;
        while ORIGINAL_vidObj.CurrentTime <= TIMING_data(k).FIRST_time(1)+DURATION
            i=i+1;
            try
                
                if mod(i,500)==0
                    fprintf('...step %i (elapsed time %s)\n',i,sec2min(i*(1/framerate)));
                end                
                
                t1=ORIGINAL_vidObj.CurrentTime;
                t2=KRON_vidObj.CurrentTime;     
                           
                data1 = single(readFrame(ORIGINAL_vidObj));
                data1 = data1(:,1:end-CUT_WIDTH_PIXELS,:);
                 
                data2 = single(readFrame(KRON_vidObj));
                
                a1 = image_zscore(squeeze(data1(:,:,1)));
                a2 = image_zscore(squeeze(data1(:,:,2)));
                a3 = image_zscore(squeeze(data1(:,:,3)));
                
                b1 = image_zscore(squeeze(data2(:,:,1)));
                b2 = image_zscore(squeeze(data2(:,:,2)));
                b3 = image_zscore(squeeze(data2(:,:,3)));
                
                val = sum(sum(a1.*b1 + a2.*b2 + a3.*b3))/(3*(N-1));                
                
                total_corr(i)=val;
                
                ORIGINAL_time(i)=t1;
                KRON_time(i)=t2;
                                
            catch err
                fprintf('cannot read frames: %s\n\n',err.message);
                break;
            end
        end                
        
        full_corrdata(k).total_corr = total_corr;
        full_corrdata(k).FIRST_time=ORIGINAL_time;
        full_corrdata(k).SECOND_time=KRON_time;
                      
        fprintf('>>>>>>>>> Segment ID %i is done\n',TIMING_data(k).ID);
        
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
save('repeating_loop_segment_data_stage2.mat','-v7.3');

close all;
for i=1:length(full_corrdata)
    
    try
        
        ID = full_corrdata(i).ID;
        IDs(i) = ID;
        plot_timecourses(ID,full_corrdata(i).total_corr,full_corrdata(i).FIRST_time,full_corrdata(i).SECOND_time);
        set(gca,'YLim',[0.5,1]);
        a=full_corrdata(i).total_corr;
        a(a==0)=[];
        a(isnan(a))=[];
        m(i)=mean(a);
        saveas(gcf,sprintf('loop_%i_correlation_ts.png',ID));
        
    catch err
        
    end
    
    
end
bar(IDs,m)
axis tight;
xlabel('Sequence ID');
ylabel('Mean correlation');



