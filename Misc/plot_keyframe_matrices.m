clc;
clear all;
close all;

addpath('/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila');

ORIGINAL_videos{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run1.avi';
ORIGINAL_videos{2} ='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run2.avi';
ORIGINAL_videos{3} ='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/Stimulus/Memento3part_run3.avi';

TIMING_file =  'loop_keyframes/memento_loop_keyframe_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,2);

START_OFFSET = 8;
END_OFFSET = 8;
%DURATION = 10;

CUT_WIDTH_PIXELS = 0;

L=length(TIMING_data);
full_corrdata(L) = struct();

for k = 1:L
    
    try
        ORIGINAL_vidObj1 = VideoReader(ORIGINAL_videos{TIMING_data(k).FIRST_ses});
        ORIGINAL_vidObj2 = VideoReader(ORIGINAL_videos{TIMING_data(k).SECOND_ses});
        
        framerate = [ORIGINAL_vidObj1.FrameRate,ORIGINAL_vidObj2.FrameRate];
        if abs(diff(framerate))>1e-6
            error('incorrect framerate!')
        end
        framerate=framerate(1);
        
        full_corrdata(k).FIRST_file = ORIGINAL_videos{TIMING_data(k).FIRST_ses};
        full_corrdata(k).SECOND_file = ORIGINAL_videos{TIMING_data(k).SECOND_ses};
        
        full_corrdata(k).FIRST_ses = TIMING_data(k).FIRST_ses;
        full_corrdata(k).SECOND_ses = TIMING_data(k).SECOND_ses;
        
        vidHeight = ORIGINAL_vidObj1.Height;
        vidWidth = ORIGINAL_vidObj1.Width-CUT_WIDTH_PIXELS;
        N = vidHeight*vidWidth;
        
        vidHeight1 = ORIGINAL_vidObj2.Height;
        vidWidth1 = ORIGINAL_vidObj2.Width;
        N1 = vidHeight1*vidWidth1;
        
        if vidHeight~=vidHeight1 || vidWidth~=vidWidth1
            warning('Video size mismatch!')
        end
        
        full_corrdata(k).ID = TIMING_data(k).ID;
        
        fprintf(' Analyzing segment ID %i\n',TIMING_data(k).ID);
        fprintf(' FIRST: %s (file %s)\n',TIMING_data(k).FIRST_time_str,ORIGINAL_videos{TIMING_data(k).FIRST_ses});
        fprintf(' SECOND (keyframe): %s (file %s)\n',TIMING_data(k).SECOND_time_str,ORIGINAL_videos{TIMING_data(k).SECOND_ses});
         
        ORIGINAL_vidObj2.CurrentTime = max(eps,TIMING_data(k).SECOND_time(1)-START_OFFSET);
        
        col = 0;
        while ORIGINAL_vidObj2.CurrentTime <= TIMING_data(k).SECOND_time(2)+END_OFFSET
            col=col+1;
            
            coltime(col)=ORIGINAL_vidObj2.CurrentTime;
            
            a = single(readFrame(ORIGINAL_vidObj2));
            
            b1 = image_zscore(squeeze(a(:,:,1)));
            b2 = image_zscore(squeeze(a(:,:,2)));
            b3 = image_zscore(squeeze(a(:,:,3)));
            
            data2(:,:,1,col)=b1;
            data2(:,:,2,col)=b2;
            data2(:,:,3,col)=b3;
            
        end        
        
        clear b1 b2 b3;
        data2(isnan(data2))=0;
        
        ORIGINAL_vidObj1.CurrentTime = max(eps,TIMING_data(k).FIRST_time(1)-START_OFFSET);       
                               
        row=0;
        while ORIGINAL_vidObj1.CurrentTime <= TIMING_data(k).FIRST_time(2)+END_OFFSET
            row=row+1;
            
            if mod(row,10)==0
                fprintf('...row %i\n',row);
            end
            
            rowtime(row) = ORIGINAL_vidObj1.CurrentTime;
            
            data1 = single(readFrame(ORIGINAL_vidObj1));
            data1 = data1(:,1:end-CUT_WIDTH_PIXELS,:);       
            
            a1 = image_zscore(squeeze(data1(:,:,1)));
            a2 = image_zscore(squeeze(data1(:,:,2)));
            a3 = image_zscore(squeeze(data1(:,:,3)));
            
            a1(isnan(a1))=0;
            a2(isnan(a2))=0;
            a3(isnan(a3))=0;
                                    
            for col=1:size(data2,4)
                
                b1 = squeeze(data2(:,:,1,col));
                b2 = squeeze(data2(:,:,2,col));
                b3 = squeeze(data2(:,:,3,col));
                
                val = sum(sum(a1.*b1 + a2.*b2 + a3.*b3))/(3*(N-1));
                corrmat(row,col)=val;
                
            end
                                
        end                
        
        clear data2 data1 a1 a2 a3 ORIGINAL_vidObj1 ORIGINAL_vidObj2;
        
        full_corrdata(k).corrmat = single(corrmat);
        full_corrdata(k).FIRST_times=rowtime;
        full_corrdata(k).SECOND_times=coltime;
        
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

save('loop_keyframe_correlation_matrices.mat','-v7.3');

% close all;
% for i=1:length(full_corrdata)
%     
%     try
%         
%         ID = full_corrdata(i).ID;
%         IDs(i) = ID;
%         plot_timecourses(ID,full_corrdata(i).total_corr,full_corrdata(i).FIRST_time,full_corrdata(i).SECOND_time);
%         set(gca,'YLim',[0.5,1]);
%         a=full_corrdata(i).total_corr;
%         a(a==0)=[];
%         a(isnan(a))=[];
%         m(i)=mean(a);
%         saveas(gcf,sprintf('loop_%i_correlation_ts.png',ID));
%         
%     catch err
%         
%     end
%     
%     
% end
% bar(IDs,m)
% axis tight;
% xlabel('Sequence ID');
% ylabel('Mean correlation');


