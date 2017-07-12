clc;
clear all;
close all;

pre_keyframes_onset = [...
    144.48,...
    170.8,...
    402.03,...
    637.077,...
    1007.191,...
    1321.705 - 1,...
    1879.798 + 10,...
    2649.467,...
    2841.695,...
    3004.893,...
    3092.618,...
    3307.149,...
    3736.151,...
    4025.67 + 2,...
    4225.648,...
    4745.739,...
    5186.296];

pre_keyframes_end = [...
    146.08,...
    177.8,...
    409.83,...
    654.877,...
    1024.891,...
    1332.605,...
    1899.098,...
    2661.567,...
    2849.295,...
    3026.928,...
    3104.818,...
    3381.549,...
    3750.051,...
    4043.77,...
    4235.348,...
    4754.439,...
    5195.396];

pre_keyframes_end = pre_keyframes_onset + 15;
pre_keyframes_onset=pre_keyframes_onset - 10;


keyframes_onset = [...
    360.00,...
    578.00 + 2,...
    919.00 + 2,...
    1272.00+1,...
    1530.00+2,...
    1818.00,...
    2582.00 + 2,...
    2990.00,...
    3078.00,...
    3249.00,...
    3410.00 + 1,...
    3599.00,...
    4195.00,...leikkautuu
    4417.00+3,...
    4696.00,...
    4970.00,...
    6341.00];

keyframes_end = [...
    366,...
    585,...
    927,...
    1279,...
    1537,...
    1822,...
    2590,...
    2995,...
    3080,...
    3256,...
    3422,...
    3605,...
    4198,... % leikkautuu
    4445,...
    4700,...
    4975,...
    6348];

keyframes_end = keyframes_onset + 10;
keyframes_onset = keyframes_onset - 10;

% pre_keyframes_onset = 173 - 4;
% pre_keyframes_end = 173 + 7;
% 
% keyframes_onset = 582 - 5;
% keyframes_end  = 582 + 6;

if length(keyframes_end)~=length(pre_keyframes_end)
   error('size does not match');
end

L = length(keyframes_end);

if any(keyframes_onset - pre_keyframes_onset < 0)
    error('causality incorrect');
end

fprintf('Reading video....');
%vidObj = VideoReader('/media/MyBook/memento/PresentationFiles/Stimulus/MEMENTO_clip.avi'); %'/media/MyBook/memento/stimulus/MEMENTO.m4v');
vidObj = VideoReader('MEMENTO.m4v');
fprintf(' done\n');

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
N = vidHeight*vidWidth;

delete('memento_videoanalysis_add.txt');
diary('memento_videoanalysis_add.txt')

for iter = [7],
    
    fprintf('--- Processing loop %i/%i -----\n',iter,L);
    
    s_first = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),...
        'colormap',[],'timepoint',-1);
    s_second = struct('cdata',zeros(vidHeight,vidWidth,3,'single'),...
        'colormap',[],'timepoint',-1);
    
    k = 1;
    vidObj.CurrentTime = pre_keyframes_onset(iter);
    while vidObj.CurrentTime < pre_keyframes_end(iter)
        s_first(k).cdata = single(readFrame(vidObj));
        s_first(k).cdata(:,:,1)=image_zscore(squeeze(s_first(k).cdata(:,:,1)));
        s_first(k).cdata(:,:,2)=image_zscore(squeeze(s_first(k).cdata(:,:,2)));
        s_first(k).cdata(:,:,3)=image_zscore(squeeze(s_first(k).cdata(:,:,3)));        
        s_first(k).timepoint = vidObj.CurrentTime;
        k = k+1;
    end
    
    for i=1:length(s_first)
       times_first{iter}(i)= s_first(i).timepoint;
    end
    
    
    k = 1;
    vidObj.CurrentTime = keyframes_onset(iter);
    while vidObj.CurrentTime < keyframes_end(iter)
        s_second(k).cdata = single(readFrame(vidObj));
        s_second(k).cdata(:,:,1)=image_zscore(squeeze(s_second(k).cdata(:,:,1)));
        s_second(k).cdata(:,:,2)=image_zscore(squeeze(s_second(k).cdata(:,:,2)));
        s_second(k).cdata(:,:,3)=image_zscore(squeeze(s_second(k).cdata(:,:,3)));   
        s_second(k).timepoint = vidObj.CurrentTime;
        k = k+1;
    end
    
    for i=1:length(s_second)
       times_second{iter}(i)= s_second(i).timepoint;
    end    
    
    success = false;
    
    fprintf(' ...running correlation analysis\n');
    corrdata = nan(length(s_first),length(s_second));
    for i=1:length(s_first)
        if mod(i,50)==0
            fprintf(' ... %i/%i\n',i,length(s_first));
        end
        a1 = squeeze(s_first(i).cdata(:,:,1));
        a2 = squeeze(s_first(i).cdata(:,:,2));
        a3 = squeeze(s_first(i).cdata(:,:,3));
        for j=1:length(s_second)
            b1 = squeeze(s_second(j).cdata(:,:,1));
            b2 = squeeze(s_second(j).cdata(:,:,2));
            b3 = squeeze(s_second(j).cdata(:,:,3));
            val = sum(sum(a1.*b1 + a2.*b2 + a3.*b3))/(3*(N-1));
            if val>0.95 && ~success
                success = true;
                first = [i,j];
                fprintf(' ... first identical frame found: %.1fs and %.1fs (frames %i and %i)\n',s_first(i).timepoint,s_second(j).timepoint,i,j);
            end
            
            corrdata(i,j)=val;
        end
    end
    fprintf(' ...done\n');
    
    all_corrdata{iter}=corrdata;
    
    maxval = max(corrdata(:));
    
    [i,j]=ind2sub(size(corrdata),find(maxval == corrdata));
    if maxval~=corrdata(i,j)
        error('does not match')
    end
    if success
        fprintf('First >0.95 correlation is %.4f is at %.2fs and %.2fs (inds %i and %i)\n', corrdata(first(1),first(2)),times_first{iter}(first(1)),times_second{iter}(first(2)),first(1),first(2));
        
        figure;
        set(gcf,'position',[10 10 2*vidObj.Width 1.05*vidObj.Height]);
        
        subplot(1,2,1);
        image(s_first(first(1)).cdata);
        set(gca,'units','pixels');
        axis equal
        axis([0,vidWidth,0,vidHeight])
        box off
        set(gca,'XTick',[],'YTick',[]);
        title(sprintf('Loop %i, first, %0.2fs (%s)',iter,times_first{iter}(first(1)),sec2min(times_first{iter}(first(1)))));
        
        subplot(1,2,2);
        image(s_second(first(2)).cdata);
        set(gca,'units','pixels');
        axis equal
        axis([0,vidWidth,0,vidHeight])
        box off
        set(gca,'XTick',[],'YTick',[]);
        title(sprintf('Loop %i, second, %0.2fs (%s)',iter,times_second{iter}(first(2)),sec2min(times_second{iter}(first(2)))));
        
        saveas(gcf,sprintf('Memento_loop_%i_images_first95.png',iter));
        
        close
        
    end
    fprintf('Best correlation is %.4f is at %.2fs and %.2fs (inds %i and %i)\n\n',corrdata(i,j),times_first{iter}(i),times_second{iter}(j),i,j);
    
    top = all_corrdata{iter};
    top(top<0.95)=0;
    handle = plot_matrix_simple(top,sprintf('Memento key-frame %i timing (>0.95)',iter),round(times_second{iter}),round(times_first{iter}));
    set(handle,'Position',[75   275   921   819]);
    saveas(handle,sprintf('Memento_keyframe_%i_matrix_over95.png',iter));
    close(handle);
    
    diary('memento_videoanalysis_add.txt')
    
end

diary off
save('repeating_frames_analysis_add.mat','all_corrdata','times_second','times_first');

%diary off
% image(s(5).cdata)
% set(gcf,'position',[150 150 vidObj.Width vidObj.Height]);
% set(gca,'units','pixels');
% set(gca,'position',[0 0 vidObj.Width vidObj.Height]);
% movie(s,1,vidObj.FrameRate);