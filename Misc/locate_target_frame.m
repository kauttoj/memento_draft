function [ allcorrs,alltimes ] = locate_target_frame(TARGET,video,interval)
%LOCATE_TARGET_FRAME Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    interval = [0,inf];
end

fprintf('Reading video....');
%vidObj = VideoReader('/media/MyBook/memento/PresentationFiles/Stimulus/MEMENTO_clip.avi'); %'/media/MyBook/memento/stimulus/MEMENTO.m4v');
vidObj = VideoReader(video);
fprintf(' done\n');

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
N = vidHeight*vidWidth;

TARGET=TARGET(1:vidHeight,1:vidWidth,:);

a1 = image_zscore(single(squeeze(TARGET(:,:,1))));
a2 = image_zscore(single(squeeze(TARGET(:,:,2))));
a3 = image_zscore(single(squeeze(TARGET(:,:,3))));

vidObj.CurrentTime = interval(1);

print_step = 1*60;
current_timepoint = interval(1) + print_step;

fprintf('...starting time %s\n',sec2min(interval(1)));

k = 1;
success=false;
while 1
    try
        alltimes(k)=vidObj.CurrentTime;
        a = readFrame(vidObj);
        b1 = image_zscore(single(squeeze(a(:,:,1))));
        b2 = image_zscore(single(squeeze(a(:,:,2))));
        b3 = image_zscore(single(squeeze(a(:,:,3))));
        val = sum(sum(a1.*b1 + a2.*b2 + a3.*b3))/(3*(N-1));        
        allcorrs(k)=val;
        
        if val>0.95 && ~success
            success=true;
            fprintf('...first >0.95 frame found! Timepoint %.2fs (%s)\n',vidObj.CurrentTime,sec2min(vidObj.CurrentTime))
        end                        
        if ~hasFrame(vidObj) || alltimes(k)>interval(2)
            break;
        end        
        if alltimes(k) > current_timepoint
            fprintf('...time %s\n',sec2min(alltimes(k)));
            current_timepoint = current_timepoint + print_step;
        end                
        k = k+1;
    catch err
        warning('iter %i failed: %s',k,err.message);
    end
end

[~,i] = max(allcorrs);
fprintf('finished! Best match %.3f was at %.2fs (%s)\n\n',allcorrs(i),alltimes(i),sec2min(alltimes(i)));

end

