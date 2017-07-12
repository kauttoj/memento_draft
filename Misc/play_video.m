function correlation = play_video(FILES,TIMES,DURATION,RATEMOD)
%PLAY_VIDEO Summary of this function goes here
%   Detailed explanation goes here

if nargin<4
    RATEMOD=1;
end

N = length(FILES);

if N>4
    error('!!')        
end

for i=1:N
    
    if ~exist(FILES{i},'file')
        error('video not found')
    end
    videoFReader{i} = VideoReader(FILES{i});
    
    videoFReader{i}.CurrentTime = TIMES(i);
    
    vidWidth(i) = videoFReader{i}.Width;
    vidHeight(i) = videoFReader{i}.Height;
       
    frameRate(i)=videoFReader{i}.FrameRate;    
    
end

if length(unique(frameRate))>1
    error('multiple framerates!')
end
    
scrsz = get(groot,'ScreenSize');

displacement = [0,0;1,0;1,0;1,1];

for i=1:N    
    pos = [displacement(i,1)*scrsz(3)/2,displacement(i,2)*scrsz(4)/2,scrsz(3)/2,scrsz(4)/2];    
    handle{i} = figure('Position',pos);
    handle_ax{i}=gca;    
    set(handle_ax{i},'XTick',[],'YTick',[],'box','on');   
    tithand{i} = title('TEST');
end

pausetime = RATEMOD*1/frameRate(1);

hand = cell(1,N);

ind = 1:N;

CUTTED_PIXELS=0;
if length(unique(vidWidth))>1
    CUTTED_PIXELS=1;
end

k=1;
tot = 0;
while tot<=DURATION
    try
        if mod(k,2)==0
            ind=fliplr(ind);
        end
        for i=ind
            data{i} = readFrame(videoFReader{i});
            delete(hand{i});
            hand{i} = image(data{i},'Parent',handle_ax{i});
            set(tithand{i},'String',[sec2min(tot+TIMES(i)),' (',num2str(round(tot+TIMES(i))),'s)']);
        end
        
        correlation(k)=corr2(squeeze(data{1}(:,(1:(end-CUTTED_PIXELS)),1)),squeeze(data{2}(:,:,1)));
        tot = tot + (1/frameRate(1));
        pause(pausetime);
        k = k+1;
    catch err
        warning('Stopped at frame %i: %s',k,err.message);
        break;
    end
end
for i=1:N
   videoFReader{i}.CurrentTime = TIMES(i);
end

% tic;
% tot = 0;
% k = 1;
% while tot<=DURATION
%     try
%         for i=1:N
%             timepoint=videoFReader{i}.CurrentTime;
%             data = readFrame(videoFReader{i});
%             delete(hand{i});
%             hand{i} = image(data,'Parent',handle_ax{i});
%             title(handle_ax{i},sprintf('frame %i: %s (%.2fs)',k,sec2min(timepoint),timepoint));
%         end
%         pause(pausetime);
%         tot = tot + (1/frameRate(1));
%         k = k+1;
%     catch err
%         warning('Stopped at frame %i: %s',k,err.message);
%         break;
%     end
% end
% 
% a = toc;
% fprintf('Loop time was %.2fs, playrate was %.3f (should be 1)\n',tot,a/k);

end

function axReturn = newplot(hsave)

axReturn = hsave;

end

