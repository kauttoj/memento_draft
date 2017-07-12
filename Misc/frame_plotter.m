function [FRAME1,FRAME2]=frame_plotter(TIMING_FILE,VIDEO_FILES,PREFIX,TYPE)

% TYPE=1 with two different versions (each element has two rows)
% TYPE=2 with one version

if nargin<4
    TYPE=2;
end

FRAME1=[];
FRAME2=[];

MOVIE_SEGMENTS = parse_timing_file(TIMING_FILE,TYPE);

if isfield(VIDEO_FILES,'SECOND')
    TYPE=1;
    for i=1:length(VIDEO_FILES.FIRST)
        vidObj_left{i} = VideoReader(VIDEO_FILES.FIRST{i});
    end
    for i=1:length(VIDEO_FILES.SECOND)
        vidObj_right{i} = VideoReader(VIDEO_FILES.SECOND{i});
    end
else
    TYPE=2;
    for i=1:length(VIDEO_FILES.FIRST)
        vidObj_left{i} = VideoReader(VIDEO_FILES.FIRST{i});
    end
end

for i=1:length(MOVIE_SEGMENTS)
    
    if TYPE==1
        
        vidObj_first = vidObj_left{MOVIE_SEGMENTS(i).FIRST_ses};
        vidObj_second = vidObj_right{MOVIE_SEGMENTS(i).SECOND_ses};
        
        HEIGHT_first = vidObj_first.Height;
        WIDTH_first = vidObj_first.Width;
        
        HEIGHT_second = vidObj_second.Height;
        WIDTH_second = vidObj_second.Width;
        
        HEIGHT=max(HEIGHT_first,HEIGHT_second);
        WIDTH=WIDTH_first+WIDTH_second;
        
        figure;
        set(gcf,'position',[10 10 WIDTH+20 HEIGHT+30],'units','pixels');
        
        subplot(1,2,1);
        vidObj_first.CurrentTime = MOVIE_SEGMENTS(i).FIRST_time(1);
        TARGET = readFrame(vidObj_first);
        image(TARGET);
        FRAME1{i}=TARGET;
        set(gca,'units','pixels','Position',[5,0,WIDTH_first-70,HEIGHT]);
        box off
        set(gca,'XTick',[],'YTick',[]);
        title(sprintf('Item %i, %0.2fs (%s)',MOVIE_SEGMENTS(i).ID,MOVIE_SEGMENTS(i).FIRST_time(1),sec2min(MOVIE_SEGMENTS(i).FIRST_time(1))));
        
        subplot(1,2,2);
        vidObj_second.CurrentTime = MOVIE_SEGMENTS(i).SECOND_time(1);
        TARGET = readFrame(vidObj_second);
        image(TARGET);
        FRAME2{i}=TARGET;
        set(gca,'units','pixels','Position',[WIDTH_first+10-70,0,WIDTH_second-70,HEIGHT]);
        box off
        set(gca,'XTick',[],'YTick',[]);
        title(sprintf('Item %i, %0.2fs (%s)',MOVIE_SEGMENTS(i).ID,MOVIE_SEGMENTS(i).SECOND_time(1),sec2min(MOVIE_SEGMENTS(i).SECOND_time(1))));
        
        set(gcf, 'PaperPositionMode','auto')
        saveas(gcf,sprintf('%s_ID_%i.png',PREFIX,MOVIE_SEGMENTS(i).ID));
        
    else
        
        vidObj_first = vidObj_left{MOVIE_SEGMENTS(i).FIRST_ses};
        
        HEIGHT_first = vidObj_first.Height;
        WIDTH_first = vidObj_first.Width;
        
        
        HEIGHT=HEIGHT_first;
        WIDTH=WIDTH_first;
        
        figure;
        set(gcf,'position',[10 10 WIDTH+20 HEIGHT+30],'units','pixels');
        
        vidObj_first.CurrentTime =  MOVIE_SEGMENTS(i).FIRST_time(1);
        TARGET = readFrame(vidObj_first);
        image(TARGET);
        FRAME1{i}=TARGET;
        set(gca,'units','pixels','Position',[5,0,WIDTH_first,HEIGHT]);
        box off
        set(gca,'XTick',[],'YTick',[]);
        title(sprintf('Item %i, %0.2fs (%s)',MOVIE_SEGMENTS(i).ID,MOVIE_SEGMENTS(i).FIRST_time(1),sec2min(MOVIE_SEGMENTS(i).FIRST_time(1))));
               
        set(gcf, 'PaperPositionMode','auto')
        saveas(gcf,sprintf('%s_ID_%i.png',PREFIX,MOVIE_SEGMENTS(i).ID));
        
    end
    
end

end

