function plot_shifted_keyframes
%PLOT_TIMELINE Summary of this function goes here
%   Detailed explanation goes here
HOME = pwd;
addpath(HOME);
% X-axis = video time (sujet)
% Y-axis = cinematic time (fabula)
MOVIE_TOTAL = 6350; %1:45:50

SHIFT_CHRONO = 0.3;

% Below timings are for file "MEMENTO.m4v"
BW_start = [...
145.8,...      % 00:02:25.800	
366.11,...      % 00:06:06.110	
585.685,...    % 00:09:45.685	
930.548,...     % 00:15:30.848	
1281.06,...    % 00:21:21.060	
1536.984,...		% 00:25:36.984	
1822.554,...		% 00:30:22.554	
2268.304,...		% 00:37:48.304	
2589.57,...     % 00:43:09.570	
2813.394,...	% 00:46:53.394		
2995.476,...	% 00:49:55.476	
3080.276,...	% 00:51:20.276		
3256.649,...	% 00:54:16.649		
3422.131,...	% 00:57:02.131		
3605.311,...		% 01:00:05.311	
3989.571,...		% 01:06:29.571	
4197.899,...		% 01:09:57.899	
4445.671,...		% 01:14:05.671
4701.495,...		% 01:18:21.495	
4790.645,...		% 01:19:50.645
4978.362,...	% 01:22:58.362	
5405.269];		% 01:30:05.269	

BW_end = [...
168.4,...	% 00:02:48.400	
402.01,...	% 00:06:42.010	
636.585,...	% 00:10:36.585	
974.348,...	% 00:16:14.348	
1317.46,...	% 00:21:57.460	
1638.384,... % 00:27:18.384	
1873.254,... % 00:31:13.254
2312.904,... % 00:38:32.904	
2648.97,...	% 00:44:08.970	
2841.394,... % 00:47:21.394	
3008.276,...	% 00:50:08.276	
3092.276,... % 00:51:32.276	
3306.049,... % 00:55:06.049	
3461.031,... % 00:57:41.031	
3735.611,... % 01:02:15.611	
4022.871,... % 01:07:02.871	
4222.899,... % 01:10:22.899	
4478.871,... %01:14:38.871
4745.795,... %01:19:05.795
4812.445,... %01:20:12.445
5187.162,...	%01:26:27.162
5738.369]; %01:35:38.369

cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/loop_keyframes');
TIMING_file = 'memento_loop_keyframe_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,2);
cd(HOME);

SESSION_OFFSET_CHRONO = [0,34*60+2,34*60+2+33*60+15];
SESSION_OFFSET_ORIG = [0,2125.5,4243.04];
for i=1:length(TIMING_data)
    keyframe_start(i)=TIMING_data(i).SECOND_time(1) + SESSION_OFFSET_ORIG(TIMING_data(i).SECOND_ses);
    loop_start(i)=TIMING_data(i).FIRST_time(1) + SESSION_OFFSET_ORIG(TIMING_data(i).FIRST_ses);
    keyframe_end(i)=TIMING_data(i).SECOND_time(2) + SESSION_OFFSET_ORIG(TIMING_data(i).SECOND_ses);
    loop_end(i)=TIMING_data(i).FIRST_time(2) + SESSION_OFFSET_ORIG(TIMING_data(i).FIRST_ses);                    
end

keyframe_start=keyframe_start+40;
keyframe_end = keyframe_end+40;

% keyframe_start = [...
%     ...%360.00,...
%     578.00,...
%     919.00,...
%     1272.00,...
%     1530.00,...
%     1818.00,...
%     2582.00,...
%     2990.00,...
%     3078.00,...
%     3249.00,...
%     3410.00,...
%     ...%3599.00,...
%     4195.00,... %leikkautuu
%     4434.00,...
%     4696.00,...
%     4970.00,...
%     6341.00];
% 
% keyframe_end = [...
%         ...%366,...
%         585,...
%         928,...
%         1279,...
%         1537,...
%         1822,...
%         2590,...
%         2995,...
%         3080,...
%         3256,...
%         3422,...
%         ...%3605,...
%         4198,...
%         4445,...
%         4700,...
%         4975,...
%         6348];
%     
% loop_start = [...
%     171.7   % 2
%     402.28  % 3
%     639.1   % 4
%     1011.6  % 5
%     1323.04 % 6
%     1889.5  % 7
%     2649.16 % 8
%     2841.52 % 9
%     3009.7  % 10
%     3092.6  % 11
%     3743.08 % 13 (+0.6 for MEM_5)
%     4032.32 % 14 (+0.6 for MEM_5)
%     4229.85% 15
%     4746    % 16
%     5190.48 % 17
%     ];
% loop_end = loop_start + 5;
% loop_duration = loop_end - loop_start;

keyframe_duration = keyframe_end-keyframe_start;
    
STORY_TOTAL = MOVIE_TOTAL-sum(keyframe_duration);
        
BW_duration = BW_end - BW_start;
COLOR_start = [0,BW_end];
COLOR_duration = BW_start(2:end) - BW_end(1:(end-1));
COLOR_duration(end+1)=MOVIE_TOTAL-BW_end(end);
COLOR_end = [BW_start,MOVIE_TOTAL];

BW_total = sum(BW_duration);

close all;
handle = figure('Position',[2         634        1233         460]);
hold on;
box off;
xlabel('Time [s]','FontSize',16)
%grid on;
axis([-50,MOVIE_TOTAL+55,-0.2,0.19]);

[Xf, Yf] = ds2nfu([-60,-60],[0,0]+0.01); 
annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','none','String','Original','HorizontalAlignment','right','FontSize',22,'VerticalAlignment','middle');    

set(gca,'YTick',[-SHIFT_CHRONO,0],'YTickLabel',{'',''},'FontSize',15);

% %leg(1) = line([0,145.8],[0,0],'Color','r','LineWidth',7);
% [Xf, Yf] = ds2nfu([0;BW_total],[-0.018,0.018]-1); 
% annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','b','facecolor','b');    

COLOR_last = STORY_TOTAL-BW_start(1);
REPEAT_total = 0;
COLOR_null = COLOR_end(1)-COLOR_start(1);

% a=0;
% for i=1:length(BW_start)
%         
%     [Xf, Yf] = ds2nfu(a+[0;BW_duration(i)],[-0.018,0.018]-SHIFT_CHRONO); 
%     annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','b');    
%     
%     a = a + BW_duration(i);
% end

a=0;
for i=1:length(COLOR_start)
    
    [Xf, Yf] = ds2nfu([COLOR_start(i);COLOR_end(i)],[-0.02,0.02]); 
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','k','facecolor','r');    
       
    if i>1
        if mod(i,2)==0
            [Xf, Yf] = ds2nfu((COLOR_start(i)+COLOR_end(i))/2, 0.04);
        else
            [Xf, Yf] = ds2nfu((COLOR_start(i)+COLOR_end(i))/2, -0.06);
        end
        annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',num2str(i-1),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');
    end

    if i>1        
%         [Xf, Yf] = ds2nfu(BW_total+a+[0,COLOR_duration(end-i+2)],[-0.02,0.02]-SHIFT_CHRONO);        
%         annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','k','facecolor','r');
%                         
%         if mod(i,2)==1
%             [Xf, Yf] = ds2nfu((2*BW_total+2*a+COLOR_duration(end-i+2))/2,0.04-SHIFT_CHRONO);
%         else
%             [Xf, Yf] = ds2nfu((2*BW_total+2*a+COLOR_duration(end-i+2))/2,-0.06-SHIFT_CHRONO);
%         end
%         annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',num2str(24-i),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');
        
        a = a + COLOR_duration(end-i+2);
    end

end

for i=1:length(keyframe_start)
    
    [Xf, Yf] = ds2nfu([keyframe_start(i);keyframe_end(i)],[-0.02,0.02]);
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','g','facecolor','g');
    
    plot(keyframe_start(i)*[1,1],[0.0,0.10],'Color','g','LineWidth',2);
    %line([COLOR_start(i),COLOR_end(i)],[0,0]+0.001,'Color','r','LineWidth',7);
    
    [Xf, Yf] = ds2nfu((keyframe_start(i)+keyframe_end(i))/2, 0.13);
    annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',[num2str(i)],'FontSize',13,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle');
    
end

% for i=1:length(loop_start)
%     
%     [Xf, Yf] = ds2nfu([loop_start(i);loop_end(i)],[-0.02,0.02]);
%     annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','g','facecolor','g');
%     
%     plot(loop_start(i)*[1,1],[-0.1,0.00],'Color','c','LineWidth',2);
%     line([COLOR_start(i),COLOR_end(i)],[0,0]+0.001,'Color','r','LineWidth',7);
%     
%     [Xf, Yf] = ds2nfu((loop_start(i)+loop_end(i))/2, -0.13);
%     annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',[num2str(i)],'FontSize',13,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle');
%     
% end


for i=1:length(BW_start)
    
    [Xf, Yf] = ds2nfu([BW_start(i);BW_end(i)],[-0.018,0.018]); 
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','k','facecolor','b');    
    

end

plot(SESSION_OFFSET_ORIG(2)*[1,1],[-0.08,0.08],'Color','k','LineWidth',2);
plot(SESSION_OFFSET_ORIG(3)*[1,1],[-0.08,0.08],'Color','k','LineWidth',2);

plot(2032*[1,1],[-0.08,0.08]-SHIFT_CHRONO,'Color','k','LineWidth',2);
plot(4140*[1,1],[-0.08,0.08]-SHIFT_CHRONO,'Color','k','LineWidth',2);

end

function ind=match_closest(keyframe_end,BW_time)

max_offset = 0;
    for i=1:length(keyframe_end)
        [offset,k]=min(abs(keyframe_end(i)-BW_time));
        if offset > max_offset
            max_offset = offset;
        end
        if offset>4
            offset,keyframe_end(i)
            error('Too large offset!');            
        end
        ind(i)=k;        
    end
    
    if length(unique(ind))<length(ind)
        error('wrong size!');
    end        
    
    fprintf('max offset was %s seconds\n',max_offset);

end
