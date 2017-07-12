function plot_timeline
%PLOT_TIMELINE Summary of this function goes here
%   Detailed explanation goes here

% X-axis = video time (sujet)
% Y-axis = cinematic time (fabula)
MOVIE_TOTAL = 6350; %1:45:50

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


keyframe_start = [...
    ...%360.00,...
    578.00,...
    919.00,...
    1272.00,...
    1530.00,...
    1818.00,...
    2582.00,...
    2990.00,...
    3078.00,...
    3249.00,...
    3410.00,...
    ...%3599.00,...
    4195.00,... %leikkautuu
    4434.00,...
    4696.00,...
    4970.00,...
    6341.00];

keyframe_end = [...
        ...%366,...
        585,...
        928,...
        1279,...
        1537,...
        1822,...
        2590,...
        2995,...
        3080,...
        3256,...
        3422,...
        ...%3605,...
        4198,...
        4445,...
        4700,...
        4975,...
        6348];
    
loop_start = [...
    171.7   % 2
    402.28  % 3
    639.1   % 4
    1011.6  % 5
    1323.04 % 6
    1889.5  % 7
    2649.16 % 8
    2841.52 % 9
    3009.7  % 10
    3092.6  % 11
    3743.08 % 13 (+0.6 for MEM_5)
    4032.32 % 14 (+0.6 for MEM_5)
    4229.85% 15
    4746    % 16
    5190.48 % 17
    ];
loop_end = loop_start + 5;
loop_duration = loop_end - loop_start;

a=BW_start(match_closest(keyframe_end(1:end-1),BW_start));
keyframe_end = [a,MOVIE_TOTAL];
keyframe_duration = keyframe_end-keyframe_start;
    
STORY_TOTAL = MOVIE_TOTAL-sum(keyframe_duration);
        
BW_duration = BW_end - BW_start;
COLOR_start = [0,BW_end];
COLOR_duration = BW_start(2:end) - BW_end(1:(end-1));
COLOR_duration(end+1)=MOVIE_TOTAL-BW_end(end);

close all;
handle = figure('Position',[2   485   814   609]);hold on;box on;grid on;
axis([0,MOVIE_TOTAL,0,MOVIE_TOTAL]);

leg(1) = line([0,145.8],STORY_TOTAL - [0,BW_start(1)],'Color','r','LineWidth',2);
line(BW_start(1)+[0,0],[0,STORY_TOTAL-BW_start(1)],'Color','k','LineWidth',1,'LineStyle','--');

BW_last = 0;
for i=1:length(BW_start)
    a = line([BW_start(i),BW_end(i)],BW_last +[0,BW_duration(i)],'Color','b','LineWidth',2);
    if i==length(BW_start)
        leg(2) = a;
    end
    BW_points{i}=BW_last +[0,BW_duration(i)];
    
    BW_last = BW_last + BW_duration(i);
    
end

COLOR_last = STORY_TOTAL-BW_start(1);
REPEAT_total = 0;
for i=1:(length(BW_start)-1)            
    
    for k=1:length(keyframe_start)
        if keyframe_start(k)<BW_start(i+1) && keyframe_start(k)>BW_start(i)
            REPEAT_total=REPEAT_total+keyframe_duration(k);            
            k=-k;
            %line([BW_start(i+1),BW_start(i+1)],[0,COLOR_last + REPEAT_total],'Color','k','LineStyle','--');
            break;
        end
    end
    

    
    line([BW_end(i),BW_start(i+1)],(COLOR_last + REPEAT_total + [-COLOR_duration(i),0]),'Color','r','LineWidth',2);
    if k<0
        k=-k;
       line([keyframe_start(k),BW_start(i+1)],[COLOR_last + REPEAT_total - keyframe_duration(k),COLOR_last + REPEAT_total],'Color','g','LineWidth',3);
    end
        
    line(BW_end(i)+[0,0],[BW_points{i}(2),COLOR_last+REPEAT_total-COLOR_duration(i)],'Color','k','LineWidth',1,'LineStyle','--');
    line(BW_start(i+1)+[0,0],[BW_points{i+1}(1),COLOR_last+REPEAT_total],'Color','k','LineWidth',1,'LineStyle','--');

    
    for kk=1:length(loop_start)
        if loop_start(kk)<BW_start(i+1) && loop_start(kk)>BW_start(i)
            if kk==length(loop_start)
                leg(4)=line([BW_end(i),loop_end(kk)],COLOR_last + REPEAT_total-COLOR_duration(i)+[0,loop_end(kk)-BW_end(i)],'Color','c','LineWidth',3);
            else
                line([BW_end(i),loop_end(kk)],COLOR_last + REPEAT_total-COLOR_duration(i)+[0,loop_end(kk)-BW_end(i)],'Color','c','LineWidth',3);
            end
            break;
        end
    end        

    COLOR_last = COLOR_last - COLOR_duration(i);

end

REPEAT_total=REPEAT_total+keyframe_duration(end);

line([BW_end(end),MOVIE_TOTAL],(COLOR_last + REPEAT_total + [-COLOR_duration(end),0]),'Color','r','LineWidth',2);

leg(3) = line([keyframe_start(end),MOVIE_TOTAL],[COLOR_last + REPEAT_total - keyframe_duration(end),COLOR_last + REPEAT_total],'Color','g','LineWidth',3);
%line(MOVIE_TOTAL+[0,0],[0,COLOR_last + REPEAT_total],'Color','k','LineStyle','--');




% for i=1:(length(BW_start)-1)
%     line([BW_end(i),BW_start(i+1)],(COLOR_last +[-COLOR_duration(i),0]),'Color','r','LineWidth',2);
%     COLOR_last = COLOR_last - COLOR_duration(i);
% end
xlabel('Movie time [s]','FontSize',16)
ylabel('Story progression time [s]','FontSize',16)
title('Memento ORIGINAL movie vs. story time')
legend(leg,{'COLOR','BW','KEY-FRAME','LOOP-START'})

annotation(handle,'textbox',...
    [0.519155272901125 0.642292930865887 0.382417358636458 0.110666046062113],...
    'String',{['Total movie duration = ',num2str(MOVIE_TOTAL,'%.0f'),'s'],['Total story duration = ',num2str(MOVIE_TOTAL-REPEAT_total,'%.0f'),'s']},...
    'LineStyle','none',...
    'FontSize',14);


figure;

plot([0,sum(BW_duration)],[0,sum(BW_duration)],'Color','b','LineWidth',2);
hold on;
plot(sum(BW_duration)+[0,STORY_TOTAL-REPEAT_total-sum(BW_duration)],sum(BW_duration)+[0,STORY_TOTAL-REPEAT_total-sum(BW_duration)],'Color','r','LineWidth',2);
xlabel('Movie time [s]','FontSize',16)
ylabel('Story progression time [s]','FontSize',16)
title('Memento CHRONOLOGICAL movie vs. story time')
axis tight;
legend('BW','COLOR')


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
