clc;
clear variables;
close all;

addpath(pwd);

% initial & keyframe onset
cd('loop_keyframes');
TIMING_file = 'memento_loop_keyframe_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,2);
cd ..

% bw onset
cd('bw_scenes');
TIMING_file = 'memento_bw_scene_timings_FINAL.txt';
TIMING_data_bw = parse_timing_file(TIMING_file,3);
cd ..

% IDs = [];
% for i=1:length(TIMING_data)
%     IDs(end+1)=TIMING_data(i).ID;
% end

key_IDs = [3,4,5,6,7,9,11,12,13,14,17,18,19,21];
initial_IDs = [1,2,3,4,5,7,9,10,11,12,15,16,17,19,21];

k=0;
kk=0;
delays=[];
for i=1:length(TIMING_data_bw)
    for j=1:length(key_IDs)
        if key_IDs(j)==TIMING_data_bw(i).ID
           if TIMING_data_bw(i).FIRST_ses == TIMING_data(j).SECOND_ses
            k=k+1;
            delays(k)=TIMING_data_bw(i).FIRST_time(1) - TIMING_data(j).SECOND_time(2);
            keyframe_duration(k)=diff(TIMING_data(j).SECOND_time);
           end
        end
    end
     for j=1:length(initial_IDs)   
        if initial_IDs(j)==TIMING_data_bw(i).ID
           if TIMING_data_bw(i).FIRST_ses == TIMING_data(j).FIRST_ses
            kk=kk+1;
            initial_delays(kk)=TIMING_data(j).FIRST_time(1) - TIMING_data_bw(i).FIRST_time(2);
            initial_duration(kk)=diff(TIMING_data(j).FIRST_time);
           end
        end
    end
end
delays(k+1)=2; % until stimulus ending
keyframe_duration(k+1)=diff(TIMING_data(15).SECOND_time);

if any(initial_delays<-0.1)
    error('!!!')
end
initial_delays(initial_delays<0)=0;

figure('Position',[596   669   644   425]);
keyframe_segments = [2,3,4,5,6,8,10,11,12,13,16,17,18,20,22];
hBar = bar(1:15,[keyframe_duration;delays]','stacked');hold on;
%plot([0,16],[1,1]*median(keyframe_duration),'k--','LineWidth',1);
plot([0,16],[1,1]*median(keyframe_duration+delays),'k--','LineWidth',1);
axis tight
text(1:15,zeros(1,15)+diff(get(gca,'YLim'))/40,num2cellstr(keyframe_segments),'HorizontalAlignment','center','Color','w','FontSize',14);

set(hBar,{'FaceColor'},{[0,0.45,0.74];'r'});
title('keyframes')
xlabel('Loop ID')
ylabel('Time [s]');
set(gca,'FontSize',16)
legend('Keyframe duration','Delay to BW','location','best')

if max(abs(initial_duration-keyframe_duration))>0.05
    error('!!!')
end

initial_segments = [1,2,3,4,5,7,9,10,11,12,15,16,17,19,21];

figure('Position',[596   669   644   425]);
hBar=bar(1:15,[initial_delays;initial_duration]','stacked');hold on;
plot([0,16],[1,1]*median(initial_delays),'k--','LineWidth',1);
text(1:15,zeros(1,15)+diff(get(gca,'YLim'))/40,num2cellstr(initial_segments),'HorizontalAlignment','center','Color','w','FontSize',14);

set(hBar,{'FaceColor'},{'r';[0,0.45,0.74]});
title('initial frames')
xlabel('Loop ID')
ylabel('Time [s]');
set(gca,'FontSize',16)
axis tight
legend('Delay to initial','Initial duration','location','best')



    