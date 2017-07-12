function plot_memento_timestructure
%PLOT_TIMELINE Summary of this function goes here
%   Detailed explanation goes here
HOME = pwd;
addpath(HOME);
% X-axis = video time (sujet)
% Y-axis = cinematic time (fabula)
MOVIE_TOTAL_ORIG = [35*60+20,35*60+17,35*60+8];
MOVIE_TOTAL_KRON = [34*60+2,33*60+15,32*60+15];

blue = [0.3,0.75,0.93];
green = [51,204,51]/255;
grey = [153,153,153]/255;

SES_INTERVAL = 60;
BAR_HEIGHT = 0.1;
KEY_HEIGHT = 0.11;
SHIFT = 0.024;

Y_POS_ORIG=0.5;
Y_POS_KRON = 0;

SES_START_ORIG = [0,MOVIE_TOTAL_ORIG(1)+SES_INTERVAL,sum(MOVIE_TOTAL_ORIG(1:2))+2*SES_INTERVAL];
SES_START_KRON = [0,MOVIE_TOTAL_KRON(1)+SES_INTERVAL,sum(MOVIE_TOTAL_KRON(1:2))+2*SES_INTERVAL];

cd('loop_keyframes');
TIMING_file = 'memento_loop_keyframe_timings_FINAL.txt';
TIMING_keys_orig = parse_timing_file(TIMING_file,2);
cd(HOME);

cd('chrono_keyframes');
TIMING_file = 'memento_chrono_keyframes_FINAL.txt';
TIMING_keys_kron = parse_timing_file(TIMING_file,3);
cd(HOME);

cd('bw_scenes');
TIMING_file = 'memento_bw_scene_timings_FINAL.txt';
TIMING_bw_orig = parse_timing_file(TIMING_file,3);
cd(HOME);

cd('chrono_original_segments');
TIMING_file = 'memento_chrono_original_color_timings_FINAL.txt';
TIMING_segments = parse_timing_file(TIMING_file,2);
cd(HOME);

TIMING_bw_kron(1).FIRST_time = [5,23*60+4];
TIMING_bw_kron(1).FIRST_ses = 1;
TIMING_bw_kron(1).ID = 1;
TIMING_bw_kron(1).FIRST_time_str = '00:00:05:00 -> 00:23:04:00';

for j=1:length(TIMING_keys_orig)
    t2 = TIMING_keys_orig(j).SECOND_time;
    ses2=TIMING_keys_orig(j).SECOND_ses;
    for i=1:length(TIMING_segments)
        t1 = TIMING_segments(i).FIRST_time;
        ses1=TIMING_segments(i).FIRST_ses;
        
        
        if ses1==ses2
            
            if diff(t1)<diff(t2)
                error('!!!')
            end
            
            if (t2(1)>t1(1) && t2(1)<t1(2))
                fprintf('ID %i found\n',TIMING_keys_orig(j).ID)
            end
        end
    end
end

%%

MOVIE_TOTAL = sum(MOVIE_TOTAL_ORIG)+2*SES_INTERVAL;

close all;
handle = figure('Position',[2         634        1233         460]);
hold on;
box off;
xlabel('Time [s]','FontSize',16)
%grid on;
axis([-50,MOVIE_TOTAL+50,-0.35,0.85]);

% [Xf, Yf] = ds2nfu([-60,-60],[0,0]+0.01); 
% annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','none','String','Original','HorizontalAlignment','right','FontSize',22,'VerticalAlignment','middle');    
% [Xf, Yf] = ds2nfu([-60,-60],[0,0]+0.01); 
% annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','none','String','Edited','HorizontalAlignment','right','FontSize',22,'VerticalAlignment','middle');    
% 
% set(gca,'YTick',[-SHIFT_CHRONO,0],'YTickLabel',{'',''},'FontSize',15);

[Xf, Yf] = ds2nfu(SES_START_ORIG(1)+[0,MOVIE_TOTAL_ORIG(1)],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','r');    
[Xf, Yf] = ds2nfu(SES_START_ORIG(1)+[0,1*60+52],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','m','facecolor','m');    

[Xf, Yf] = ds2nfu(SES_START_ORIG(2)+[0,MOVIE_TOTAL_ORIG(2)],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','r');    

[Xf, Yf] = ds2nfu(SES_START_ORIG(3)+[0,MOVIE_TOTAL_ORIG(3)],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','r');    

[Xf, Yf] = ds2nfu(SES_START_KRON(1)+[0,MOVIE_TOTAL_KRON(1)],Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','r');    
[Xf, Yf] = ds2nfu(SES_START_KRON(1)+[0,5],Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','m','facecolor','m');    

[Xf, Yf] = ds2nfu(SES_START_KRON(2)+[0,MOVIE_TOTAL_KRON(2)],Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','r');    

[Xf, Yf] = ds2nfu(SES_START_KRON(3)+[0,MOVIE_TOTAL_KRON(3)],Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]);
annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','r');    

a=5;
for i=1:length(TIMING_bw_orig)
        
    ses = TIMING_bw_orig(i).FIRST_ses;
    times =  TIMING_bw_orig(i).FIRST_time;
    
    times(2)=times(2)-1.6;
            
    TIMING_bw_kron(i).FIRST_ses=1;
    TIMING_bw_kron(i).FIRST_time=a+[0,diff(times)];
    a = a + diff(times);
    
    if mod(i,2)==0
        s = SHIFT;
    else        
        s = -SHIFT;
    end
        
    [Xf, Yf] = ds2nfu(SES_START_ORIG(ses)+times,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor',grey,'facecolor',grey);        
    
%     if i>1 && TIMING_bw_orig(i-1).FIRST_ses == TIMING_bw_orig(i).FIRST_ses 
%         [Xf1, Yf] = ds2nfu(SES_START_ORIG(ses)+old_times,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + s); 
%         [Xf2, Yf] = ds2nfu(SES_START_ORIG(ses)+times,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + s);               
%         Xf = [Xf1(2),Xf2(1)];            
%         annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(i-1),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');           
%     end
    old_times = times;

end
% [Xf1, Yf] = ds2nfu(TIMING_bw_orig(7).FIRST_time+1.6,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
% [Xf2, Yf] = ds2nfu(MOVIE_TOTAL_ORIG(2)*[1,1],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + SHIFT);               
% Xf = [Xf1(2),Xf2(1)];            
% annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(7),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');           
% [Xf2, Yf] = ds2nfu(SES_START_ORIG(3)+TIMING_bw_orig(18).FIRST_time+1.6,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
% [Xf1, Yf] = ds2nfu(SES_START_ORIG(3)*[1,1],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + SHIFT);               
% Xf = [Xf1(2),Xf2(1)];            
% annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(17),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');           
% [Xf1, Yf] = ds2nfu(SES_START_ORIG(3)+TIMING_bw_orig(22).FIRST_time+1.6,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
% [Xf2, Yf] = ds2nfu(MOVIE_TOTAL*[1,1],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] - SHIFT);               
% Xf = [Xf1(2),Xf2(1)];            
% annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(22),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');           

for i=1:length(TIMING_bw_kron)
        
    ses = TIMING_bw_kron(i).FIRST_ses;
    times =  TIMING_bw_kron(i).FIRST_time;
    
    [Xf, Yf] = ds2nfu(SES_START_KRON(ses)+times,Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor',grey,'LineWidth',1.5);    

end


for i=1:length(TIMING_keys_kron)
        
    ses = TIMING_keys_kron(i).FIRST_ses;
    times =  TIMING_keys_kron(i).FIRST_time;
    
    [Xf, Yf] = ds2nfu(SES_START_KRON(ses)+times,Y_POS_KRON*[1,1] + [-BAR_HEIGHT - KEY_HEIGHT,BAR_HEIGHT]); 
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','b','facecolor','b','LineWidth',2);             
    annotation(handle,'textbox',[Xf(2),Yf(1)-diff(Yf)/6,diff(Xf),eps],'String',num2str(i),'FontSize',12,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');    

end

for i=1:length(TIMING_keys_orig)
        
    ses = TIMING_keys_orig(i).FIRST_ses;
    times =  TIMING_keys_orig(i).FIRST_time;
    
    [Xf, Yf] = ds2nfu(SES_START_ORIG(ses)+times,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT + KEY_HEIGHT]); 
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','b','facecolor','b','LineWidth',2);    
    annotation(handle,'textbox',[Xf(2),Yf(2)+diff(Yf)/6,diff(Xf),eps],'String',num2str(i),'FontSize',12,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');    
    
    ses = TIMING_keys_orig(i).SECOND_ses;
    times =  TIMING_keys_orig(i).SECOND_time;
    
    [Xf, Yf] = ds2nfu(SES_START_ORIG(ses)+times,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT - KEY_HEIGHT,BAR_HEIGHT]); 
    annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor',green,'facecolor',green,'LineWidth',2);   
    annotation(handle,'textbox',[Xf(2),Yf(1)-diff(Yf)/6,diff(Xf),eps],'String',num2str(i),'FontSize',12,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');       
    
end


TIMING_segments(21).FIRST_time(1)=1495;
TIMING_segments(21).SECOND_time(1) = TIMING_segments(21).SECOND_time(1) +  331.8400;
for i=1:length(TIMING_segments)
        
    ses1 = TIMING_segments(i).FIRST_ses;
    times1 =  TIMING_segments(i).FIRST_time;  
    ID = TIMING_segments(i).ID;
    
    ses2 = TIMING_segments(i).SECOND_ses;
    times2 =  TIMING_segments(i).SECOND_time;    
    
    if mod(i,2)==0
        s = SHIFT;
    else        
        s = -SHIFT;
    end    
    
%     [Xf, Yf] = ds2nfu(SES_START_ORIG(ses1)+times1,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
%     annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','w','facecolor','r');    
%     
     if ID==8   
         [Xf, Yf] = ds2nfu(SES_START_KRON(ses2)+times2,Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
         annotation(handle,'rectangle',[Xf(2),Yf(1),eps,diff(Yf)],'edgecolor','w','facecolor','r');                
     end
     [Xf, Yf] = ds2nfu(SES_START_KRON(ses2)+times2,Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT]); 
     annotation(handle,'rectangle',[Xf(1),Yf(1),eps,diff(Yf)],'edgecolor','w','facecolor','r');       
    
        %[Xf2, Yf] = ds2nfu(SES_START_KRON(ses2)+old_times,Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + s); 
        [Xf, Yf] = ds2nfu(SES_START_KRON(ses2)+times2,Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + s);                          
        annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(ID ),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');           
    
        [Xf, Yf] = ds2nfu(SES_START_ORIG(ses1)+times1,Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + s);                          
        annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(ID ),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');           
    
    old_times = times;    

end
[Xf, Yf] = ds2nfu([SES_START_ORIG(1)+TIMING_bw_orig(7).FIRST_time(2),SES_START_ORIG(1)+MOVIE_TOTAL_ORIG(1)],Y_POS_ORIG*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + s);
annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(7),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');
[Xf, Yf] = ds2nfu(SES_START_KRON(3)+[TIMING_segments(7).SECOND_time(2),TIMING_segments(6).SECOND_time(1)],Y_POS_KRON*[1,1] + [-BAR_HEIGHT,BAR_HEIGHT] + s);
annotation(handle,'textbox',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'String',num2str(7),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');


% 
% a=0;
% for i=1:length(COLOR_start)
%     
%     [Xf, Yf] = ds2nfu([COLOR_start(i);COLOR_end(i)],[-0.02,0.02]); 
%     annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','k','facecolor','r');    
%        
%     if i>1
%         if mod(i,2)==0
%             [Xf, Yf] = ds2nfu((COLOR_start(i)+COLOR_end(i))/2, 0.04);
%         else
%             [Xf, Yf] = ds2nfu((COLOR_start(i)+COLOR_end(i))/2, -0.06);
%         end
%         annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',num2str(i-1),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');
%     end
% 
%     if i>1        
%         [Xf, Yf] = ds2nfu(BW_total+a+[0,COLOR_duration(end-i+2)],[-0.02,0.02]-SHIFT_CHRONO);        
%         annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','k','facecolor','r');
%                         
%         if mod(i,2)==1
%             [Xf, Yf] = ds2nfu((2*BW_total+2*a+COLOR_duration(end-i+2))/2,0.04-SHIFT_CHRONO);
%         else
%             [Xf, Yf] = ds2nfu((2*BW_total+2*a+COLOR_duration(end-i+2))/2,-0.06-SHIFT_CHRONO);
%         end
%         annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',num2str(24-i),'FontSize',14,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle','FontWeight','bold');
%         
%         a = a + COLOR_duration(end-i+2);
%     end
% 
% end
% 
% for i=1:length(keyframe_start)
%     
%     [Xf, Yf] = ds2nfu([keyframe_start(i);keyframe_end(i)],[-0.02,0.02]);
%     annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','g','facecolor','g');
%     
%     plot(keyframe_start(i)*[1,1],[0.0,0.10],'Color','g','LineWidth',2);
%     %line([COLOR_start(i),COLOR_end(i)],[0,0]+0.001,'Color','r','LineWidth',7);
%     
%     [Xf, Yf] = ds2nfu((keyframe_start(i)+keyframe_end(i))/2, 0.13);
%     annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',[num2str(i)],'FontSize',13,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle');
%     
% end
% 
% for i=1:length(loop_start)
%     
%     [Xf, Yf] = ds2nfu([loop_start(i);loop_end(i)],[-0.02,0.02]);
%     annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','g','facecolor','g');
%     
%     plot(loop_start(i)*[1,1],[-0.1,0.00],'Color','c','LineWidth',2);
%     %line([COLOR_start(i),COLOR_end(i)],[0,0]+0.001,'Color','r','LineWidth',7);
%     
%     [Xf, Yf] = ds2nfu((loop_start(i)+loop_end(i))/2, -0.13);
%     annotation(handle,'textbox',[Xf,Yf,eps,eps],'String',[num2str(i)],'FontSize',13,'HorizontalAlignment','center','LineStyle','none','VerticalAlignment','middle');
%     
% end
% 
% 
% for i=1:length(BW_start)
%     
%     [Xf, Yf] = ds2nfu([BW_start(i);BW_end(i)],[-0.018,0.018]); 
%     annotation(handle,'rectangle',[Xf(1),Yf(1),diff(Xf),diff(Yf)],'edgecolor','k','facecolor','b');    
%     
% 
% end
% 
% plot(SESSION_OFFSET_ORIG(2)*[1,1],[-0.08,0.08],'Color','k','LineWidth',2);
% plot(SESSION_OFFSET_ORIG(3)*[1,1],[-0.08,0.08],'Color','k','LineWidth',2);
% 
% plot(2032*[1,1],[-0.08,0.08]-SHIFT_CHRONO,'Color','k','LineWidth',2);
% plot(4140*[1,1],[-0.08,0.08]-SHIFT_CHRONO,'Color','k','LineWidth',2);

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
