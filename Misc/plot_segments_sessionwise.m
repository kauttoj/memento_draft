clc;
clear variables;
close all;

HOME = pwd;

kk=0;

TIMING_file = 'bw_scenes/memento_bw_scene_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,3);

figure;
for i=1:length(TIMING_data)
    
    subplot(3,1,TIMING_data(i).FIRST_ses);hold on;
    if i==1
        kk=kk+1;
        lin(kk)=plot([TIMING_data(i).FIRST_time(1)-eps,TIMING_data(i).FIRST_time,TIMING_data(i).FIRST_time(2)+eps],[0,1,1,0],'b');
    else
        plot([TIMING_data(i).FIRST_time(1)-eps,TIMING_data(i).FIRST_time,TIMING_data(i).FIRST_time(2)+eps],[0,1,1,0],'b');
    end
    box on;
    
end

TIMING_file = 'loop_keyframes/memento_loop_keyframe_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,2);

%cd('chrono_keyframes');
TIMING_file = 'chrono_keyframes/memento_chrono_keyframes_FINAL.txt';
TIMING_keys_kron = parse_timing_file(TIMING_file,3);
%cd(HOME);

clear tempdata time_matrix label;
k=0;
for i=1:length(TIMING_data)
    k=k+1;
    tempdata(k).SES=TIMING_data(i).FIRST_ses;
    tempdata(k).ID=TIMING_data(i).ID;
    tempdata(k).type='I';
    tempdata(k).time=TIMING_data(i).FIRST_time(1);
    label{k}=[tempdata(k).type,num2str(tempdata(k).ID)];
    k=k+1;
    tempdata(k).SES=TIMING_data(i).SECOND_ses;
    tempdata(k).ID=TIMING_data(i).ID;
    tempdata(k).type='K';    
    tempdata(k).time=TIMING_data(i).SECOND_time(1);
    label{k}=[tempdata(k).type,num2str(tempdata(k).ID)];
end

time_matrix = zeros(length( tempdata),length( tempdata));
for i=1:length(tempdata)
    for j=1:length(tempdata)
        if i~=j
            if tempdata(i).SES == tempdata(j).SES
                time_matrix(i,j)=abs(tempdata(i).time - tempdata(j).time);
            end
        end
    end
end

ind = [1:2:30,2:2:30];
time_matrix=time_matrix(ind,ind);
label = label(ind);

plot_matrix_simple(round(time_matrix),'',label,label);hold on;
a=time_matrix(16:end,16:end);
a=triu(a,1);a=a(:);a=a(a>0);[min(a),median(a),max(a)]
plot(15.5*[1,1],[0,30],'k--')
plot([0,30],15.5*[1,1],'k--')

times = 0:0.01:36*60;
envelopes{1}=0*times;
envelopes{2}=0*times;
envelopes{3}=0*times;

for i=1:length(TIMING_data)
    
    a=find(TIMING_data(i).FIRST_time(1)<times,1,'first'):find(TIMING_data(i).FIRST_time(2)>times,1,'last');
    envelopes{TIMING_data(i).FIRST_ses}(a)=1;
    a=find(TIMING_data(i).SECOND_time(1)<times,1,'first'):find(TIMING_data(i).SECOND_time(2)>times,1,'last');
    envelopes{TIMING_data(i).SECOND_ses}(a)=1;
    
    subplot(3,1,TIMING_data(i).FIRST_ses);hold on;
    
    if i==1,
        kk=kk+1;
        lin(kk)=plot([TIMING_data(i).FIRST_time(1)-eps,TIMING_data(i).FIRST_time,TIMING_data(i).FIRST_time(2)+eps],[0,0.95,0.95,0],'g');
    else
        plot([TIMING_data(i).FIRST_time(1)-eps,TIMING_data(i).FIRST_time,TIMING_data(i).FIRST_time(2)+eps],[0,0.95,0.95,0],'g');
    end
    
    subplot(3,1,TIMING_data(i).SECOND_ses);hold on;   
    if i==1,
        kk=kk+1;
        lin(kk)=plot([TIMING_data(i).SECOND_time(1)-eps,TIMING_data(i).SECOND_time,TIMING_data(i).SECOND_time(2)+eps],[0,1.05,1.05,0],'c');
    else        
        plot([TIMING_data(i).SECOND_time(1)-eps,TIMING_data(i).SECOND_time,TIMING_data(i).SECOND_time(2)+eps],[0,1.05,1.05,0],'c');
    end
    
end

nulltimes{1} = [110,310,764,1130,13,2003]+5;
nulltimes{2} = [99,252,605,1410,1793]+5;
nulltimes{3} = [47,315,525,1097,1551,1829]+5;

t = linspace(0,35*60,1000);
for ses=1:3,
    subplot(3,1,ses);
    plot(t,0*t,'k');
    
    for i=1:length(nulltimes{ses})
        a=nulltimes{ses}(i)+[0,4];
        plot([a(1)-eps,a,a(2)+eps],[0,0.90,0.9,0],'m');
        
        if ses==1 && i==1
            kk=kk+1;
            lin(kk)=plot([a(1)-eps,a,a(2)+eps],[0,0.90,0.9,0],'m');
        else
            plot([a(1)-eps,a,a(2)+eps],[0,0.90,0.9,0],'m');
        end
    end
    
    box on;
    axis tight;
    set(gca,'YLim',[-0.05,1.05]);
    title(sprintf('Session %i',ses));
end
xlabel('Time [s]','FontSize',14);
legend(lin,'BW segment','Initial-frame','Key-frame','Null-frame')

figure;

for ses=1:3
    
    duration{ses}=[];
    s=find(envelopes{ses}>0,1,'first');
    l=find(envelopes{ses}>0,1,'last');
    k=0;
    for i=s:l
        if envelopes{ses}(i)==0
            k=k+1;
        else
            if k>0
                duration{ses}(end+1)=k*0.01;
                k=0;
            end
        end
    end
    subplot(3,1,ses);
    bar(duration{ses},'r');
    hold on;
    plot([0.5,length(duration{ses})+0.5],[1,1]*mean(duration{ses}),'k--');
    axis tight;
    text(1:length(duration{ses}),0*duration{ses}+max(duration{ses})/20,num2cellstr(round(duration{ses})),'HorizontalAlignment','center','Color','k','FontSize',12);
    xlabel('Interval')
    ylabel('Duration [s]')
end

%%


%cd('chrono_keyframes');
TIMING_file = 'chrono_keyframes/memento_chrono_keyframes_FINAL.txt';
TIMING_keys_kron = parse_timing_file(TIMING_file,3);
%cd(HOME);

clear tempdata time_matrix label;
k=0;
for i=1:length(TIMING_data)
    k=k+1;
    tempdata(k).SES=TIMING_data(i).FIRST_ses;
    tempdata(k).ID=k;
    tempdata(k).type='I';
    tempdata(k).time=TIMING_data(i).FIRST_time(1);
    label{k}=[tempdata(k).type,num2str(tempdata(k).ID)];
end

time_matrix = zeros(length( tempdata),length( tempdata));
for i=1:length(tempdata)
    for j=1:length(tempdata)
        if i~=j
            if tempdata(i).SES == tempdata(j).SES
                time_matrix(i,j)=abs(tempdata(i).time - tempdata(j).time);
            end
        end
    end
end

plot_matrix_simple(round(time_matrix),'',label,label);hold on;
a=triu(time_matrix,1);a=a(:);a=a(a>0);[min(a),median(a),max(a)]
