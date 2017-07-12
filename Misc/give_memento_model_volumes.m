function [fmri_volumes,handles] = give_memento_model_volumes(movie_offset_times,movie_event_durations,keepold)

if nargin<3
    keepold=0;
end
%GIVE_MODEL_VOLUMES Summary of this function goes here
%   Detailed explanation goes here

fmri_volumes = [];
TR=1.56;

% first pulse is always ~0.18 seconds behind video start
% ALSO 
% slice-timing correction adds another delay (currently half-TR)
first_volume_movie_difference = 0.19 + TR/2;

ses1_offset=0;
ses2_offset=2125.5; % pitäisi olla 2125.16
ses3_offset=4243.04;% % pitäisi olla 4242.24
total_movie_time = 6350;

ses1_vols = 1369;
ses2_vols = 1369;
ses3_vols = 1369;

ses1_scanner_movie_times = ((1:ses1_vols)-1)*TR;
ses2_scanner_movie_times = ((1:ses2_vols)-1)*TR;
ses3_scanner_movie_times = ((1:ses3_vols)-1)*TR;

ses1_scanner_times = ((1:ses1_vols)-1)*TR + ses1_offset + first_volume_movie_difference;
ses2_scanner_times = ((1:ses2_vols)-1)*TR + ses2_offset + first_volume_movie_difference;
ses3_scanner_times = ((1:ses3_vols)-1)*TR + ses3_offset + first_volume_movie_difference;

if length(movie_event_durations)==1 && length(movie_offset_times)>1
    a = ones(1,length(movie_offset_times));
    movie_event_durations = a*movie_event_durations;
    clear a;
end

dt = 0.01;
movie_time = 0:dt:ses3_scanner_times(end);
movie_envelope = 0*movie_time;
all_inds=[];
for i=1:length(movie_offset_times)
   if movie_offset_times(i)==ses2_scanner_times(1)
       movie_offset_times(i)=movie_offset_times(i)+dt;
       warning('Offset shifted (at boundary)')
   end
   if movie_offset_times(i)==ses3_scanner_times(1)
       movie_offset_times(i)=movie_offset_times(i)+dt;
       warning('Offset shifted (at boundary)')
   end          
       
   [~,k]=min(abs( movie_offset_times(i)-movie_time));
   kk = round(movie_event_durations(i)/dt);
   ind = k:min(k+kk,length(movie_envelope));
   movie_envelope(ind)=1;
   if nnz(ismember(ind,all_inds))>0
       error('Envelope overlap detected!')
   end
   all_inds = [all_inds,ind];
end
clear all_inds

[~,k]=min(abs( ses2_scanner_times(1)-movie_time));
if movie_envelope(k)==1
    c=0;
    while movie_envelope(k)==1
        movie_envelope(k)=0;
        k=k+1;
        c=c+1;
    end
    warning('Boundary leak (ses2)! Deleted %f seconds',c*dt);
end

[~,k]=min(abs( ses3_scanner_times(1)-movie_time));
if movie_envelope(k)==1
    c=0;
    while movie_envelope(k)==1
        movie_envelope(k)=0;
        k=k+1;
        c=c+1;
    end
    warning('Boundary leak (ses3)! Deleted %f seconds',c*dt);
end

hrf = spm_hrf(dt);

[~,k1]=min(abs( ses1_scanner_times(end)-movie_time));
[~,k11]=min(abs(ses2_offset-movie_time));
[~,k2]=min(abs( ses2_scanner_times(1)-movie_time));
[~,k3]=min(abs( ses2_scanner_times(end)-movie_time));
[~,k33]=min(abs(ses3_offset-movie_time));
[~,k4]=min(abs( ses3_scanner_times(1)-movie_time));
[~,k5]=min(abs( ses3_scanner_times(end)-movie_time));
[~,k55]=min(abs(total_movie_time-movie_time));
ses1_envelope = movie_envelope(1:k1);
ses1_envelope(k11:end)=0;
ses2_envelope = movie_envelope(k2:k3);
ses2_envelope(k33:end)=0;
ses3_envelope = movie_envelope(k4:k5);
ses3_envelope(k55:end)=0;

ses1_envelope_time = movie_time(1:k1);
ses1_envelope_time(1)=ses1_envelope_time(1)-1e-6;
ses1_envelope_time(end)=ses1_envelope_time(end)+1e-6;
ses2_envelope_time = movie_time(k2:k3);
ses2_envelope_time(1)=ses2_envelope_time(1)-1e-6;
ses2_envelope_time(end)=ses2_envelope_time(end)+1e-6;
ses3_envelope_time = movie_time(k4:k5);
ses3_envelope_time(1)=ses3_envelope_time(1)-1e-6;
ses3_envelope_time(end)=ses3_envelope_time(end)+1e-6;

ses1_model = conv(hrf',ses1_envelope);
ses2_model = conv(hrf',ses2_envelope);
ses3_model = conv(hrf',ses3_envelope);

ses1_model=ses1_model(1:length(ses1_envelope));
ses2_model=ses2_model(1:length(ses2_envelope));
ses3_model=ses3_model(1:length(ses3_envelope));

if ses1_envelope_time(1)>ses1_scanner_times(1) || ses1_envelope_time(end)<ses1_scanner_times(end)
    error('Bad ranges! (ses1)')
end
ses1_fmri_model = interp1(ses1_envelope_time,ses1_model,ses1_scanner_times);
if ses2_envelope_time(1)>ses2_scanner_times(1) || ses2_envelope_time(end)<ses2_scanner_times(end)
    error('Bad ranges! (ses2)')
end
ses2_fmri_model = interp1(ses2_envelope_time,ses2_model,ses2_scanner_times);
if ses3_envelope_time(1)>ses3_scanner_times(1) || ses3_envelope_time(end)<ses3_scanner_times(end)
    error('Bad ranges! (ses3)')
end
ses3_fmri_model = interp1(ses3_envelope_time,ses3_model,ses3_scanner_times);

if nnz(isnan(ses1_fmri_model))>0
    error('bad model values! (ses1)')
end
if nnz(isnan(ses2_fmri_model))>0
    error('bad model values! (ses2)')
end
if nnz(isnan(ses3_fmri_model))>0
    error('bad model values! (ses3)')
end

fmri_volumes{1,1}=ses1_scanner_times-ses1_offset;
fmri_volumes{1,2}=ses1_fmri_model;
%fmri_volumes{1,3}=ses1_envelope;
if max(fmri_volumes{1,2})>1e-6
   fmri_volumes{1,2}=fmri_volumes{1,2}/max(fmri_volumes{1,2}); 
end
fmri_volumes{2,1}=ses2_scanner_times-ses2_offset;
fmri_volumes{2,2}=ses2_fmri_model;
%fmri_volumes{2,3}=ses2_envelope;
if max(fmri_volumes{2,2})>1e-6
   fmri_volumes{2,2}=fmri_volumes{2,2}/max(fmri_volumes{2,2}); 
end
fmri_volumes{3,1}=ses3_scanner_times-ses3_offset;
fmri_volumes{3,2}=ses3_fmri_model;
%fmri_volumes{3,3}=ses3_envelope;
if max(fmri_volumes{3,2})>1e-6
   fmri_volumes{3,2}=fmri_volumes{3,2}/max(fmri_volumes{3,2}); 
end

col='r';
%close all;
if keepold==1
    hold on;
    col='g';
else
    figure;
end

handles(1)=subplot(3,1,1);hold on;
plot(fmri_volumes{1,1},fmri_volumes{1,2},col,'linewidth',2);
plot(ses1_envelope_time-ses1_offset,ses1_envelope,'b');
axis tight;
ylim([-0.2,1.1]);
box on;
title('Session 1');
xlabel('Time [s]')
ylabel('Model and fmri response')

handles(2)=subplot(3,1,2);hold on;
plot(fmri_volumes{2,1},fmri_volumes{2,2},col,'linewidth',2);
plot(ses2_envelope_time-ses2_offset,ses2_envelope,'b');
axis tight;
ylim([-0.2,1.1]);
box on;
title('Session 2');
xlabel('Time [s]')
ylabel('Model and fmri response')

handles(3)=subplot(3,1,3);hold on;
plot(fmri_volumes{3,1},fmri_volumes{3,2},col,'linewidth',2);
plot(ses3_envelope_time-ses3_offset,ses3_envelope,'b');
axis tight;
ylim([-0.2,1.1]);
box on;
title('Session 3');
xlabel('Time [s]')
ylabel('Model and fmri response')

end

