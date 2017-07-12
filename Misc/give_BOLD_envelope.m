function [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(onset_time,onset_duration,TR,first_volume_offset,THRESHOLD)

if nargin<5
    THRESHOLD=0.4;
end

MAX_VOLUMES = ceil((2*60*60)/TR); % 2h

volume_times = ((1:MAX_VOLUMES)-1)*TR + first_volume_offset;

dt = 0.01;
movie_times = 0:dt:volume_times(end);
movie_envelope = 0*movie_times;

for i=1:length(onset_time)
    ind = find(movie_times>onset_time(i),1,'first'):find(movie_times<onset_time(i)+onset_duration(i),1,'last');
    movie_envelope(ind)=1;        
end

hrf = spm_hrf(dt);
BOLD_envelope = conv(hrf',movie_envelope);
BOLD_envelope=BOLD_envelope(1:length(movie_envelope));
BOLD_envelope = BOLD_envelope/max(BOLD_envelope);
BOLD_envelope_fmri = interp1(movie_times,BOLD_envelope,volume_times);
fmri_volumes = find(BOLD_envelope_fmri>THRESHOLD);

fprintf('\n');
for i=1:length(onset_time)
    fprintf('---- computed BOLD response envelope (onset=%.1fs, duration=%.2fs, dt=%.3fs, offset=%.3fs, th=%.2f) -----\n',onset_time(i),onset_duration(i),dt,first_volume_offset,THRESHOLD);
end

end

