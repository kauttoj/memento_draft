clear all;
close all;

TR=1.56;

BW_onset = [...
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

BW_duration = BW_end - BW_onset;

COLOR_onset = BW_end(1:end-1);
COLOR_duration = BW_onset(2:end);

pre_keyframes_onset = [...
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

pre_keyframes_end = [...
    177.8
    409.83
    654.877
    1024.891
    1332.605
    1899.098
    2661.567
    2849.295
    3026.928
    3104.818
    3750.051
    4043.77
    4235.348
    4754.439
    5195.396];
    
pre_keyframes_duration = pre_keyframes_end-pre_keyframes_onset;

keyframes_onset = [...
    580.7   % 2
    920.96  % 3
    1273    % 4
    1532.1  % 5
    1820.4  % 6
    2582.8  % 7
    2990.7  % 8
    3078.12 % 9
    3248.9  % 10
    3411    % 11
    4194.15 % 13 (+0.6 for MEM_5)
    4435.2  % 14
    4695.8  % 15
    4970.05 % 16
    6346.53 % 17
    ];

keyframes_end = [...
        585
        927
        1279
        1537
        1822
        2590
        2995
        3080
        3256
        3422
        4198 % leikkautuu
        4445
        4700
        4975
        6348];
    
keyframes_duration = keyframes_end-keyframes_onset;


final_models{1}.ID = 'BW_scenes';
final_models{1}.onset = BW_onset;
final_models{1}.duration = BW_duration;

final_models{2}.ID = 'COLOR_scenes';
final_models{2}.onset = COLOR_onset;
final_models{2}.duration = COLOR_duration;

final_models{3}.ID = 'key_frames';
final_models{3}.onset = keyframes_onset;
final_models{3}.duration = keyframes_duration;

N_models = 3;

FILE = 'selected_regressors.txt';

formatSpec = '%s%s%s%f%f%s%f%f%s%f%f%s%[^\n\r]';
fileID = fopen(FILE,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t',  'ReturnOnError', false);
fclose(fileID);

N=length(dataArray{1});

ID = dataArray{1};
startTime = dataArray{4};
endTime = dataArray{7};
durationTime = endTime - startTime;

fmri_volumes = [];
model_ind = zeros(1,N);

models(1).ID=ID{1};
model_ind(1)=1;
k=1;

for i=2:N
    if ~strcmp(models(k).ID,ID{i})
        k=k+1;
        models(k).ID=ID{i};
        model_ind(i)=model_ind(i-1)+1;
    else
        model_ind(i)=model_ind(i-1);
    end
end

home = pwd;

ind = unique(model_ind);
bad_models = [];
for i=1:k
    fprintf('model: %s\n',models(i).ID);
    models(i).onset = startTime(find(model_ind==ind(i)));
    models(i).duration = durationTime(find(model_ind==ind(i)));
    models(i).ending = endTime(find(model_ind==ind(i)));
    
    ii = find(models(i).duration<1);
    if ~isempty(ii)
        fprintf('!!!!!!!! deleted too short duration events: %s\n',models(i).ID);
    end
    models(i).onset(ii)=[];
    models(i).duration(ii)=[];
    models(i).ending(ii) = [];
    
    min_dist = 1.5;
    models = join_overlapping(models,min_dist);
    
    if any( models(i).onset(2:end) - models(i).ending(1:end-1) < min_dist)             
       fprintf('!!!!!!!! overlap: %s\n',models(i).ID);
       bad_models=[bad_models,i];
    elseif length(models(i).onset)<12
       fprintf('!!!!!!!! too few events: %s\n',models(i).ID);
       bad_models=[bad_models,i];        
    else
%        cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2014/analysis');
%        [other_models{i},handles] = give_memento_model_volumes(models(i).onset,models(i).duration,0); 
%        close all;
%        cd(home);
    end
end
models(bad_models)=[];
% other_models(bad_models)=[];
% other_models_info = models;
for i=1:length(models)
   final_models{N_models+i}.ID=models(i).ID;
   final_models{N_models+i}.onset=models(i).onset;
   final_models{N_models+i}.duration=models(i).duration;   
end

for i=1:length(final_models)      
   [onset,inds] = split_per_session(final_models{i}.onset);   
   duration{1} = final_models{i}.duration(inds{1});
   duration{2} = final_models{i}.duration(inds{2});
   duration{3} = final_models{i}.duration(inds{3});
   
   %final_models{i} = rmfield(final_models{i},'onset');
   %final_models{i} = rmfield(final_models{i},'duration');
   
   final_models{i}.onset_splitted=onset;
   final_models{i}.duration_splitted=duration;
end

memento_models = final_models;
save('memento_models.mat','memento_models');

%save('other_models_data.mat','other_models','other_models_info');