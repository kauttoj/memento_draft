clc;
clear all;
close all;

nworker = 6;
NEED_TO_STANDARDIZE = 1;
% 
% create local cluster
try
    myCluster = gcp('nocreate');
    if isempty(myCluster)
        %delete(gcp)
        myCluster = parcluster('local');
        myCluster.NumWorkers=nworker;
        parpool(myCluster);
    end
    N_workers = myCluster.NumWorkers;
catch err % old matlab?
    if ~matlabpool('size')
       eval(['matlabpool local ',num2str(nworker)]);
    end
    N_workers = matlabpool('size');
end

data_path{1}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session1/post_processed_normal/bramila/';
data_path{2}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session2/post_processed_normal/bramila/';
data_path{3}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session3/post_processed_normal/bramila/';

home = pwd;

S = {'KRON_2','KRON_3','KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13'};

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA');
EPI_mask = compute_analysis_mask(data_path,S);
EPI_mask_ind = find(EPI_mask);
save_nii_oma(EPI_mask,'EPI_mask.nii');
fprintf('EPI mask has %i voxels\n',nnz(EPI_mask));

%nii = load_nii('/triton/becs/scratch/braindata/kauttoj2/code/conn_latest/rois/HO_atlas.img');
%PATTERN_mask = EPI_mask.*(nii.img>0);

nii = load_nii('/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox/Templates/grey.nii');
PATTERN_mask = EPI_mask.*(nii.img>0.25);

maskID=find(PATTERN_mask>0);
save_nii_oma(PATTERN_mask,'PATTERN_mask.nii');
fprintf('PATTERN mask has %i voxels\n',nnz(PATTERN_mask));
      
MAX_IMAGES = [1319,1319,1275];

TR = 1.56;
first_volume_offset = 0.19+TR/2;
EVENT_DURATION = 5;

cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/chrono_keyframes');
TIMING_file = 'memento_chrono_keyframes_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,3);
%TIMING_data = TIMING_data(randperm(length(TIMING_data)));

% Ses1: 24min, 26min
% Ses2: 2min, 
k=0;

k=k+1;TIMING_data_null(k).FIRST_ses=3;
TIMING_data_null(k).FIRST_time = 60*2;
k=k+1;TIMING_data_null(k).FIRST_ses=3;
TIMING_data_null(k).FIRST_time = 60*6;
k=k+1;TIMING_data_null(k).FIRST_ses=3;
TIMING_data_null(k).FIRST_time = 60*10;
k=k+1;TIMING_data_null(k).FIRST_ses=3;
TIMING_data_null(k).FIRST_time = 60*15;
k=k+1;TIMING_data_null(k).FIRST_ses=3;
TIMING_data_null(k).FIRST_time = 60*20;
k=k+1;TIMING_data_null(k).FIRST_ses=3;
TIMING_data_null(k).FIRST_time = 60*24;

k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*3;
k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*7;
k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*11;
k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*17;
k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*19;
k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*22;
k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*25;
k=k+1;TIMING_data_null(k).FIRST_ses=2;
TIMING_data_null(k).FIRST_time = 60*29;

k=k+1;TIMING_data_null(k).FIRST_ses=1;
TIMING_data_null(k).FIRST_time = 60*29;

for i=1:length(TIMING_data_null)
    TIMING_data_null(i).ID=i;
    TIMING_data_null(i).FIRST_time_str = sec2min(TIMING_data_null(i).FIRST_time);
end

%eka 23-loppuun

averaged_volumes = struct();
k=0;
for i=1:length(TIMING_data) 
    k=k+1;    
    [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data(i).FIRST_time(1),EVENT_DURATION,TR,first_volume_offset);    
    averaged_volumes(k).volumes = fmri_volumes;
    averaged_volumes(k).session = TIMING_data(i).FIRST_ses;
    averaged_volumes(k).id=TIMING_data(i).ID;
    averaged_volumes(k).type = 'keyframe';
    ind1 = averaged_volumes(k).volumes<=MAX_IMAGES(averaged_volumes(k).session);
    m=length(ind1);
    averaged_volumes(k).volumes=averaged_volumes(k).volumes(1:m);    
end
    
for i=1:length(TIMING_data_null)
    k=k+1;
    [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data_null(i).FIRST_time(1),EVENT_DURATION,TR,first_volume_offset);
    averaged_volumes(k).volumes = fmri_volumes;
    averaged_volumes(k).session = TIMING_data(i).FIRST_ses;
    averaged_volumes(k).id=TIMING_data_null(i).ID;  
    averaged_volumes(k).type = 'nullframe';
    ind2 = averaged_volumes(k).volumes<=MAX_IMAGES(averaged_volumes(k).session);   
    m=length(ind2);        
    averaged_volumes(k).volumes=averaged_volumes(k).volumes(1:m);    
end

cd(home);

selected_volumes=1:length(averaged_volumes);
averaged_volumes=averaged_volumes(selected_volumes);

if NEED_TO_STANDARDIZE == 1
    fprintf('\n-- Running standardization step ---\n')
    for ses=1:3
        fprintf('Session %i\n',ses)
        cd(data_path{ses});
        parfor subj_ind = 1:length(S);
            subj = S{subj_ind};
            fprintf('...subject %s\n',subj);
            nii = load_nii([subj,'_mask_detrend_fullreg_filtered.nii']);
            data=nii.img;
            siz=size(data);
            data=reshape(data,[],siz(4))';
            data = zscore(data);
            data=data';
            data=reshape(data,siz);
            nii.img = data;
            save_nii(nii,[subj,'_mask_detrend_fullreg_filtered_zscored.nii']);
        end
    end
end
cd(home);

fprintf('\n-- Starting data extraction ---\n')

for subj_ind = 1:length(S);
    
    subj = S{subj_ind};
    fprintf('Subject %s\n',subj);
    
    avg_img = nan(91,109,91,length(averaged_volumes));
    
    BAD_FRAMES = [];
        
    fprintf('...reading volumes ')

    for i = 1:length(averaged_volumes)
        ses = averaged_volumes(i).session;
        cd(data_path{ses});
        
        ind = averaged_volumes(i).volumes;
        
        if max(ind)-min(ind) > 30 || any( diff(ind)~=1)
            max(ind),min(ind)
            error('Too large distance between indices!!!')
        end
        
        MAX_TIMEPOINTS = get_nii_frame([subj,'_mask_detrend_fullreg_filtered_zscored.nii']);
        
        if any(ind>MAX_TIMEPOINTS)
            BAD_FRAMES(end+1)=i;
        else
            nii = load_nii([subj,'_mask_detrend_fullreg_filtered_zscored.nii'],ind);
            img=nii.img;
            avg_img(:,:,:,i)=mean(img,4);
        end
        
    end

    avg_img(:,:,:,BAD_FRAMES)=[];
    
    subj_averaged_volumes{subj_ind} = averaged_volumes;
    subj_averaged_volumes{subj_ind}(BAD_FRAMES)=[];
    
    if size(avg_img,4)~=length(subj_averaged_volumes{subj_ind})
        error('Wrong number of volumes!');
    end
    if nnz(isnan(avg_img))>0
        error('NaN values found!')
    end
    
    fprintf('\n done, saving file...\n')
    
    cd(home);

    save_nii_oma(avg_img,[subj,'_averaged_patterns.nii']); 
        
end

fprintf('saving summary file...\n')

cd(home);
save('averaged_patterns_summary.mat','averaged_volumes','subj_averaged_volumes','S');

fprintf('\n--- ALL DONE!! ----\n')
