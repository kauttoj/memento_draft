clc;
clear all;
close all;

nworker = 6;
NEED_TO_STANDARDIZE = 0;
 

data_path{1}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/post_processed_normal/bramila/';
data_path{2}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/post_processed_normal/bramila/';
data_path{3}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/post_processed_normal/bramila/';

home = pwd;

S = {'S7','S8', 'S9', 'S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'}; % S13 incomplete

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
      
MAX_IMAGES = [1369,1369,1369];

TR = 1.56;
first_volume_offset = 0.19+TR/2;
EVENT_DURATION = 5;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');
cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/bw_scenes');
TIMING_file = 'memento_bw_scene_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,3);
%TIMING_data = TIMING_data(randperm(length(TIMING_data)));

fprintf('\n-- Processing timing data ---\n')

averaged_volumes = struct();
k=0;

SHIFT = 20 + 20*rand(1,length(TIMING_data));

for i=1:length(TIMING_data)     
    k=k+1;    
    [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data(i).FIRST_time(1),EVENT_DURATION,TR,first_volume_offset);    
    averaged_volumes(k).volumes = fmri_volumes;
    averaged_volumes(k).session = TIMING_data(i).FIRST_ses;
    averaged_volumes(k).id=TIMING_data(i).ID;
    averaged_volumes(k).type = 'BW';
    ind1 = averaged_volumes(k).volumes<=MAX_IMAGES(averaged_volumes(k).session);
        
    k=k+1;
    [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data(i).FIRST_time(1)+SHIFT(i),EVENT_DURATION,TR,first_volume_offset);    
    averaged_volumes(k).volumes = fmri_volumes;
    averaged_volumes(k).session = TIMING_data(i).FIRST_ses;
    averaged_volumes(k).id=TIMING_data(i).ID;
    averaged_volumes(k).type = 'null';
    ind2 = averaged_volumes(k).volumes<=MAX_IMAGES(averaged_volumes(k).session);
    
    m=min(length(ind1),length(ind2));
    
    averaged_volumes(k-1).volumes=averaged_volumes(k-1).volumes(1:m);
    averaged_volumes(k).volumes=averaged_volumes(k).volumes(1:m);    
end

cd(home);

selected_volumes=1:length(averaged_volumes);
averaged_volumes=averaged_volumes(selected_volumes);

if NEED_TO_STANDARDIZE == 1    
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
    BAD_IDs = [];
        
    fprintf('...reading volumes\n ')

    for i = 1:length(averaged_volumes)
        ses = averaged_volumes(i).session;
        cd(data_path{ses});
        
        ind = averaged_volumes(i).volumes;
        
        if max(ind)-min(ind) > 20 || any( diff(ind)~=1)
            max(ind),min(ind)
            error('Too large distance between indices!!!')
        end
        
        MAX_TIMEPOINTS = get_nii_frame([subj,'_mask_detrend_fullreg_filtered_zscored.nii']);
        
        if any(ind>MAX_TIMEPOINTS)
            BAD_FRAMES(end+1)=i;
            BAD_IDs(end+1)=averaged_volumes(i).id;
        else
            nii = load_nii([subj,'_mask_detrend_fullreg_filtered_zscored.nii'],ind);
            img=nii.img;
            avg_img(:,:,:,i)=mean(img,4);
        end
        
    end

    avg_img(:,:,:,BAD_FRAMES)=[];
    
    subj_averaged_volumes{subj_ind} = averaged_volumes;
    subj_averaged_volumes{subj_ind}(BAD_FRAMES)=[];
    
    BAD_FRAMES = [];    
    for i=1:length(subj_averaged_volumes{subj_ind})
       if ismember(subj_averaged_volumes{subj_ind}(i).id,BAD_IDs)>0
           BAD_FRAMES(end+1)=i;
       end
    end    
    subj_averaged_volumes{subj_ind}(BAD_FRAMES)=[];
    avg_img(:,:,:,BAD_FRAMES)=[];    
    
    if size(avg_img,4)~=length(subj_averaged_volumes{subj_ind})
        error('Wrong number of volumes!');
    end
    if nnz(isnan(avg_img))>0
        error('NaN values found!')
    end        
    
    cd(home);
    
    fprintf('...subject has %i averaged volumes. Saving data.\n ',size(avg_img,4))

    save_nii_oma(avg_img,[subj,'_averaged_patterns.nii']); 
        
end

fprintf('\nsaving summary file...\n')

cd(home);
save('averaged_patterns_summary.mat','averaged_volumes','subj_averaged_volumes','S');

fprintf('\n--- ALL DONE!! ----\n')
