function prepare_data_ver3_looper(FOLDER,DELAY,DURATION)


NEED_TO_STANDARDIZE = 0;

%% create local cluster
% nworker = 6;
% try
%     myCluster = gcp('nocreate');
%     if isempty(myCluster)
%         %delete(gcp)
%         myCluster = parcluster('local');
%         myCluster.NumWorkers=nworker;
%         parpool(myCluster);
%     end
%     N_workers = myCluster.NumWorkers;
% catch err % old matlab?
%     if ~matlabpool('size')
%        eval(['matlabpool local ',num2str(nworker)]);
%     end
%     N_workers = matlabpool('size');
% end

data_path{1}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/post_processed_normal/bramila/';
data_path{2}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/post_processed_normal/bramila/';
data_path{3}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/post_processed_normal/bramila/';

cd(FOLDER);

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
first_volume_offset = 0.19 + TR/2;
EVENT_DURATION = DURATION;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');
cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/bw_scenes');
TIMING_file = 'memento_bw_scene_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,3);

fprintf('\n-- Processing timing data ---\n')

SHIFT = [42.0,26.2,42.6,25.4,41.0,40.8,38.1,40.6,40.1,43.2,38.5,35.0,29.0,41.2,33.8,31.8,40.6,27.8,25.1,29.1,43.0,43.3];

averaged_volumes = struct();
k=0;
for i=1:length(TIMING_data)     
    k=k+1;    
    [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data(i).FIRST_time(1)+DELAY,EVENT_DURATION,TR,first_volume_offset);    
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

cd(FOLDER);

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
cd(FOLDER);

fprintf('\n-- Starting data extraction ---\n')

subj_averaged_volumes=cell(1,length(S));
for subj_ind = 1:length(S);
    
    subj = S{subj_ind};
    fprintf('Subject %s\n',subj);
    
    avg_img = nan(91,109,91,length(averaged_volumes));
    
    BAD_FRAMES = [];
    BAD_IDs = [];
        
    fprintf('...reading volumes ')

    for i = 1:length(averaged_volumes)
        ses = averaged_volumes(i).session;
        cd(data_path{ses});
        
        ind = averaged_volumes(i).volumes;
        
        if max(ind)-min(ind) > 25 || any( diff(ind)~=1)
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
    
    cd(FOLDER);
    
    fprintf('...subject has %i averaged volumes\n',size(avg_img,4));

    save_nii_oma(avg_img,[subj,'_averaged_patterns.nii']); 
        
end

fprintf('saving summary file...\n')

cd(FOLDER);

save('averaged_patterns_summary.mat','averaged_volumes','subj_averaged_volumes','S');

fprintf('\n--- ALL DONE!! ----\n')
