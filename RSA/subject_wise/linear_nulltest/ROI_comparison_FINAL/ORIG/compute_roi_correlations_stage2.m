function compute_roi_correlations_stage2(DELAYS,DURATION)

diary('ROI_correlation_diary.txt')

S = {'S7','S8', 'S9','S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'}; % S13 incomplete

MAX_IMAGES = [1369,1369,1369];

HOME = pwd;

TR = 1.56;
first_volume_offset = 0.19 + TR/2;
EVENT_DURATION = DURATION;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');
cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/loop_keyframes');
TIMING_file = 'memento_loop_keyframe_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,2);

fprintf('\n-- Processing timing data ---\n')

for l=1:length(DELAYS)
    DELAY=DELAYS(l);
    averaged_volumes = struct();
    k=0;
    for i=1:length(TIMING_data)
        
        k=k+1;
        [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data(i).SECOND_time(1)+DELAY,EVENT_DURATION,TR,first_volume_offset,0.45);
        averaged_volumes(k).volumes = fmri_volumes;
        averaged_volumes(k).session = TIMING_data(i).SECOND_ses;
        averaged_volumes(k).id=TIMING_data(i).ID;
        averaged_volumes(k).type = 'keyframe';
        ind = averaged_volumes(k).volumes<=MAX_IMAGES(averaged_volumes(k).session);
        averaged_volumes(k).volumes=averaged_volumes(k).volumes(ind);
               
    end
    averaged_volumes_all{l}=averaged_volumes;
    
end

cd(HOME);
load('roi_and_mask_data.mat');

save('volume_selection_data.mat','averaged_volumes_all','EVENT_DURATION','DELAYS','S','TIMING_data','-v7.3');

N_rois = length(my_masks);
for i=1:N_rois
    mask_voxels(i)=length(my_masks{i});
end

fprintf('\n--- Number of ROIs is %i ---\n',N_rois);

fprintf('\n--- Starting data extraction ---\n');

for subj_ind = 1:length(S);
    
    subj = S{subj_ind};
    fprintf('\nSubject %s\n',subj);
               
    clear data;
    load([subj,'_roi_data.mat']);
    
    if length(data{1})~=N_rois
        error('ROI count does not match');
    end
    
    subj_averaged_volumes=cell(length(averaged_volumes_all),1);
        
    for l=1:length(averaged_volumes_all)
        
        fprintf('... delay %.2fs (%i of %i):',DELAYS(l),l,length(averaged_volumes_all));
        
        averaged_volumes = averaged_volumes_all{l};
        
        avg_img=cell(1,N_rois);
        for i=1:N_rois
            avg_img{i}=nan(length(averaged_volumes),mask_voxels(i));
        end
        
        BAD_FRAMES = [];
        
        for i = 1:length(averaged_volumes)
            ses = averaged_volumes(i).session;
            ind = averaged_volumes(i).volumes;
            
            if isempty(ind)
                BAD_FRAMES(end+1)=i;
                break;
            end
            
            if max(ind)-min(ind) > 20 || any( diff(ind)~=1)
                max(ind),min(ind)
                error('Too large distance between indices!!!')
            end
            
            MAX_TIMEPOINTS = size(data{ses}{1},1);
            
            if any(ind>MAX_TIMEPOINTS)
                BAD_FRAMES(end+1)=i;
            else
                for j=1:N_rois
                    avg_img{j}(i,:)=mean(data{ses}{j}(ind,:),1);
                end
            end
        end
        
        for j=1:N_rois
            avg_img{j}(BAD_FRAMES,:)=[];
            if nnz(isnan(avg_img{j}))>0
                error('NaN values found!')
            end
        end
        
        averaged_volumes(BAD_FRAMES)=[];
        subj_averaged_volumes{l} = averaged_volumes;
        
        if size(avg_img{1},1)~=length(subj_averaged_volumes{l})
            error('Wrong number of volumes!');
        end                
        
        fprintf(' subject has %i averaged volumes\n',size(avg_img{1},1));
        
        avg_img_all{l}=avg_img;
        
    end
    
    save([subj,'_averaged_patterns.mat'],'avg_img_all','subj_averaged_volumes','subj','-v7.3');  
    
    clear data;
    fprintf('Computing correlations\n');
    
    all_correlations = nan(length(avg_img_all),N_rois);
    all_cormats = cell(length(avg_img_all),N_rois);
    for i=1:length(avg_img_all)
        dim = size(avg_img_all{i}{1},1);
        inds = triu(ones(dim,dim),1)>0;
        for j=1:N_rois
            data=avg_img_all{i}{j};
            cormat = corr(data');
            all_cormats{i,j}=cormat;
            if dim~=size(cormat,1)
                error('Dimension mismatch!');
            end
            all_correlations(i,j)=mean(cormat(inds));
        end        
        if nnz(isnan(all_correlations(i,:)))>0
            error('NaN correlations found!');
        end
        fprintf('...correlations for delay %.2fs: mean=%.3f, low=%.3f, high=%.3f\n',DELAYS(i),mean(all_correlations(i,:)),prctile(all_correlations(i,:),5),prctile(all_correlations(i,:),95));
    end
    save([subj,'_correlation_results.mat'],'all_correlations','all_cormats','-v7.3');
end

diary off

fprintf('\n--- ALL DONE!! ----\n')
