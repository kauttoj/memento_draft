function prepare_data_ver6_looper(FOLDER,DELAY,HARD_LIMIT)

nworker = 7;
NEED_TO_STANDARDIZE = 0;
 
data_path{1}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/post_processed_normal/bramila/';
data_path{2}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/post_processed_normal/bramila/';
data_path{3}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/post_processed_normal/bramila/';

home = FOLDER;

cd(FOLDER);

S = {'S7','S8', 'S9', 'S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'}; % S13 incomplete

cd ..
cd ..
nii=load_nii('ORIG_KRON_mask.nii');
EPI_mask=nii.img;
fprintf('EPI mask has %i voxels\n',nnz(EPI_mask));

cd(FOLDER);

nii = load_nii('/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox/Templates/grey.nii');
PATTERN_mask = EPI_mask.*(nii.img>0.15);

maskID=find(PATTERN_mask>0);
save_nii_oma(PATTERN_mask,'PATTERN_mask.nii');
fprintf('PATTERN mask has %i voxels\n',nnz(PATTERN_mask));
      
MAX_IMAGES = [1369,1369,1369];

TR = 1.56;
first_volume_offset = 0.18+TR/2;
EVENT_DURATION = (HARD_LIMIT-1)*TR;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');
cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/loop_keyframes');
TIMING_file = 'memento_loop_keyframe_timings_FINAL.txt';
TIMING_data = parse_timing_file(TIMING_file,2);

cd(FOLDER);

fprintf('\n-- Processing timing data ---\n')

averaged_volumes = struct();
k=0;
for i=1:length(TIMING_data)     
    
%     k=k+1;    
%     [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data(i).FIRST_time(1)+DELAY,EVENT_DURATION,TR,first_volume_offset,0.49);    
%     %averaged_volumes(k).volumes = fmri_volumes;
%     averaged_volumes(k).session = TIMING_data(i).FIRST_ses;
%     averaged_volumes(k).id=TIMING_data(i).ID;
%     averaged_volumes(k).type = 'initial';
%    
%     [BOLD_envelope_fmri,ind]=sort(BOLD_envelope_fmri,'descend');
%     averaged_volumes(k).volumes = sort(ind(1:HARD_LIMIT));
        
    k=k+1;
    [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(TIMING_data(i).SECOND_time(1)+DELAY,EVENT_DURATION,TR,first_volume_offset,0.49);
    %averaged_volumes(k).volumes = fmri_volumes;
    averaged_volumes(k).session = TIMING_data(i).SECOND_ses;
    averaged_volumes(k).id=TIMING_data(i).ID;
    averaged_volumes(k).type = 'keyframe';

    [BOLD_envelope_fmri,ind]=sort(BOLD_envelope_fmri,'descend');
    averaged_volumes(k).volumes = sort(ind(1:HARD_LIMIT));
    
end

% nulltimes{1} = [110,310,764,1130,13,2003] + [2.7472    2.3766    2.8785    1.9672    0.1071    2.5474];
% nulltimes{2} = [99,252,605,1410,1793] + [2.8020    2.0362    2.2732    2.2294    1.1767];
% nulltimes{3} = [47,315,525,1097,1551,1829] + [1.9664    0.5136    2.1181    0.0955    0.8308    0.1385];
% 
% for s=1:3,
%     for i=1:length(nulltimes{s})
%         k=k+1;
%         [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(nulltimes{s}(i),EVENT_DURATION,TR,first_volume_offset,0.49);
%         %averaged_volumes(k).volumes = fmri_volumes;
%         averaged_volumes(k).session = s;
%         averaged_volumes(k).id=i;
%         averaged_volumes(k).type = 'null';
%         
%         [BOLD_envelope_fmri,ind]=sort(BOLD_envelope_fmri,'descend');
%         averaged_volumes(k).volumes = sort(ind(1:HARD_LIMIT));
%     end
% end

cd(home);

selected_volumes=1:length(averaged_volumes);
averaged_volumes=averaged_volumes(selected_volumes);

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
