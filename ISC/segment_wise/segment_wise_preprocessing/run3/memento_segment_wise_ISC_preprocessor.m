clc;
clear all;
close all;

addpath('/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila');
addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');

HOMEPATH = pwd;
CURRPATH = HOMEPATH;

%% SEGMENT TIMING DATA

TIMING_FILE = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/memento_final_timings.txt';
MOVIE_SEGMENTS = parse_timing_file(TIMING_FILE);
%load('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/repeating_segment_data_stage1.mat','full_corrdata');
%SEGMENT_TIMINGS = load('???????');

START_SHIFT = 20; % seconds
END_SHIFT = 3; % seconds

SELECTED_IDs=1:22;
SELECTED_IDs(SELECTED_IDs==16)=[];

ind=[];
k=0;
for i=1:length(MOVIE_SEGMENTS)
    MOVIE_SEGMENTS(i).ORIGINAL_time = MOVIE_SEGMENTS(i).ORIGINAL_time + [START_SHIFT,END_SHIFT];
    MOVIE_SEGMENTS(i).KRON_time=MOVIE_SEGMENTS(i).KRON_time + [START_SHIFT,END_SHIFT];
    if ismember(MOVIE_SEGMENTS(i).ID,SELECTED_IDs)        
        if diff(MOVIE_SEGMENTS(i).ORIGINAL_time)>130                       
            ind(end+1)=i;
            if abs(diff(MOVIE_SEGMENTS(i).ORIGINAL_time)-diff(MOVIE_SEGMENTS(i).KRON_time))>1e-10
                error('Bad timing!')
            end
        end
    end
end
MOVIE_SEGMENTS = MOVIE_SEGMENTS(ind);

%% INPUT DATA PATHS
%Memento3part_run1.avi
%Memento3part_run2.avi
%Memento3part_run3.avi

ORIGINAL_SUBJ = {'S7','S9','S10','S12','S13','S15','S16','S17','S21','S22','S23'};

ORIGINAL_EPI_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/FunImgNormNii/';
ORIGINAL_EPI_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/FunImgNormNii/';
ORIGINAL_EPI_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/FunImgNormNii/';
ORIGINAL_MOTION_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/RealignParameter/';
ORIGINAL_MOTION_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/RealignParameter/';
ORIGINAL_MOTION_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/RealignParameter/';
ORIGINAL_MASK_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/post_processed_normal/bramila/';
ORIGINAL_MASK_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/post_processed_normal/bramila/';
ORIGINAL_MASK_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/post_processed_normal/bramila/';
ORIGINAL_FILENAME = '/4D.nii';
ORIGINAL_MASK_FILENAME = '_analysis_mask.nii';
ORIGINAL_MOTION_FILENAME = '/4D_motion_params.txt';

KRON_SUBJ = {'KRON_2','KRON_3','KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13'};

KRON_EPI_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session1/FunImgNormNii/';
KRON_EPI_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session2/FunImgNormNii/';
KRON_EPI_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session3/FunImgNormNii/';
KRON_MOTION_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session1/RealignParameter/';
KRON_MOTION_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session2/RealignParameter/';
KRON_MOTION_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session3/RealignParameter/';
KRON_MASK_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session1/post_processed_normal/bramila/';
KRON_MASK_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session2/post_processed_normal/bramila/';
KRON_MASK_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session3/post_processed_normal/bramila/';
KRON_FILENAME = '/4D.nii';
KRON_MASK_FILENAME = '_analysis_mask.nii';
KRON_MOTION_FILENAME = '/4D_motion_params.txt';


%% Parameters

SELECTED_SEGMENTS = 1:length(MOVIE_SEGMENTS);
INTERPOLATION_METHOD = 'linear';
TR = 1.56;
DIARY_NAME = 'memento_segment_wise_ISC_analysis_diary_preprocessor.txt';
N_WORKER = 6;
TR_isc = (TR/3);
DELETE_SEGMENT_NII = 0;

CFG.HP_limit = 0.01;
CFG.FWHM = 6;
CFG.TR = TR;
CFG.DETREND_type = 'polynomial-nodemean';

delete(DIARY_NAME);
diary(DIARY_NAME)

%% Create group mask

group_mask = true(91,109,91);
for SES=1:3
    k=0;
    for SUB = ORIGINAL_SUBJ
        k=k+1;
        nam = [ORIGINAL_MASK_PATH{SES},SUB{1},ORIGINAL_MASK_FILENAME];
        nii=load_nii(nam);
        mask = nii.img>0;
        fprintf('mask size %i (%s)\n',nnz(mask),nam);
        group_mask = group_mask.*mask;
    end
    k=0;
    for SUB = KRON_SUBJ
        k=k+1;
        nam = [KRON_MASK_PATH{SES},SUB{1},KRON_MASK_FILENAME];
        nii=load_nii(nam);
        mask = nii.img>0;
        fprintf('mask size %i (%s)\n',nnz(mask),nam);        
        group_mask = group_mask.*mask;
    end
    fprintf('\n');
end
group_mask=group_mask>0;
fprintf('>>>>>> GROUP MASK SIZE %i \n',nnz(group_mask));

nii=load_nii('/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila/external/grey.nii');
grey_mask = nii.img>0.25;
group_mask=group_mask.*grey_mask;

save_nii_oma(group_mask,'memento_segment_wise_ISC_mask_preprocessor.nii');
group_mask_ind = find(group_mask);

diary(DIARY_NAME)

%% Check number of volumes
original_volumes=[];
for SES=1:3
    k=0;
    for SUB = ORIGINAL_SUBJ
        k=k+1;
        original_volumes(SES,k)=get_nii_frame([ORIGINAL_EPI_PATH{SES},SUB{1},ORIGINAL_FILENAME]);
    end
end
ORIGINAL_max_volumes = min(original_volumes');
fprintf('\n>>>>>> Max available volumes (ORIGINAL): %i %i %i\n',ORIGINAL_max_volumes(1),ORIGINAL_max_volumes(2),ORIGINAL_max_volumes(3));

kron_volumes=[];
for SES=1:3
    k=0;
    for SUB = KRON_SUBJ
        k=k+1;
        kron_volumes(SES,k)=get_nii_frame([KRON_EPI_PATH{SES},SUB{1},KRON_FILENAME]);
    end
end
KRON_max_volumes = min(kron_volumes');
fprintf('\n>>>>>> Max available volumes (KRON): %i %i %i\n\n',KRON_max_volumes(1),KRON_max_volumes(2),KRON_max_volumes(3));

diary(DIARY_NAME)

%% Test non-neural group differences (full-sessions)
clear cfg;
cfg.prepro_suite = 'spm';

for SES=1:3
    original_full_fwd{SES} = nan(ORIGINAL_max_volumes(SES),length(ORIGINAL_SUBJ));
    k=0;
    for SUB = ORIGINAL_SUBJ
        k=k+1;
        cfg.motionparam = [ORIGINAL_MOTION_PATH{SES},SUB{1},ORIGINAL_MOTION_FILENAME];                        
        [fwd,rms]=bramila_framewiseDisplacement(cfg);
        ORIGINAL_full_rms(SES,k)=mean(rms);
        ORIGINAL_full_fwd{SES}(:,k)=fwd(1:ORIGINAL_max_volumes(SES));
        
        A = load(cfg.motionparam);
        ORIGINAL_motiondata{SES,k} = A;
    end
end

for SES=1:3    
    kron_full_fwd{SES} = nan(KRON_max_volumes(SES),length(KRON_SUBJ));
    k=0;
    for SUB = KRON_SUBJ
        k=k+1;
        cfg.motionparam = [KRON_MOTION_PATH{SES},SUB{1},KRON_MOTION_FILENAME]; 
        [fwd,rms]=bramila_framewiseDisplacement(cfg);
        KRON_full_rms(SES,k)=mean(rms);
        KRON_full_fwd{SES}(:,k)=fwd(1:KRON_max_volumes(SES));
        
        A = load(cfg.motionparam);
        KRON_motiondata{SES,k} = A;   
    end
end

ORIGINAL_mat_ind = find(triu(ones(length(ORIGINAL_SUBJ),length(ORIGINAL_SUBJ)),1));
KRON_mat_ind = find(triu(ones(length(KRON_SUBJ),length(KRON_SUBJ)),1));

fprintf('\n\nTesting group differences of RMS and FWD ISC\n')
for SES=1:3
    testdata = [ORIGINAL_full_rms(SES,:),KRON_full_rms(SES,:)];
    rms_full_difference_stats{SES} = bramila_two_sample_test(testdata,...
        [ones(1,size(ORIGINAL_full_rms,2)),2*ones(1,size(KRON_full_rms,2))],...
        10000,N_WORKER,0,0);
    fprintf('\n>>>>>> Session %i motion RMS: pval = %.3f\n\n',SES,2*min(rms_full_difference_stats{SES}.pvals)); 
    
    m = min(size(ORIGINAL_full_fwd{SES},1),size(KRON_full_fwd{SES},1));
    
    isc1 = corr(ORIGINAL_full_fwd{SES}(1:m,:));
    isc1 = isc1(ORIGINAL_mat_ind);
    isc2 = corr(KRON_full_fwd{SES}(1:m,:));
    isc2 = isc2(KRON_mat_ind);
    
    testdata = [isc1',isc2'];
    testdata = atanh(testdata); % convert to the full real axis
    
    fwd_full_ISC_difference_stats{SES} = bramila_two_sample_test(testdata,...
        [ones(1,length(isc1)),2*ones(1,length(isc2))],...
        10000,N_WORKER,0,0);
    
    fprintf('\n>>>>>> Session %i motion ISC (%i vols): pval = %.3f\n\n',SES,m,2*min(fwd_full_ISC_difference_stats{SES}.pvals));    
end

clear isc1 isc2 testdata

%% Start analysis pipeline
kk=0;
fprintf('----Volume selection data (NOTE: requested times and volumes are shifted by HRF!)------\n');
for segment_nr = SELECTED_SEGMENTS   
    kk=kk+1;
    
    % ORIGINAL
    ORIGINAL_SES = MOVIE_SEGMENTS(segment_nr).ORIGINAL_ses;    
    ORIGINAL_movie_times = MOVIE_SEGMENTS(segment_nr).ORIGINAL_time;    
    ORIGINAL_fmri_times = ((1:ORIGINAL_max_volumes(ORIGINAL_SES)) -1)*TR;   
    vol_min = find(ORIGINAL_movie_times(1)>ORIGINAL_fmri_times,1,'last');
    vol_max = find(ORIGINAL_movie_times(2)<ORIGINAL_fmri_times,1,'first');       
    ORIGINAL_vol = (vol_min:vol_max);         
    if ORIGINAL_movie_times(1)-ORIGINAL_fmri_times(ORIGINAL_vol(1)) > -ORIGINAL_movie_times(2)+ORIGINAL_fmri_times(ORIGINAL_vol(end))
        ORIGINAL_additional_vol = 0;
    else
        ORIGINAL_additional_vol = 1;
    end
    
    % KRON
    KRON_SES = MOVIE_SEGMENTS(segment_nr).KRON_ses;    
    KRON_movie_times = MOVIE_SEGMENTS(segment_nr).KRON_time;
    KRON_fmri_times = ((1:KRON_max_volumes(KRON_SES)) -1)*TR;
    vol_min = find(KRON_movie_times(1)>KRON_fmri_times,1,'last');
    vol_max = find(KRON_movie_times(2)<KRON_fmri_times,1,'first');
    KRON_vol = (vol_min:vol_max);             
    if KRON_movie_times(1)-KRON_fmri_times(KRON_vol(1)) > -KRON_movie_times(2)+KRON_fmri_times(KRON_vol(end))
        KRON_additional_vol = 0;
    else
        KRON_additional_vol = 1;
    end
      
    if length(ORIGINAL_vol)>length(KRON_vol)
        if KRON_additional_vol==0
            KRON_vol=[KRON_vol(1)-1,KRON_vol];
        else
            KRON_vol=[KRON_vol,KRON_vol(end)+1];
        end
    elseif length(ORIGINAL_vol)<length(KRON_vol)
        if ORIGINAL_additional_vol==0
            ORIGINAL_vol=[ORIGINAL_vol(1)-1,ORIGINAL_vol];
        else
            ORIGINAL_vol=[ORIGINAL_vol,ORIGINAL_vol(end)+1];
        end        
    end
    
    clear fmri_times requested_times movie_times SES vols;
    
    ORIGINAL_fmri_times = ORIGINAL_fmri_times(ORIGINAL_vol);
    ORIGINAL_requested_times = ORIGINAL_movie_times(1):TR_isc:ORIGINAL_movie_times(2);
    ORIGINAL_requested_times(ORIGINAL_requested_times>=ORIGINAL_fmri_times(end))=[];    
    N_vols1=length(ORIGINAL_requested_times);
    if ~(ORIGINAL_fmri_times(1)<ORIGINAL_requested_times(1) && ORIGINAL_fmri_times(end)>ORIGINAL_requested_times(end))
        error('times are incorrect!')
    end
    ORIGINAL_requested_time{kk}=ORIGINAL_requested_times;
    ORIGINAL_vols{kk} = ORIGINAL_vol;  
    ORIGINAL_fmri_time{kk}=ORIGINAL_fmri_times;
    ORIGINAL_ses(kk)=ORIGINAL_SES;        
    
    fprintf('ORIGINAL film, ID %i (scan session %i):\n    sequence time %.2f-%.2f (requested %.2f-%.2f)\n    fmri times %.2f-%.2f (vols %i-%i)\n',...
        MOVIE_SEGMENTS(segment_nr).ID,...
        ORIGINAL_ses(kk),...
        ORIGINAL_movie_times(1)-START_SHIFT,...
        ORIGINAL_movie_times(end)-END_SHIFT,...  
        ORIGINAL_requested_time{kk}(1),...
        ORIGINAL_requested_time{kk}(end),...
        ORIGINAL_fmri_time{kk}(1),...
        ORIGINAL_fmri_time{kk}(end),...        
        ORIGINAL_vols{kk}(1),...
        ORIGINAL_vols{kk}(end));
       
    KRON_fmri_times = KRON_fmri_times(KRON_vol);
    KRON_requested_times = KRON_movie_times(1):TR_isc:KRON_movie_times(2);
    KRON_requested_times(KRON_requested_times>=KRON_fmri_times(end))=[];    
    N_vols2=length(KRON_requested_times);
    if ~(KRON_fmri_times(1)<KRON_requested_times(1) && KRON_fmri_times(end)>KRON_requested_times(end))
        error('times are incorrect!')
    end
    KRON_requested_time{kk}=KRON_requested_times;
    KRON_vols{kk} = KRON_vol;  
    KRON_fmri_time{kk}=KRON_fmri_times;
    KRON_ses(kk)=KRON_SES;
    fprintf('KRON film, ID %i (scan session %i):\n    sequence time %.2f-%.2f (requested %.2f-%.2f)\n    fmri times %.2f-%.2f (vols %i-%i)\n',...
        MOVIE_SEGMENTS(segment_nr).ID,...
        KRON_ses(kk),...
        KRON_movie_times(1)-START_SHIFT,...
        KRON_movie_times(end)-END_SHIFT,...  
        KRON_requested_time{kk}(1),...
        KRON_requested_time{kk}(end),...
        KRON_fmri_time{kk}(1),...
        KRON_fmri_time{kk}(end),...        
        KRON_vols{kk}(1),...
        KRON_vols{kk}(end));    
    if N_vols1~=N_vols2
       error('Volume count is different!')
    end                    
    
    
end

kk=0;
for segment_nr = SELECTED_SEGMENTS   
    kk=kk+1;
    
    fprintf('\n\n------------ Starting ISC calculations for segment %i/%i (ID %i) ---------------\n',kk,length(SELECTED_SEGMENTS),MOVIE_SEGMENTS(segment_nr).ID);
        
    CURRPATH = [HOMEPATH,filesep,sprintf('Seg_ID_%i',MOVIE_SEGMENTS(segment_nr).ID)];
    mkdir(CURRPATH);
                      
    funAlign_createpool(N_WORKER);
    
    % RUN PREPROCESSOR
    tempdata = cell(1,length(KRON_SUBJ)); 
    parfor k=1:length(KRON_SUBJ)
        SUB = KRON_SUBJ{k};                 
        nam = [KRON_EPI_PATH{KRON_ses(kk)},SUB,KRON_FILENAME];
        fprintf('...subj %s: %s\n',SUB,nam);
        nii=load_nii(nam,KRON_vols{kk});
        save_nii(nii,[CURRPATH,filesep,SUB,'.nii']);  
        save_text([CURRPATH,filesep,SUB,'_motion.txt'],KRON_motiondata{KRON_ses(kk),k}(KRON_vols{kk},:));
        data=run_segment_preprocessor([CURRPATH,filesep,SUB,'.nii'],[CURRPATH,filesep,SUB,'_motion.txt'],CURRPATH,CFG);
        delete([CURRPATH,filesep,SUB,'_motion.txt']);
        if DELETE_SEGMENT_NII==1
            delete([CURRPATH,filesep,SUB,'.nii']);
        end
        data=reshape(data,[],size(data,4))';
        data=data(:,group_mask_ind);
        if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
            error('data seems bad!')
        end                
        data = interp1(KRON_fmri_time{kk},data,KRON_requested_time{kk},INTERPOLATION_METHOD);
        tempdata{k}=data;
    end      
    KRON_data = zeros(length(KRON_requested_time{kk}),length(group_mask_ind),length(KRON_SUBJ),'single');
    for k=1:length(KRON_SUBJ)
        KRON_data(:,:,k)=tempdata{k};
    end   
    clear tempdata data nii;
    
    diary(DIARY_NAME);
    
    tempdata = cell(1,length(ORIGINAL_SUBJ)); 
    parfor k=1:length(ORIGINAL_SUBJ)
        SUB = ORIGINAL_SUBJ{k};                 
        nam = [ORIGINAL_EPI_PATH{ORIGINAL_ses(kk)},SUB,ORIGINAL_FILENAME];
        fprintf('...subj %s: %s\n',SUB,nam);
        nii=load_nii(nam,ORIGINAL_vols{kk});
        save_nii(nii,[CURRPATH,filesep,SUB,'.nii']);             
        save_text([CURRPATH,filesep,SUB,'_motion.txt'],ORIGINAL_motiondata{ORIGINAL_ses(kk),k}(ORIGINAL_vols{kk},:));
        data=run_segment_preprocessor([CURRPATH,filesep,SUB,'.nii'],[CURRPATH,filesep,SUB,'_motion.txt'],CURRPATH,CFG);
        delete([CURRPATH,filesep,SUB,'_motion.txt']);
        if DELETE_SEGMENT_NII==1
            delete([CURRPATH,filesep,SUB,'.nii']);
        end
        data=reshape(data,[],size(data,4))';
        data=data(:,group_mask_ind);
        if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
            error('data seems bad!')
        end                
        data = interp1(ORIGINAL_fmri_time{kk},data,ORIGINAL_requested_time{kk},INTERPOLATION_METHOD);
        tempdata{k}=data;
    end      
    ORIGINAL_data = zeros(length(ORIGINAL_requested_time{kk}),length(group_mask_ind),length(ORIGINAL_SUBJ),'single');
    for k=1:length(ORIGINAL_SUBJ)
        ORIGINAL_data(:,:,k)=tempdata{k};
    end   
    clear tempdata data nii;
        
    %% Compute fwd difference stats
    
    isc1 = corr(ORIGINAL_full_fwd{ORIGINAL_ses(kk)}(ORIGINAL_vols{kk},:));
    isc1 = isc1(ORIGINAL_mat_ind);
    isc2 = corr(KRON_full_fwd{KRON_ses(kk)}(KRON_vols{kk},:));
    isc2 = isc2(KRON_mat_ind);
    
    testdata = [isc1',isc2'];
    testdata = atanh(testdata); % convert to the full real axis
    
    fprintf('\n');
    fwd_ISC_difference_stats{kk} = bramila_two_sample_test(testdata,...
        [ones(1,length(isc1)),2*ones(1,length(isc2))],...
        10000,N_WORKER,0,0);
    
    fprintf('\n>>>>>> Motion fwd ISC test (%i vols, %i vs. %i): pval = %.3f\n',length(ORIGINAL_vols{kk}),length(isc1),length(isc2),2*min(fwd_ISC_difference_stats{kk}.pvals));       
    
    rms1=sqrt(mean((ORIGINAL_full_fwd{ORIGINAL_ses(kk)}(ORIGINAL_vols{kk},:)).^2,1));
    rms2=sqrt(mean((KRON_full_fwd{KRON_ses(kk)}(KRON_vols{kk},:)).^2,1));
    
    testdata = [rms1,rms2];
    rms_difference_stats{kk} = bramila_two_sample_test(testdata,...
        [ones(1,length(rms1)),2*ones(1,length(rms2))],...
        10000,N_WORKER,0,0);    
    
    fprintf('\n>>>>>> Motion rms test (%i vols, %i vs. %i): pval = %.3f\n\n',length(ORIGINAL_vols{kk}),length(rms1),length(rms2),2*min(rms_difference_stats{kk}.pvals));
        
    clear isc1 isc2 testdata rms1 rms2 minvols
    %% Compute ISC with full matrices
    
    % element_indices==1 for ORIG_data
    % element_indices==2 for KRON_data
    % element_indices==3 for mixed data
    [ISC_matrices{kk},element_indices,ORIGINAL_ISC_vals,KRON_ISC_vals,ORIGINAL_ISC_permdist,KRON_ISC_permdist] = compute_ISC(ORIGINAL_data,KRON_data,2e5,N_WORKER);
    
    if nnz(isnan(ISC_matrices{kk}(:)))>0
        warning('NaN values found!')
    end
    clear ORIGINAL_data KRON_data;
    
    fprintf('>>>>>> ISC results: ORIGINAL mean ISC %.3f (null mean %.3f), KRON mean ISC %.3f (null mean %.3f)\n\n',mean(ORIGINAL_ISC_vals),mean(ORIGINAL_ISC_permdist),mean(KRON_ISC_vals),mean(KRON_ISC_permdist));
    
    diary(DIARY_NAME);
    
    %% ISC difference test (two separate groups, label-mixing permutation)
    
    testdata = ISC_matrices{kk}';
    testdata = [testdata(:,element_indices==1),testdata(:,element_indices==2)];
    testdata = atanh(testdata); % convert to the full real axis
    if nnz(isnan(testdata(:)))>0
        error('Bad data!');
    end
    
    ISC_difference_stats{kk} = bramila_two_sample_test(testdata,...
        [ones(1,nnz(element_indices==1)),2*ones(1,nnz(element_indices==2))],...
        10000,N_WORKER,0,0);
    
    ISC_difference_stats{kk}.twotail_pvals =  2*min(ISC_difference_stats{kk}.pvals');
    ISC_difference_stats{kk}.twotail_pvals_cor = mafdr(ISC_difference_stats{kk}.twotail_pvals,'BHFDR',true);    
    
    fprintf('\n>>>>>> ISC group-difference result: %i/%i significant at p<0.05\n',nnz(ISC_difference_stats{kk}.twotail_pvals_cor<0.05),length(ISC_difference_stats{kk}.twotail_pvals_cor));
    
    diary(DIARY_NAME);
    
    %% ISC perspective difference test (mixed group, Mantel test with model matrix)
    
    testdata = ISC_matrices{kk};
    testdata = 1 - testdata; % dissimilarity measure!
    
    perspective_model = nan(length(element_indices),1);
    perspective_model(element_indices==1)=0; % high similarity between ORIGINAL viewers
    perspective_model(element_indices==2)=0; % high similarity between KRON viewers
    perspective_model(element_indices==3)=1; % low similarity between mixed viewers
    
    % testdata = pairs x voxels
    % perspective_model = pairs x 1
    [ISC_perspective_corr{kk},ISC_perspective_pval{kk}]=bramila_mantel_vector(testdata,perspective_model,2e+5,'spearman',N_WORKER);        
    ISC_perspective_pval_cor{kk} = mafdr(ISC_perspective_pval{kk},'BHFDR',true);
    
    fprintf('\n>>>>>> ISC Mantel test result: %i/%i significant at p<0.05\n',nnz(ISC_perspective_pval_cor{kk}<0.05),length(ISC_perspective_pval_cor{kk}));    
    
    diary(DIARY_NAME);
    
    clear testdata;
    
    %% Subject classification test (quess which group)
    % the last cell-element is the mean results (for averaged ISC corrmat)
    group_ID = [ones(1,length(ORIGINAL_SUBJ)),2*ones(1,length(KRON_SUBJ))];
    [perspective_knn_accuracy{kk},perspective_knn_results{kk},perspective_knn_pval{kk}] = perspective_kNN_classification(ISC_matrices{kk},[5,7,9],group_ID,60000);   
    perspective_knn_pval_cor{kk}=mafdr(perspective_knn_pval{kk},'BHFDR',true);
    
    %% Save all results so far

    diary(DIARY_NAME);
    
    save('memento_segment_wise_ISC_preprocessor.mat','-v7.3');
    
end

diary(DIARY_NAME)

diary off;

fprintf('\nALL DONE!\n\n')



