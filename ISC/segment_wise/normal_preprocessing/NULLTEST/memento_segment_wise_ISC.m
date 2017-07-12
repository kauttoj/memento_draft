clc;
clear all;
close all;

HOMEPATH = pwd;

addpath('/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila');
addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/ISC/segment_wise');
addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');

HOMEPATH = pwd;
CURRPATH = HOMEPATH;
%% SEGMENT TIMING DATA

TIMING_FILE = [pwd,filesep,'memento_null_timings.txt'];
MOVIE_SEGMENTS = parse_timing_file(TIMING_FILE);

START_SHIFT = 10; % seconds
END_SHIFT = 5; % seconds

SELECTED_IDs=1:5;

ind=[];
k=0;
for i=1:length(MOVIE_SEGMENTS)
    MOVIE_SEGMENTS(i).ORIGINAL_time = MOVIE_SEGMENTS(i).ORIGINAL_time + [START_SHIFT,END_SHIFT];
    MOVIE_SEGMENTS(i).KRON_time=MOVIE_SEGMENTS(i).KRON_time + [START_SHIFT,END_SHIFT];
    if ismember(MOVIE_SEGMENTS(i).ID,SELECTED_IDs)        
        if diff(MOVIE_SEGMENTS(i).ORIGINAL_time)>150                       
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

%ORIGINAL_SUBJ = {'S5','S7','S8','S9','S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'};
ORIGINAL_SUBJ = {'S7','S9','S10','S12','S13','S15','S16','S17','S21','S22','S23'};
%ORIGINAL_SUBJ = {'S8','S9','S12','S15','S17','S21','S22','S23'};
%ORIGINAL_SUBJ = {'S8','S9','S12'};

ORIGINAL_EPI_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/post_processed_normal/bramila/';
ORIGINAL_EPI_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/post_processed_normal/bramila/';
ORIGINAL_EPI_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/post_processed_normal/bramila/';
ORIGINAL_MOTION_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/RealignParameter/';
ORIGINAL_MOTION_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/RealignParameter/';
ORIGINAL_MOTION_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/RealignParameter/';
ORIGINAL_FILENAME = '_mask_detrend_fullreg_filtered_smoothed.nii';
ORIGINAL_MASK_FILENAME = '_analysis_mask.nii';
ORIGINAL_MOTION_FILENAME = '/4D_motion_params.txt';

%KRON_SUBJ = {'KRON_2','KRON_3','KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13'};
KRON_SUBJ = {'KRON_2','KRON_3','KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13'};
%KRON_SUBJ = {'KRON_5','KRON_6','KRON_7','KRON_8'};

KRON_EPI_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session1/post_processed_normal/bramila/';
KRON_EPI_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session2/post_processed_normal/bramila/';
KRON_EPI_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session3/post_processed_normal/bramila/';
KRON_MOTION_PATH{1} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session1/RealignParameter/';
KRON_MOTION_PATH{2} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session2/RealignParameter/';
KRON_MOTION_PATH{3} = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/KRON/Session3/RealignParameter/';
KRON_FILENAME = '_mask_detrend_fullreg_filtered_smoothed.nii';
KRON_MASK_FILENAME = '_analysis_mask.nii';
KRON_MOTION_FILENAME = '/4D_motion_params.txt';

%% Parameters

SELECTED_SEGMENTS = 1:length(MOVIE_SEGMENTS);
INTERPOLATION_METHOD = 'linear';
TR = 1.56;
DIARY_NAME = 'memento_segment_wise_ISC_analysis_diary.txt';
N_WORKER = 3;
TR_isc = (TR/2);

delete(DIARY_NAME);
diary(DIARY_NAME)

%% Create group mask

group_mask = true(91,109,91);
for SES=1:3
    k=0;
    for SUB = ORIGINAL_SUBJ
        k=k+1;
        nam = [ORIGINAL_EPI_PATH{SES},SUB{1},ORIGINAL_MASK_FILENAME];
        nii=load_nii(nam);
        mask = nii.img>0;
        fprintf('mask size %i (%s)\n',nnz(mask),nam);
        group_mask = group_mask.*mask;
    end
    k=0;
    for SUB = KRON_SUBJ
        k=k+1;
        nam = [KRON_EPI_PATH{SES},SUB{1},KRON_MASK_FILENAME];
        nii=load_nii(nam);
        mask = nii.img>0;
        fprintf('mask size %i (%s)\n',nnz(mask),nam);        
        group_mask = group_mask.*mask;
    end
    fprintf('\n');
end
group_mask=group_mask>0;
%---TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%group_mask(:,35:end,:)=0;
%group_mask(:,:,40:end)=0;
%---TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fprintf('>>>>>> GROUP MASK SIZE %i \n',nnz(group_mask));
save_nii_oma(group_mask,'memento_segment_wise_ISC_mask.nii');
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
        original_full_rms(SES,k)=mean(rms);
        original_full_fwd{SES}(:,k)=fwd(1:ORIGINAL_max_volumes(SES));
    end
end

for SES=1:3    
    kron_full_fwd{SES} = nan(KRON_max_volumes(SES),length(KRON_SUBJ));
    k=0;
    for SUB = KRON_SUBJ
        k=k+1;
        cfg.motionparam = [KRON_MOTION_PATH{SES},SUB{1},KRON_MOTION_FILENAME]; 
        [fwd,rms]=bramila_framewiseDisplacement(cfg);
        kron_full_rms(SES,k)=mean(rms);
        kron_full_fwd{SES}(:,k)=fwd(1:KRON_max_volumes(SES));
    end
end

mat_ind_original = find(triu(ones(length(ORIGINAL_SUBJ),length(ORIGINAL_SUBJ)),1));
mat_ind_kron = find(triu(ones(length(KRON_SUBJ),length(KRON_SUBJ)),1));

fprintf('\n\nTesting group differences of RMS and FWD ISC\n')
for SES=1:3
    testdata = [original_full_rms(SES,:),kron_full_rms(SES,:)];    
    rms_full_difference_stats{SES} = bramila_two_sample_test(testdata,...
        [ones(1,size(original_full_rms,2)),2*ones(1,size(kron_full_rms,2))],...
        10000,N_WORKER,0,0);
    fprintf('\n>>>>>> Session %i motion RMS: pval = %.3f\n\n',SES,2*min(rms_full_difference_stats{SES}.pvals)); 
    
    m = min(size(original_full_fwd{SES},1),size(kron_full_fwd{SES},1));
    
    isc1 = corr(original_full_fwd{SES}(1:m,:));
    isc1 = isc1(mat_ind_original);
    isc2 = corr(kron_full_fwd{SES}(1:m,:));
    isc2 = isc2(mat_ind_kron);
    
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
    SES = MOVIE_SEGMENTS(segment_nr).ORIGINAL_ses;    
    movie_times = MOVIE_SEGMENTS(segment_nr).ORIGINAL_time;    
    fmri_times = ((1:ORIGINAL_max_volumes(SES)) -1)*TR;            
    vol_min = find(movie_times(1)>fmri_times,1,'last');
    vol_max = find(movie_times(2)<fmri_times,1,'first');       
    vols = (vol_min:vol_max);         
    fmri_times = fmri_times(vols);
    requested_times = movie_times(1):TR_isc:movie_times(2);
    requested_times(requested_times>=fmri_times(end))=[];    
    N_vols1=length(requested_times);
    if ~(fmri_times(1)<requested_times(1) && fmri_times(end)>requested_times(end))
        error('times are incorrect!')
    end
    ORIGINAL_requested_times{kk}=requested_times;
    ORIGINAL_vols{kk} = vols;  
    ORIGINAL_fmri_time{kk}=fmri_times;
    fprintf('ORIGINAL film, ID %i (scan session %i):\n    sequence time %.2f-%.2f (requested %.2f-%.2f)\n    fmri times %.2f-%.2f (vols %i-%i)\n',...
        MOVIE_SEGMENTS(segment_nr).ID,...
        SES,...
        movie_times(1)-START_SHIFT,...
        movie_times(end)-END_SHIFT,...  
        requested_times(1),...
        requested_times(end),...
        fmri_times(1),...
        fmri_times(end),...        
        vols(1),...
        vols(end));
       
    % KRON
    SES = MOVIE_SEGMENTS(segment_nr).KRON_ses;    
    movie_times = MOVIE_SEGMENTS(segment_nr).KRON_time;    
    fmri_times = ((1:KRON_max_volumes(SES)) -1)*TR;            
    vol_min = find(movie_times(1)>fmri_times,1,'last');
    vol_max = find(movie_times(2)<fmri_times,1,'first');       
    vols = (vol_min:vol_max);         
    fmri_times = fmri_times(vols);
    requested_times = movie_times(1):TR_isc:movie_times(2);
    requested_times(requested_times>=fmri_times(end))=[];    
    N_vols2=length(requested_times);
    if ~(fmri_times(1)<requested_times(1) && fmri_times(end)>requested_times(end))
        error('times are incorrect!')
    end
    KRON_requested_times{kk}=requested_times;
    KRON_vols{kk} = vols;  
    KRON_fmri_time{kk}=fmri_times;
    fprintf('KRON film, ID %i (scan session %i):\n    sequence time %.2f-%.2f (requested %.2f-%.2f)\n    fmri times %.2f-%.2f (vols %i-%i)\n',...
        MOVIE_SEGMENTS(segment_nr).ID,...
        SES,...
        movie_times(1)-START_SHIFT,...
        movie_times(end)-END_SHIFT,...  
        requested_times(1),...
        requested_times(end),...
        fmri_times(1),...
        fmri_times(end),...        
        vols(1),...
        vols(end));    
    if N_vols1~=N_vols2
       error('Volume count is different!')
    end                    
    clear fmri_times requested_times movie_times SES vols;
end


kk=0;
for segment_nr = SELECTED_SEGMENTS   
    kk=kk+1;
    
    fprintf('\n\n------------ Starting ISC calculations for segment %i/%i ---------------\n',kk,length(SELECTED_SEGMENTS));
    
    SES = MOVIE_SEGMENTS(segment_nr).ORIGINAL_ses;   
    fprintf('\nReading ORIG data\n'); 
    k=0;
    for SUB = ORIGINAL_SUBJ
        k=k+1;          
        nam = [ORIGINAL_EPI_PATH{SES},SUB{1},ORIGINAL_FILENAME];
        fprintf('...subj %s: %s\n',SUB{1},nam);
        nii=load_nii(nam,ORIGINAL_vols{kk});
        data = nii.img;
        siz=size(data);
        data=reshape(data,[],size(data,4))';
        data=data(:,group_mask_ind);
        if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
            error('data seems bad!')
        end                
        data = interp1(ORIGINAL_fmri_time{kk},data,ORIGINAL_requested_times{kk},INTERPOLATION_METHOD);
        if k==1
            ORIGINAL_data=zeros(length(ORIGINAL_requested_times{kk}),length(group_mask_ind),length(ORIGINAL_SUBJ),'single');
        end        
        ORIGINAL_data(:,:,k)=data;
    end   
    original_fwd = original_full_fwd{SES}(ORIGINAL_vols{kk},:); 
    diary(DIARY_NAME);
    
    % KRON
    SES = MOVIE_SEGMENTS(segment_nr).KRON_ses;  
    fprintf('\nReading KRON data\n');           
    k=0;
    for SUB = KRON_SUBJ
        k=k+1;          
        nam = [KRON_EPI_PATH{SES},SUB{1},KRON_FILENAME];
        fprintf('...subj %s: %s\n',SUB{1},nam);
        nii=load_nii(nam,KRON_vols{kk});
        data = nii.img;
        siz=size(data);
        data=reshape(data,[],size(data,4))';
        data=data(:,group_mask_ind);
        if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
            error('data seems bad!')
        end                
        data = interp1(KRON_fmri_time{kk},data,KRON_requested_times{kk},INTERPOLATION_METHOD);
        if k==1
            KRON_data=zeros(length(KRON_requested_times{kk}),length(group_mask_ind),length(KRON_SUBJ),'single');
        end        
        KRON_data(:,:,k)=data;        
    end              
    kron_fwd = kron_full_fwd{SES}(KRON_vols{kk},:);
    diary(DIARY_NAME);
    
    minvols=min(size(original_fwd,1),size(kron_fwd,1));
    kron_fwd = kron_fwd(1:minvols,:);
    original_fwd = original_fwd(1:minvols,:);  
            
    clear data;
    clear nii;
    
    %% Compute fwd difference stats
    
    isc1 = corr(original_fwd);
    isc1 = isc1(mat_ind_original);
    isc2 = corr(kron_fwd);
    isc2 = isc2(mat_ind_kron);  
    
    testdata = [isc1',isc2'];
    testdata = atanh(testdata); % convert to the full real axis
    
    fprintf('\n');
    fwd_ISC_difference_stats{kk} = bramila_two_sample_test(testdata,...
        [ones(1,length(isc1)),2*ones(1,length(isc2))],...
        10000,N_WORKER,0,0);
    
    fprintf('\n>>>>>> Motion fwd ISC test (%i vols, %i vs. %i): pval = %.3f\n',minvols,length(isc1),length(isc2),2*min(fwd_ISC_difference_stats{kk}.pvals));       
    
    rms1=sqrt(mean(original_fwd.^2,1));
    rms2=sqrt(mean(kron_fwd.^2,1));
    
    testdata = [rms1,rms2];
    rms_difference_stats{kk} = bramila_two_sample_test(testdata,...
        [ones(1,length(rms1)),2*ones(1,length(rms2))],...
        10000,N_WORKER,0,0);    
    
    fprintf('\n>>>>>> Motion rms test (%i vols, %i vs. %i): pval = %.3f\n\n',minvols,length(rms1),length(rms2),2*min(rms_difference_stats{kk}.pvals));
        
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
    [ISC_perspective_corr{kk},ISC_perspective_pval{kk}]=bramila_mantel_vector(testdata,perspective_model,200000,'spearman',N_WORKER);        
    ISC_perspective_pval_cor{kk} = mafdr(ISC_perspective_pval{kk},'BHFDR',true);
    
    fprintf('\n>>>>>> ISC Mantel test result: %i/%i significant at p<0.05\n',nnz(ISC_perspective_pval_cor{kk}<0.05),length(ISC_perspective_pval_cor{kk}));    
    
    diary(DIARY_NAME);
    
    clear testdata;
    
    %% Save all results so far
    
    save('memento_segment_wise_ISC.mat','-v7.3');
    
end

diary(DIARY_NAME)

fprintf('\n\nALL DONE!\n')

diary off;


