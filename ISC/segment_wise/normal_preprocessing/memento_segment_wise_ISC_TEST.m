clc;
clear all;
close all;

HOMEPATH = pwd;

addpath('/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila');
addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/ISC/segment_wise');
%% SEGMENT TIMING DATA

%load('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/repeating_segment_data_stage2.mat','full_corrdata');

cd('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');
TIMING_file = 'jelena_segment_timings_corrected_stage3.txt';
TIMING_data = parse_timing_file(TIMING_file);
cd(HOMEPATH);

%SEGMENT_TIMINGS = load('???????');

%selected_segments = [

% original:	636 - 930	(video1)
% 		201 - 462	(video2)
% 		1498 - 2107	(video3)
% kron:		1312 -1606	(video3)
% 		0 - 261		(video3)
% 		1385 - 1994	(video1)

START_SHIFT = 8; % seconds
END_SHIFT = 4; % seconds

SELECTED_IDs=16;

MOVIE_SEGMENTS(length(SELECTED_IDs))=struct();

k=0;
for i=1:length(full_corrdata)
    if ismember(full_corrdata(i).ID,SELECTED_IDs)
        k=k+1;
        MOVIE_SEGMENTS(k).ID = full_corrdata(i).ID;
        
        MOVIE_SEGMENTS(k).ORIGINAL_times = [full_corrdata(i).beginning_best_time(1)+START_SHIFT,full_corrdata(i).ending_best_time(1)+END_SHIFT];
        MOVIE_SEGMENTS(k).ORIGINAL_ses = full_corrdata(i).ORIGINAL_ses;
        
        MOVIE_SEGMENTS(k).KRON_times = [full_corrdata(i).beginning_best_time(2)+START_SHIFT,full_corrdata(i).ending_best_time(2)+END_SHIFT];
        MOVIE_SEGMENTS(k).KRON_ses = full_corrdata(i).KRON_ses;
        
        if abs(diff(diff(MOVIE_SEGMENTS(k).ORIGINAL_times)-diff(MOVIE_SEGMENTS(k).KRON_times)))>0.25
            error('Bad timing!')
        end
        
    end    
end

clear full_corrdata;

% MOVIE_SEGMENTS(1).ORIGINAL_times = [636,930];
% MOVIE_SEGMENTS(1).ORIGINAL_ses = 1;
% MOVIE_SEGMENTS(1).KRON_times = [1312,1606];
% MOVIE_SEGMENTS(1).KRON_ses = 3;
% 
% MOVIE_SEGMENTS(2).ORIGINAL_times = [201,462];
% MOVIE_SEGMENTS(2).ORIGINAL_ses = 2;
% MOVIE_SEGMENTS(2).KRON_times = [1e-10,261];
% MOVIE_SEGMENTS(2).KRON_ses = 3;
% 
% MOVIE_SEGMENTS(3).ORIGINAL_times = [1498,2107];
% MOVIE_SEGMENTS(3).ORIGINAL_ses = 3;
% MOVIE_SEGMENTS(3).KRON_times = [1385,1994];
% MOVIE_SEGMENTS(3).KRON_ses = 1;

%% INPUT DATA PATHS
%Memento3part_run1.avi
%Memento3part_run2.avi
%Memento3part_run3.avi

%ORIGINAL_SUBJ = {'S5','S7','S8','S9','S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'};
ORIGINAL_SUBJ = {'S8','S9'};%,'S12','S13','S15','S16','S17','S19','S21','S22','S23'};
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
KRON_SUBJ = {'KRON_2','KRON_3'};%,'KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13'};
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
N_WORKER = 4;
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
group_mask(:,35:end,:)=0;
group_mask(:,:,40:end)=0;
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

kk=0;
for segment_nr = SELECTED_SEGMENTS   
    kk=kk+1;
    
    fprintf('\n\n------------ Starting ISC calculations for segment %i/%i ---------------\n',kk,length(SELECTED_SEGMENTS));
    
%     CURRPATH = [HOMEPATH,filesep,'Segment_',num2str(kk)];
%     if ~exist(CURRPATH,'dir')
%        mkdir(CURRPATH);
%     end    
    
    % ORIGINAL
    fprintf('\nReading ORIGINAL data\n');     
    SES = MOVIE_SEGMENTS(segment_nr).ORIGINAL_ses;
    movie_times = MOVIE_SEGMENTS(segment_nr).ORIGINAL_times;    
    vol_min = floor(movie_times(1)/TR);
    vol_max = ceil(movie_times(2)/TR);    
    vols = (vol_min:vol_max) +1;    
    if vols(end)>ORIGINAL_max_volumes(SES)
        warning('Number of volumes too large!')
        vols(vols>ORIGINAL_max_volumes(SES))=[];
    end    
    fmri_time = (vols-1)*TR;
    requested_times = movie_times(1):TR_isc:movie_times(2);
    requested_times(requested_times>=fmri_time(end))=[];
    N_vols1=length(requested_times);
    if ~(fmri_time(1)<requested_times(1) && fmri_time(end)>requested_times(end))
        error('times are incorrect!')
    end
    
  
    fprintf('ORIGINAL: ses %i, vols %i-%i, fmri time %.2f-%.2f, seq. time %.2f-%.2f, interp. times %.2f-%.2f\n',SES,vols(1),vols(end),fmri_time(1),fmri_time(end),movie_times(1),movie_times(end),requested_times(1),requested_times(end));
    
    k=0;
    for SUB = ORIGINAL_SUBJ
        k=k+1;          
        nam = [ORIGINAL_EPI_PATH{SES},SUB{1},ORIGINAL_FILENAME];
        fprintf('...subj %s: %s\n',SUB{1},nam);
        nii=load_nii(nam,vols);
        data = nii.img;
        siz=size(data);
        data=reshape(data,[],size(data,4))';
        data=data(:,group_mask_ind);
        if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
            error('data seems bad!')
        end                
        data = interp1(fmri_time,data,requested_times,INTERPOLATION_METHOD);
        if k==1
            ORIGINAL_data=zeros(length(requested_times),length(group_mask_ind),length(ORIGINAL_SUBJ));
        end        
        ORIGINAL_data(:,:,k)=data;
    end   
    diary(DIARY_NAME);
    
    % KRON
    fprintf('\nReading KRON data\n');
    SES = MOVIE_SEGMENTS(segment_nr).KRON_ses;
    movie_times = MOVIE_SEGMENTS(segment_nr).KRON_times;
    vol_min = floor(movie_times(1)/TR);
    vol_max = ceil(movie_times(2)/TR);    
    vols = (vol_min:vol_max) +1;        
    if vols(end)>KRON_max_volumes(SES)
        warning('Number of volumes too large!')
        vols(vols>KRON_max_volumes(SES))=[];
    end
    fmri_time = (vols-1)*TR;
    requested_times = movie_times(1):TR_isc:movie_times(2);
    requested_times(requested_times>=fmri_time(end))=[];
    N_vols2=length(requested_times);            
    if ~(fmri_time(1)<requested_times(1) && fmri_time(end)>requested_times(end))
        error('times are incorrect!')
    end
    
    
    minvols=min(size(original_fwd,1),size(kron_fwd,1));
    kron_fwd = kron_fwd(1:minvols,:);
    original_fwd = original_fwd(1:minvols,:);
    
    if N_vols1~=N_vols2
        error('Number of volumes does not match (boundary effect?)!')
    end        
    
    fprintf('KRON: ses %i, vols %i-%i, fmri time %.2f-%.2f, seq. time %.2f-%.2f, interp. time %.2f-%.2f\n',SES,vols(1),vols(end),fmri_time(1),fmri_time(end),movie_times(1),movie_times(end),requested_times(1),requested_times(end));
    
    k=0;
    for SUB = KRON_SUBJ
        k=k+1;          
        nam = [KRON_EPI_PATH{SES},SUB{1},KRON_FILENAME];
        fprintf('...subj %s: %s\n',SUB{1},nam);
        nii=load_nii(nam,vols);
        data = nii.img;
        siz=size(data);
        data=reshape(data,[],size(data,4))';
        data=data(:,group_mask_ind);
        if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
            error('data seems bad!')
        end                
        data = interp1(fmri_time,data,requested_times,INTERPOLATION_METHOD);
        if k==1
            KRON_data=zeros(length(requested_times),length(group_mask_ind),length(KRON_SUBJ));
        end        
        KRON_data(:,:,k)=data;        
    end             
    diary(DIARY_NAME);
    
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
    [ISC_matrices{kk},element_indices,ORIGINAL_ISC_vals,KRON_ISC_vals,ORIGINAL_ISC_permdist,KRON_ISC_permdist] = compute_ISC(ORIGINAL_data,KRON_data,1e5,N_WORKER);
    
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
    [ISC_perspective_corr{kk},ISC_perspective_pval{kk}]=bramila_mantel_vector(testdata,perspective_model,100000,'spearman',N_WORKER);        
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


