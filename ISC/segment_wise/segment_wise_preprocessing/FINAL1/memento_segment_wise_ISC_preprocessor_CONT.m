clear all;
clc;

fprintf('Continuing computations, loading data\n');

load('memento_segment_wise_ISC_preprocessor.mat');

INDS = SELECTED_SEGMENTS(segment_nr:end);
kk=kk-1;

for segment_nr = INDS
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
    [perspective_knn_accuracy{kk},perspective_knn_results{kk},perspective_knn_pval{kk}] = perspective_kNN_classification(ISC_matrices{kk},[5,7,9],group_ID,50000);   
    perspective_knn_pval_cor{kk}=mafdr(perspective_knn_pval{kk},'BHFDR',true);
    
    %% Save all results so far
    
    save('memento_segment_wise_ISC_preprocessor.mat','-v7.3');
    
end

diary(DIARY_NAME)

diary off;

fprintf('\nALL DONE!\n\n')

% if DELETE_SEGMENT_NII==0
%    
% fprintf('\n Starting segment joining process\n');
% kk=0;
% for segment_nr = SELECTED_SEGMENTS   
%     kk=kk+1;    
%     fprintf('\n\n------------ Adding segment %i/%i (ID %i) ---------------\n',kk,length(SELECTED_SEGMENTS),MOVIE_SEGMENTS(segment_nr).ID);
%                      
%     for k=1:length(ORIGINAL_SUBJ)
%         if kk==1
%             
%         else
%         
%         end                   
%         SUB = ORIGINAL_SUBJ{k};
%         nii = load_nii([CURRPATH,filesep,SUB,'_JOINED.nii']);
%     end
%     
%     
% end     
% end


