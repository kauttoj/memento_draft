clc;
clear all;
close all;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/ISC/segment_wise');

%%
% load('memento_segment_wise_ISC_preprocessor.mat',...
% 'SELECTED_SEGMENTS',...
% 'MOVIE_SEGMENTS',...
% 'INTERPOLATION_METHOD',...
% 'N_WORKER',...
% 'KRON_SUBJ',...
% 'KRON_FILENAME',...
% 'KRON_fmri_time',...
% 'KRON_requested_time',...
% 'ORIG_SUBJ',...
% 'ORIG_FILENAME',...
% 'ORIG_fmri_time',...
% 'ORIG_requested_time',...
% 'group_mask_ind');
% 
% kk=0;
% for segment_nr = SELECTED_SEGMENTS   
%     kk=kk+1;
%     
%     fprintf('\n\n------------ Starting ISC calculations for segment %i/%i (ID %i) ---------------\n',kk,length(SELECTED_SEGMENTS),MOVIE_SEGMENTS(segment_nr).ID);
%         
%     CURRPATH = [HOMEPATH,filesep,sprintf('Seg_ID_%i',MOVIE_SEGMENTS(segment_nr).ID)];
%                       
%     funAlign_createpool(N_WORKER);
%     
%     % RUN PREPROCESSOR
%     tempdata = cell(1,length(KRON_SUBJ)); 
%     parfor k=1:length(KRON_SUBJ)
%         SUB = KRON_SUBJ{k};                 
%         nam = [KRON_EPI_PATH{KRON_ses(kk)},SUB,KRON_FILENAME];
%         fprintf('...subj %s: %s\n',SUB,nam);
%         nii=load_nii(nam,KRON_vols{kk});
%         data = nii.img;
%         data=reshape(data,[],size(data,4))';
%         data=data(:,group_mask_ind);
%         if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
%             error('data seems bad!')
%         end                
%         data = interp1(KRON_fmri_time{kk},data,KRON_requested_time{kk},INTERPOLATION_METHOD);
%         tempdata{k}=data;
%     end      
%     KRON_data = zeros(length(KRON_requested_time{kk}),length(group_mask_ind),length(KRON_SUBJ),'single');
%     for k=1:length(KRON_SUBJ)
%         KRON_data(:,:,k)=tempdata{k};
%     end   
%     clear tempdata data nii;
%     
%     diary(DIARY_NAME);
%     
%     tempdata = cell(1,length(ORIGINAL_SUBJ)); 
%     parfor k=1:length(ORIGINAL_SUBJ)
%         SUB = ORIGINAL_SUBJ{k};                 
%         nam = [ORIGINAL_EPI_PATH{ORIGINAL_ses(kk)},SUB,ORIGINAL_FILENAME];
%         fprintf('...subj %s: %s\n',SUB,nam);
%         nii=load_nii(nam,ORIGINAL_vols{kk});
%         save_nii(nii,[CURRPATH,filesep,SUB,'.nii']);             
%         save_text([CURRPATH,filesep,SUB,'_motion.txt'],ORIGINAL_motiondata{ORIGINAL_ses(kk),k}(ORIGINAL_vols{kk},:));
%         data=run_segment_preprocessor([CURRPATH,filesep,SUB,'.nii'],[CURRPATH,filesep,SUB,'_motion.txt'],CURRPATH,CFG);
%         delete([CURRPATH,filesep,SUB,'_motion.txt']);
%         if DELETE_SEGMENT_NII==1
%             delete([CURRPATH,filesep,SUB,'.nii']);
%         end
%         data=reshape(data,[],size(data,4))';
%         data=data(:,group_mask_ind);
%         if nnz(isnan(data))>0 || nnz(sum(data)==0)>0
%             error('data seems bad!')
%         end                
%         data = interp1(ORIGINAL_fmri_time{kk},data,ORIGINAL_requested_time{kk},INTERPOLATION_METHOD);
%         tempdata{k}=data;
%     end      
%     ORIGINAL_data = zeros(length(ORIGINAL_requested_time{kk}),length(group_mask_ind),length(ORIGINAL_SUBJ),'single');
%     for k=1:length(ORIGINAL_SUBJ)
%         ORIGINAL_data(:,:,k)=tempdata{k};
%     end   
%     clear tempdata data nii;
%         
%     %% Compute ISC with full matrices
%     
%     [ISC_matrices{kk},element_indices,ORIGINAL_ISC_vals,KRON_ISC_vals,ORIGINAL_ISC_permdist,KRON_ISC_permdist] = compute_ISC(ORIGINAL_data,KRON_data,2e5,N_WORKER);
% end



load('memento_segment_wise_ISC_preprocessor.mat','group_mask');
group_mask_ind = find(group_mask);

load('memento_segment_wise_ISC_preprocessor.mat','perspective_knn_accuracy','ISC_perspective_corr');

T = 0;
for i=1:length(ISC_matrices)
    
    t1 = length(ORIGINAL_fmri_time{i});
    t2 = length(KRON_fmri_time{i});
    
    if t1~=t2
        error('Segment has inconsistent length!')
    end
    
    T = T + t1;
    
    
    ISC_perspective_corr{i}(ind);
    perspective_knn_accuracy{i}(ind);
    
end

N_WORKER = 7;
try
    myCluster = gcp('nocreate');
    if isempty(myCluster)
        %delete(gcp)
        myCluster = parcluster('local');
        myCluster.NumWorkers=N_WORKER;
        parpool(myCluster);
    end
    N_workers = myCluster.NumWorkers;
catch err % old matlab?
    if ~matlabpool('size')
        eval(['matlabpool local ',num2str(N_WORKER)]);
    end
    N_workers = matlabpool('size');
end

load('memento_segment_wise_ISC_preprocessor.mat');

T = 0;
for i=1:length(ISC_matrices)
    
    t1 = length(ORIGINAL_fmri_time{i});
    t2 = length(KRON_fmri_time{i});
    
    if t1~=t2
        error('Segment has inconsistent length!')
    end
    
    T = T + t1;
    
    if i==1        
        mean_ISC_mat = t1*ISC_matrices{i};
    else
        mean_ISC_mat = mean_ISC_mat + t1*ISC_matrices{i};
    end
    
end

mean_ISC_mat = mean_ISC_mat/T;

testdata = mean_ISC_mat;
testdata = 1 - testdata; % dissimilarity measure!

perspective_model = nan(length(element_indices),1);
perspective_model(element_indices==1)=0; % high similarity between ORIGINAL viewers
perspective_model(element_indices==2)=0; % high similarity between KRON viewers
perspective_model(element_indices==3)=1; % low similarity between mixed viewers

% testdata = pairs x voxels
% perspective_model = pairs x 1
fprintf('Starting mean perspective ISC permutations...\n');
[mean_perspective_knn_accuracy,mean_perspective_knn_results,mean_perspective_knn_pval] = perspective_kNN_classification(mean_ISC_mat,[3,5,7],group_ID,1e+6);
mean_perspective_knn_pval_cor=mafdr(mean_perspective_knn_pval,'BHFDR',true);

[mean_ISC_perspective_corr,mean_ISC_perspective_pval]=bramila_mantel_vector(testdata,perspective_model,1e+6,'spearman',N_WORKER);
mean_ISC_perspective_pval_cor = mafdr(mean_ISC_perspective_pval,'BHFDR',true);


save('mean_perspective_ISC_results.mat','mean_ISC_mat','mean_ISC_perspective_corr','mean_ISC_perspective_pval','mean_ISC_perspective_pval_cor','element_indices','group_mask','mean_perspective_knn_accuracy','mean_perspective_knn_pval','mean_perspective_knn_pval_cor')

%%
addpath('/triton/becs/scratch/braindata/kauttoj2/code/spm12/');

load('memento_segment_wise_ISC_preprocessor.mat','group_mask');
group_mask_ind = find(group_mask);
%%

a=0*group_mask;
i=mean_ISC_perspective_pval_cor<0.01;
a(group_mask_ind(i))=mean_ISC_perspective_corr(i);
a=extentThreshold(a,20);
save_nii_oma(a,'mean_ISC_p001fdr_k20.nii');

a=0*group_mask;
i=mean_perspective_knn_pval_cor<0.01;
a(group_mask_ind(i))=mean_perspective_knn_accuracy(i);
a=extentThreshold(a,20);
save_nii_oma(a,'mean_kNN_p001fdr_k20.nii');

%%

a=0*group_mask;
i=mean_ISC_perspective_pval_cor<0.05;
a(group_mask_ind(i))=mean_ISC_perspective_corr(i);
a=extentThreshold(a,20);
save_nii_oma(a,'mean_ISC_p005fdr_k20.nii');

a=0*group_mask;
i=mean_perspective_knn_pval_cor<0.05;
a(group_mask_ind(i))=mean_perspective_knn_accuracy(i);
a=extentThreshold(a,20);
save_nii_oma(a,'mean_kNN_p005fdr_k20.nii');

%%


fprintf('ALL DONE!\n');
