clc;
clear all;
close all;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/ISC/segment_wise');


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

[mean_perspective_knn_accuracy,mean_perspective_knn_results,mean_perspective_knn_pval] = perspective_kNN_classification(mean_ISC_mat,[3,5,7],group_ID,2e+6);
mean_perspective_knn_pval_cor=mafdr(mean_perspective_knn_pval,'BHFDR',true);

[mean_ISC_perspective_corr,mean_ISC_perspective_pval,mean_ISC_perspective_nulldist]=bramila_mantel_vector(testdata,perspective_model,2e+6,'spearman',N_WORKER);
mean_ISC_perspective_pval_cor = mafdr(mean_ISC_perspective_pval,'BHFDR',true);

save('mean_perspective_ISC_results.mat','mean_ISC_mat','mean_ISC_perspective_corr','mean_ISC_perspective_pval','mean_ISC_perspective_pval_cor','element_indices','group_mask','mean_perspective_knn_accuracy','mean_perspective_knn_pval','mean_perspective_knn_pval_cor','mean_ISC_perspective_nulldist','-v7.3');

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
save_nii_oma(a,'mean_ISC_p005fdr_k1.nii');
a=extentThreshold(a,20);
save_nii_oma(a,'mean_ISC_p005fdr_k50.nii');

a=0*group_mask;
i=mean_perspective_knn_pval_cor<0.05;
a(group_mask_ind(i))=mean_perspective_knn_accuracy(i);
save_nii_oma(a,'mean_kNN_p005fdr_k1.nii');
a=extentThreshold(a,20);
save_nii_oma(a,'mean_kNN_p005fdr_k50.nii');

%%


fprintf('ALL DONE!\n');
