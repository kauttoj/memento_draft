clc;
clear all;
close all;

N_WORKER = 5;

load('memento_segment_wise_ISC_preprocessor.mat');

for i=1:length(ISC_matrices)
    
    if i==1        
        mean_ISC_mat = ISC_matrices{i};
    else
        mean_ISC_mat = mean_ISC_mat + ISC_matrices{i};
    end
    
end

mean_ISC_mat = mean_ISC_mat/length(ISC_matrices);

testdata = mean_ISC_mat;
testdata = 1 - testdata; % dissimilarity measure!

perspective_model = nan(length(element_indices),1);
perspective_model(element_indices==1)=0; % high similarity between ORIGINAL viewers
perspective_model(element_indices==2)=0; % high similarity between KRON viewers
perspective_model(element_indices==3)=1; % low similarity between mixed viewers

% testdata = pairs x voxels
% perspective_model = pairs x 1
fprintf('Starting mean perspective ISC permutations...\n');
[mean_ISC_perspective_corr,mean_ISC_perspective_pval]=bramila_mantel_vector(testdata,perspective_model,1e+6,'spearman',N_WORKER);
mean_ISC_perspective_pval_cor = mafdr(mean_ISC_perspective_pval,'BHFDR',true);

[mean_perspective_knn_accuracy,mean_perspective_knn_results,mean_perspective_knn_pval] = perspective_kNN_classification(mean_ISC_mat,[5,7,9],group_ID,100000);
mean_perspective_knn_pval_cor=mafdr(mean_perspective_knn_pval,'BHFDR',true);

save('mean_perspective_ISC_results.mat','mean_ISC_mat','mean_ISC_perspective_corr','mean_ISC_perspective_pval','mean_ISC_perspective_pval_cor','element_indices','group_mask','mean_perspective_knn_accuracy','mean_perspective_knn_pval','mean_perspective_knn_pval_cor')

fprintf('ALL DONE!\n');
