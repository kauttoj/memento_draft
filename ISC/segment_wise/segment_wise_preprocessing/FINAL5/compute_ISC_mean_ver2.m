clc;
clear all;
close all;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/ISC/segment_wise');
addpath('/m/nbe/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/ISC/segment_wise');

load('memento_segment_wise_ISC_preprocessor.mat','group_mask');
group_mask_ind = find(group_mask);

COORD = [...
    -10,30,24;...
    12,-66,34;...
    46,-58,52;...
    -38,22,36;...
    -6,-94,10;...
    60,-16,-2;...
    ];
COORD1 = COORD;
COORD = mni2ind(COORD);
COORD_IND = sub2ind(size(group_mask),COORD(:,1),COORD(:,2),COORD(:,3));

temp = 0*group_mask;
temp(COORD_IND)=1;
findClusterLabels_ver2(temp,1,20)

LABELS{1}='L cingulate gyrus';
LABELS{2}='R precuneous cortex';
LABELS{3}='R lat. occipital cortex';
LABELS{4}='L mid. front. gyrus';
LABELS{5}='L occipital pole';
LABELS{6}='R sup. temporal gyrus';

load('memento_segment_wise_ISC_preprocessor.mat',...
    'perspective_knn_accuracy','ISC_perspective_corr',...
    'ORIGINAL_fmri_time','KRON_fmri_time','MOVIE_SEGMENTS',...
    'ORIGINAL_requested_time','KRON_requested_time','fwd_ISC_difference_stats','rms_difference_stats',...
    'ISC_perspective_pval','perspective_knn_pval');

T = 0;
for i=1:length(ISC_perspective_corr)
    
    t11 = length(ORIGINAL_fmri_time{i});
    t22 = length(KRON_fmri_time{i});
        
    if t11~=t22
        error('Segment has inconsistent length!')
    end        
    
    t1 = ORIGINAL_requested_time{i};
    t2 = KRON_requested_time{i};
    
    t1 = t1(end)-t1(1);
    t2 = t2(end)-t2(1);
    
    if abs(t1 - t2)>1e-4
        error('!!!!')
    end
    
    time(i)=mean([t1,t2]);           
    
    for j=1:length(COORD_IND)                            
        
        if length(ISC_perspective_corr{i})~=length(group_mask_ind)
            error('!!!!')
        end
        
        k = find(COORD_IND(j)==group_mask_ind);
        
        RSA_val(i,j)=ISC_perspective_corr{i}(k);
        kNN_val(i,j)=perspective_knn_accuracy{i}(k);
        
        RSA_val_pval(i,j)=ISC_perspective_pval{i}(k);
        kNN_val_pval(i,j)=perspective_knn_pval{i}(k);
        
    end
    
    ID(i)=MOVIE_SEGMENTS(i).ID;
    
end

TIMELIMITS = [100,150,200,250];


for BIG_LOOP = 1:length(TIMELIMITS)
    
    close all;
    
    LIMIT = TIMELIMITS(BIG_LOOP);
    ind = time>=LIMIT;
    handles=tight_subplot(2,3, [0.05,0.04],[0.15 0.05],[0.15 0.05]);
    
    for j=1:length(COORD_IND)
        set(gcf, 'CurrentAxes',handles(j));
        
        a = RSA_val(ind,j);
        a(:,2) = kNN_val(ind,j);
        a(:,3)=a(:,1).*a(:,2);
        
        b = RSA_val_pval(ind,j);
        b(:,2) = kNN_val_pval(ind,j);
        b(:,3)=b(:,1) < 0.05 & b(:,2) < 0.05;        
        
        %a(:,2) = a(:,2)-0.5;
        
        bar(a(:,3));
        set(gca,'XTick',1:i,'XTickLabel',ID(ind));
        axis tight;
        set(gca,'XTickLabel',[]);
        
        x = [get(gca,'XLim'),get(gca,'YLim')];
        set(gca,'YLim',[x(3),x(4)+0.015]);        
        
        for k=1:nnz(ind)
            if b(k,3)>0
                text(k-0.25,a(k,3)+0.006,'*');
            end
        end
        
        s = sprintf('%s\n(%i,%i,%i)',LABELS{j},COORD1(j,1),COORD1(j,2),COORD1(j,3));
        text(x(1)+1,x(4)-(x(4)-x(3))*0.1,s,'FontSize',12);
        
        set(gca,'XTick',1:i,'XTickLabel',ID(ind),'Fontsize',10);
        %xticklabel_rotate([],45,[],'Fontsize',14);
        
    end
    
    set(gcf, 'CurrentAxes',handles(1));
    ylabel('kNN\timesRSA')
    set(gcf, 'CurrentAxes',handles(4));
    ylabel('kNN\timesRSA')
    
    for k=[4,5,6]
        set(gcf, 'CurrentAxes',handles(k));
        xlabel('Segment ID','FontSize',12)
    end
    %set(gca,'XTick',1:i,'XTickLabel',ID(ind));
      
    
    disp('');
    
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

TIMELIMITS = [100,150,200,250];
BIG_LOOP = 0;

for BIG_LOOP = 1:length(TIMELIMITS)
    
    LIMIT = TIMELIMITS(BIG_LOOP);
    
    STR = [num2str(LIMIT),'_'];
    
    T = 0;
    mean_ISC_mat=0;
    for i=1:length(ISC_perspective_corr)
        
        t11 = length(ORIGINAL_fmri_time{i});
        t22 = length(KRON_fmri_time{i});
        
        if t11~=t22
            error('Segment has inconsistent length!')
        end
        
        t1 = ORIGINAL_requested_time{i};
        t2 = KRON_requested_time{i};
        
        t1 = t1(end)-t1(1);
        t2 = t2(end)-t2(1);
        
        if abs(t1 - t2)>1e-4
            error('!!!!')
        end
        
        time(i)=mean([t1,t2]);
        
        if time(i)>= LIMIT
            
            T = T + time(i);
            
            mean_ISC_mat = mean_ISC_mat + time(i)*ISC_matrices{i};
            
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
    
    [mean_perspective_knn_accuracy,mean_perspective_knn_results,mean_perspective_knn_pval,mean_perspective_knn_nulldist] = perspective_kNN_classification(mean_ISC_mat,[3,5,7],group_ID,1e+6);
    mean_perspective_knn_pval_cor=mafdr(mean_perspective_knn_pval,'BHFDR',true);
    
    [mean_ISC_perspective_corr,mean_ISC_perspective_pval,mean_ISC_perspective_nulldist]=bramila_mantel_vector(testdata,perspective_model,1e+6,'spearman',N_WORKER);
    mean_ISC_perspective_pval_cor = mafdr(mean_ISC_perspective_pval,'BHFDR',true);
    
    save([STR,'mean_perspective_ISC_results.mat'],'mean_ISC_mat','mean_ISC_perspective_corr','mean_ISC_perspective_pval','mean_ISC_perspective_pval_cor','element_indices','group_mask','mean_perspective_knn_accuracy','mean_perspective_knn_pval','mean_perspective_knn_pval_cor','mean_ISC_perspective_nulldist','mean_perspective_knn_nulldist','-v7.3')
    
    %%
    addpath('/triton/becs/scratch/braindata/kauttoj2/code/spm12/');
    
    %%
    
    a=0*group_mask;
    i=mean_ISC_perspective_pval_cor<0.01;
    a(group_mask_ind(i))=mean_ISC_perspective_corr(i);
    save_nii_oma(a,[STR,'mean_ISC_p001fdr_k1.nii']);
    a=extentThreshold(a,20);
    save_nii_oma(a,[STR,'mean_ISC_p001fdr_k20.nii']);
    
    a=0*group_mask;
    i=mean_perspective_knn_pval_cor<0.01;
    a(group_mask_ind(i))=mean_perspective_knn_accuracy(i);
    save_nii_oma(a,[STR,'mean_kNN_p001fdr_k1.nii']);
    a=extentThreshold(a,20);
    save_nii_oma(a,[STR,'mean_kNN_p001fdr_k20.nii']);
    
    %%
    
    a=0*group_mask;
    i=mean_ISC_perspective_pval_cor<0.05;
    a(group_mask_ind(i))=mean_ISC_perspective_corr(i);
    save_nii_oma(a,[STR,'mean_ISC_p005fdr_k1.nii']);
    a1=a;
    a=extentThreshold(a,20);
    save_nii_oma(a,[STR,'mean_ISC_p005fdr_k20.nii']);
    
    a=0*group_mask;
    i=mean_perspective_knn_pval_cor<0.05;
    a(group_mask_ind(i))=mean_perspective_knn_accuracy(i);
    save_nii_oma(a,[STR,'mean_kNN_p005fdr_k1.nii']);
    a2=a;
    a=extentThreshold(a,20);
    save_nii_oma(a,[STR,'mean_kNN_p005fdr_k20.nii']);
    
    save_nii_oma(a1.*a2,[STR,'mean_kNN_and_ISC_p005fdr_k1.nii']);
    a=extentThreshold(a1.*a2,20);
    save_nii_oma(a,[STR,'mean_kNN_and_ISC_p005fdr_k20.nii']);
    
    %%
end

fprintf('ALL DONE!\n');
