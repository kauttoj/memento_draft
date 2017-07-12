clc;
clear all;
close all;

HOME = pwd;

cd('/m/nbe/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA/subject_wise/linear_nulltest/ROI_comparison_FINAL/KRON');

KRON_vols = load('volume_selection_data.mat');
KRON_res = load('roi_and_mask_data.mat');
%load combined_nullvals.mat

cd('/m/nbe/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA/subject_wise/linear_nulltest/ROI_comparison_FINAL/ORIG');

ORIG_vols = load('volume_selection_data.mat');
ORIG_res = load('roi_and_mask_data.mat');

cd(HOME);

common_map = zeros(91,109,91);
ORIG_map = common_map;
KRON_map = common_map;

KRON_rois_count = length(KRON_res.my_rois);
ORIG_rois_count = length(ORIG_res.my_rois);

for i = 1:length(ORIG_res.my_rois)
    x = ORIG_res.my_rois(i).centroid(1);
    y = ORIG_res.my_rois(i).centroid(2);
    z = ORIG_res.my_rois(i).centroid(3);
    ORIG_size(i)= length(ORIG_res.my_masks{i});
    ORIG_map(x,y,z)=i;
    common_map(x,y,z)=common_map(x,y,z)+1;
end

for i = 1:length(KRON_res.my_rois)
    x = KRON_res.my_rois(i).centroid(1);
    y = KRON_res.my_rois(i).centroid(2);
    z = KRON_res.my_rois(i).centroid(3);
    KRON_size(i)= length(KRON_res.my_masks{i});
    KRON_map(x,y,z)=i;
    common_map(x,y,z)=common_map(x,y,z)+1;
end

ORIG_keep = ORIG_map(common_map==2);
KRON_keep = KRON_map(common_map==2);

ORIG_res.my_rois = ORIG_res.my_rois(ORIG_keep);
ORIG_res.my_masks = ORIG_res.my_masks(ORIG_keep);
ORIG_size = ORIG_size(ORIG_keep);

KRON_res.my_rois = KRON_res.my_rois(KRON_keep);
KRON_res.my_masks = KRON_res.my_masks(KRON_keep);
KRON_size = KRON_size(KRON_keep);

for i = 1:length(KRON_res.my_rois)
    if nnz(KRON_res.my_rois(i).map - ORIG_res.my_rois(i).map)>0
        error('!!!!')
    end
end

val = abs(KRON_size-ORIG_size)./((KRON_size+ORIG_size)/2);
keep =val<0.10;

ORIG_res.my_rois = ORIG_res.my_rois(keep);
ORIG_res.my_masks = ORIG_res.my_masks(keep);
ORIG_size = ORIG_size(keep);

KRON_res.my_rois = KRON_res.my_rois(keep);
KRON_res.my_masks = KRON_res.my_masks(keep);
KRON_size = KRON_size(keep);

for i = 1:length(KRON_res.my_rois)
    if nnz(KRON_res.my_rois(i).map - ORIG_res.my_rois(i).map)>0
        error('!!!!')
    end
end

N_rois = length(ORIG_res.my_rois);

if nnz(ORIG_vols.DELAYS - KRON_vols.DELAYS)>0
    error('!!!')
end
DELAYS = ORIG_vols.DELAYS;
% for i=1:N_rois
%    threshold(i) = prctile(all_nullvals(i,:),99);
%    threshold1(i) = prctile(all_nullvals(i,:),99.99);
%    n_voxels(i) = length(my_masks{i});
% end
% 
% mean_threshold = mean(threshold);
% scale=threshold/mean_threshold;
% 

cd('/m/nbe/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA/subject_wise/linear_nulltest/ROI_comparison_FINAL/ORIG');

correlations = nan(length(ORIG_vols.S),ORIG_rois_count,length(ORIG_vols.DELAYS));

for subj_ind = 1:length(ORIG_vols.S);
    subj = ORIG_vols.S{subj_ind};
    load([subj,'_correlation_results.mat']);
    for j=1:ORIG_rois_count
        correlations(subj_ind,j,:) = all_correlations(:,j);
    end
end

if nnz(isnan(correlations))>0
   error('!!!') 
end

all_correlation=correlations(:,ORIG_keep,:);
all_correlation=all_correlation(:,keep,:);

ORIG_correlations = all_correlation;
clear correlation;

cd('/m/nbe/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA/subject_wise/linear_nulltest/ROI_comparison_FINAL/KRON');

correlations = nan(length(KRON_vols.S),KRON_rois_count,length(KRON_vols.DELAYS));

for subj_ind = 1:length(KRON_vols.S);
    subj = KRON_vols.S{subj_ind};
    load([subj,'_correlation_results.mat']);
    for j=1:KRON_rois_count
        correlations(subj_ind,j,:) = all_correlations(:,j);
    end
end

if nnz(isnan(correlations))>0
   error('!!!') 
end

all_correlation=correlations(:,KRON_keep,:);
all_correlation=all_correlation(:,keep,:);

KRON_correlations = all_correlation;
clear correlation;

cd(HOME);

%%

ORIG_correlations = atanh(ORIG_correlations);
ORIG_mean_correlations=squeeze(mean(ORIG_correlations,1));

KRON_correlations = atanh(KRON_correlations);
KRON_mean_correlations=squeeze(mean(KRON_correlations,1));

clear p_vals t_vals;
for i=1:size(ORIG_correlations,3)
   [h,p_vals(:,i),~,b]=ttest2(squeeze(ORIG_correlations(:,:,i)),squeeze(KRON_correlations(:,:,i)),'Vartype','unequal','tail','both');
   t_vals(:,i)=b.tstat;
end

DF = size(ORIG_correlations,1) + size(KRON_correlations,1) - 2;
th_001 = fzero(@(x) 2 * tcdf(-abs(x),DF) -0.001,10);
th_001_ala = fzero(@(x) 2 * tcdf(-abs(x),DF) -0.001,-10);

TIME_ind = find(DELAYS>=-8 & DELAYS<=6);
SELECTED_DELAYS = DELAYS(TIME_ind);

p = p_vals(:,TIME_ind);
s = size(p);
p = p(:);
[p_cor,q_cor]=mafdr(p);
p_cor=reshape(p_cor,s);
q_cor=reshape(q_cor,s);
%good_rois = find(sum(p_vals(:,TIME_ind)<0.001,2)>0);
good_rois = find(sum(p_cor<0.01,2)>0);

plot(DELAYS,t_vals(good_rois,:)); hold on;
plot(DELAYS([1,end]),th_001*ones(1,2),'b--');
plot(DELAYS([1,end]),th_001_ala*ones(1,2),'b--');
axis tight;

leg = [];
k=0;
for i=good_rois'
    k=k+1;
    leg{k}=KRON_res.my_rois(i).better_label;
    leg{k}=[leg{k},' ',num2str(max(t_vals(i,TIME_ind)))];
end
legend(leg);
title(sprintf('ROIs with p<0.001 between %.1fs - %.1fs (total %i)',DELAYS(TIME_ind(1)),DELAYS(TIME_ind(end)),length(good_rois)));
xlabel('Onset delay [s]')
ylabel(sprintf('t_{%i}-value',DF))
selected_data = t_vals(good_rois,TIME_ind);

my_rois = KRON_res.my_rois(good_rois);

%%

[~,~,MNI_id,MNI_label] = embed_rois(my_rois,ones(91,109,91));
%MNI_id(ismember(MNI_id,[2,11]))=0;

figure;
leg={};
for u = unique(MNI_id)
    if u>0
        i = find(MNI_id==u);        
        
        if length(i)>0
                       
            a=nan(1,length(DELAYS));
            for j=1:length(DELAYS)           
                data1=squeeze(mean(ORIG_correlations(:,good_rois(i),j),2));
                data2=squeeze(mean(KRON_correlations(:,good_rois(i),j),2));                
                [~,~,~,b]=ttest2(data1,data2,'Vartype','unequal');                                
                a(j)=b.tstat;
            end
            
            plot(DELAYS,a);hold on;
            leg=[leg,[MNI_label{u},' (',num2str(length(i)),')']];
        end
        
    end
end
legend(leg,'location','nw');
xlabel('Onset delay [s]','FontSize',14);
ylabel(sprintf('t_{%i}-value',DF),'FontSize',14);
plot(DELAYS([1,end]),th_001*[1,1],'k-');
axis tight;
a = get(gca,'ylim');
plot([0,0],[a(1),a(2)],'k');

%%

N_comps = 2;

figure;
Z = linkage(selected_data,'average','spearman');
w = optimalleaforder(Z,pdist(selected_data,'spearman'));
[h,nodes] = dendrogram(Z,length(good_rois),'Reorder',w,'ColorThreshold','default');
T2 = cluster(Z,'maxclust',N_comps);
T2 = T2(w);
%silval=silval(wmy_rois=my_rois(w);

ORIG_correlations = ORIG_correlations(:,good_rois,:);
ORIG_correlations = ORIG_correlations(:,w,:);

KRON_correlations = KRON_correlations(:,good_rois,:);
KRON_correlations = KRON_correlations(:,w,:);

selected_data = selected_data(w,:);

cormat = corr(selected_data','type','spearman');
%cormat = cormat(w,w);
cormat=cormat-diag(diag(cormat));
plot_matrix_simple(cormat,'cormat');

tvals_arr = t_vals(good_rois,:);
tvals_arr = tvals_arr(w,:);

figure;

for n=1:N_comps
    
    subplot(N_comps,1,n);hold on;
    
    plot(DELAYS,tvals_arr(T2==n,:));hold on;
    plot(DELAYS,th_001*ones(1,length(DELAYS)),'k--','LineWidth',2);
    
    %[~,vec] = pca(a','NumComponents',1);
    a=tvals_arr(T2==n,:);
    plot(DELAYS,mean(a),'k','LineWidth',2);
    plot([0,0],[min(a(:)),max(a(:))],'k');
    
    axis tight;
    title(sprintf('Cluster %i (%i ROIs)',n,nnz(T2==n)));
    %xlabel('Delay')
    ylabel(sprintf('t_{%i}-value',DF),'FontSize',14)
    set(gca,'XGrid','on')
    box on;
    
end
xlabel('Onset delay [s]','FontSize',14)

mean_corvals = 0*T2;
new_ind = 0*T2;
peaks = 0*T2;
my_rois1=my_rois;
for i=1:max(T2)
    ind = find(T2==i);
    c = corr(selected_data(ind,:)','type','spearman');
    for j=1:length(ind)
        a=true(1,size(c,1));
        a(j)=false;
        mean_corvals(ind(j))=mean(c(j,a));
        [~,peaks(ind(j))] = max(selected_data(ind(j),:));
    end
    [~,ii]=sort(mean_corvals(ind),'descend');    
    mean_corvals(ind)=mean_corvals(ind(ii));
    my_rois1(ind)=my_rois(ind(ii));
    peaks(ind)=peaks(ind(ii));
end
fprintf('Region of interest\tx\ty\tz\tcluster\tISC\tpeak\n')
for i=1:length(T2)    
    fprintf('%s\t%i\t%i\t%i\t%i\t%.2f\t%.1f\n',my_rois1(i).better_label,my_rois1(i).centroidMNI(1),my_rois1(i).centroidMNI(2),my_rois1(i).centroidMNI(3),T2(i),mean_corvals(i),SELECTED_DELAYS(peaks(i)));
end


figure;
addpath('/m/nbe/scratch/braindata/kauttoj2/code/BrainNetViewer');
%H_BrainNet = BrainNet_MapCfg_custom(my_rois,a,{cmap,col},[],[],[]);
H_BrainNet = BrainNet_MapCfg_custom(my_rois,ones(1,length(T2)),T2,[],[],[]);



