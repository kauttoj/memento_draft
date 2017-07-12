clc;
clear all;
close all;

load volume_selection_data.mat
load roi_and_mask_data.mat
load combined_nullvals_deconvolved.mat

N_rois = length(my_masks);
ROI_id = 1:N_rois;

selection_threshold = 99.9;

for i=1:N_rois
   threshold(i) = prctile(all_nullvals_deconvolved(i,:),99);
   threshold1(i) = prctile(all_nullvals_deconvolved(i,:),selection_threshold);
   n_voxels(i) = length(my_masks{i});
end

mean_threshold = mean(threshold);
scale=threshold/mean_threshold;

mean_correlations = zeros(length(DELAYS),N_rois);

for subj_ind = 1:length(S);
    subj = S{subj_ind};
    load([subj,'_correlation_results.mat']);
    for j=1:N_rois
        mean_correlations(:,j) = mean_correlations(:,j) + all_correlations(:,j);
    end
end
mean_correlations=mean_correlations/length(S);

good_rois=[];
for j=1:N_rois
    %if nnz((mean_correlations(:,j))>0.04)>1%nnz(mean_correlations(DELAYS<7,j)>threshold1(j))>1
    if nnz(mean_correlations(DELAYS>-5 & DELAYS<6,j)>threshold1(j))>0
        good_rois(end+1)=j;
    end
    mean_correlations(:,j)=mean_correlations(:,j)/scale(j);
    threshold_check(j)=threshold(j)/scale(j);
end
plot(DELAYS,ones(size(DELAYS))*mean_threshold,'k--','LineWidth',2);hold on;
plot(DELAYS,mean_correlations(:,good_rois));
axis tight;
title(sprintf('ROIs with p<%.2s correlations (%i)',(100-selection_threshold)/100,length(good_rois)));
xlabel('Delay')
ylabel('Correlation')
box on;
    
ind = find(DELAYS>-5 & DELAYS<6);
SELECTED_DELAYS = DELAYS(ind);
selected_data = mean_correlations(ind,good_rois)';
figure;
plot(SELECTED_DELAYS,ones(size(SELECTED_DELAYS))*mean_threshold,'k--','LineWidth',2);hold on;
plot(SELECTED_DELAYS,selected_data');
axis tight;
title(sprintf('ROIs with p<%.2s correlations (%i)',(100-selection_threshold)/100,length(good_rois)));
xlabel('Delay')
ylabel('Correlation')
box on;

%%

% figure;
%T1 = clusterdata(selected_data,'maxclust',2,'distance','correlation','linkage','average');
%[silhCos,h] = silhouette(selected_data,T1,'correlation');

figure;
[T2,cmeans,sumd] = kmeans(selected_data,2,'replicates',200,'display','final','distance','correlation');
[silval,h] = silhouette(selected_data,T2,'correlation');
w=[];
for i=unique(T2')
    j=find(T2==i);
    v = silval(j);
    [v,k]=sort(v);
    w=[w;j(k)];
end
T2=T2(w);
silval=silval(w);

ROI_id=ROI_id(good_rois);
ROI_id=ROI_id(w);
my_rois=my_rois(good_rois);
my_rois=my_rois(w);
selected_data = selected_data(w,:);

[~,ind]=sort(T2);
cormat = corr(selected_data');
%cormat = cormat(w,w);
cormat=cormat-diag(diag(cormat));
plot_matrix_simple(cormat,'cormat');

figure;
subplot(2,1,1);hold on;
plot(DELAYS,ones(size(DELAYS))*mean_threshold,'k--','LineWidth',2);hold on;
plot(DELAYS,mean_correlations(:,good_rois(w(T2==1))));
axis tight;
title(sprintf('ROIs with p<%.2s correlations, cluster 1 (%i)',(100-selection_threshold)/100,nnz(T2==1)));
%xlabel('Delay')
ylabel('Correlation')
set(gca,'XGrid','on')
box on;

subplot(2,1,2);hold on;
plot(DELAYS,ones(size(DELAYS))*mean_threshold,'k--','LineWidth',2);hold on;
plot(DELAYS,mean_correlations(:,good_rois(w(T2==2))));
axis tight;
title(sprintf('ROIs with p<%.2s correlations, cluster 2 (%i)',(100-selection_threshold)/100,nnz(T2==2)));
xlabel('Delay')
ylabel('Correlation')
set(gca,'XGrid','on')
box on;

peak_delay=[];
for i=1:size(selected_data,1)
    [~,k]=max(selected_data(i,:));
    peak_delay(i)=SELECTED_DELAYS(k);
end
bad = [];%peak_delay==6;

peak_delay(bad)=[];
selected_data(bad,:)=[];
ROI_id(bad)=[];
my_rois(bad)=[];

figure;
[yy,xx]=hist(peak_delay,SELECTED_DELAYS);
bar(xx,yy);
xlabel('Delay')
ylabel('Maximum peak count')
axis tight;
title(sprintf('Maxima location for ROIs with p<%.2s correlations (%i)',(100-selection_threshold)/100,length(ROI_id)));
set(gca,'XTick',SELECTED_DELAYS)

l=find(yy,1,'first');
h=find(yy,1,'last');

cmap = jet(h-l+1);
col = peak_delay - (SELECTED_DELAYS(l)) + 1;
addpath('/triton/becs/scratch/braindata/kauttoj2/code/BrainNetViewer');
H_BrainNet = BrainNet_MapCfg_custom(my_rois,ones(length(my_rois),1),{cmap,col},[],[],[]);
%H_BrainNet = BrainNet_MapCfg_custom(my_rois,ones(length(my_rois),1),T2,[],[],[]);

figure;
plot(1:size(cmap,1),l:h);
cbar = colorbar;
colormap(cmap);
set(cbar,'YTick',linspace(1/(2*size(cmap,1)),(2*size(cmap,1)-1)/(2*size(cmap,1)),size(cmap,1)));
set(cbar,'FontSize',14)
a=num2cellstr(xx(l):xx(h));
for i=1:length(a)
    a{i}=[a{i},'s'];
end
set(cbar,'YTickLabel',a)

% Gaussian mixture
% data=[randn(300,1);randn(200,1)*2-8];GMModel = fitgmdist(data,2);
% [y1,x1]=hist(data,30);
% x2=linspace(min(data),max(data),100);y2=pdf(GMModel,x2');y2=y2*max(y1)/max(y2);
% plot(x1,y1,x2,y2)

%%
% 
% Dist = pdist(selected_data','correlation');
% Tree = linkage(Dist,'average');
% leafOrder = optimalleaforder(Tree,Dist);
% [h,nodes] = dendrogram(Tree,length(good_rois),'Reorder',leafOrder,'ColorThreshold','default');
% 
% cormat = corr(selected_data);
% cormat = cormat(leafOrder,leafOrder);

    %subplot(length(ROIs_of_interest),1,i);hold on;