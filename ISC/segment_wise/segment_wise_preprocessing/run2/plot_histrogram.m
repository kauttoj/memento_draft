clc;
clear all;
close all;

load('memento_segment_wise_ISC_preprocessor.mat','ISC_perspective_pval','ISC_perspective_pval_cor','MOVIE_SEGMENTS');

clear count_corr count_noncorr seg

for k = 1:length(ISC_perspective_pval)        
    count_noncorr(k)=nnz(ISC_perspective_pval{k}<0.001);
    count_corr(k)=nnz(ISC_perspective_pval_cor{k}<0.05);
    seg(k)=MOVIE_SEGMENTS(k).ID;
end

load('mean_perspective_ISC_results.mat');

count_corr(end+1)=nnz(mean_ISC_perspective_pval_cor<0.05);
count_noncorr(end+1)=nnz(mean_ISC_perspective_pval<0.001);

a=num2cellstr(seg);
a{end+1}='mean';

figure;
bar(count_corr);
axis tight;
set(gca,'XTick',1:length(count_corr),'XTickLabel',a);
title('ISC perspective differences (p<0.05 FDR)')

figure;
bar(count_noncorr);
axis tight;
set(gca,'XTick',1:length(count_noncorr),'XTickLabel',a);
title('ISC perspective differences (p<0.001 uncorr.)')