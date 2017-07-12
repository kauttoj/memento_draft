clear all;close all;clc;

load('memento_segment_wise_ISC_preprocessor.mat','ISC_perspective_pval','ISC_perspective_pval_cor');

N=length(ISC_perspective_pval);

cor=zeros(N,N);
for i=1:N,
    map1=ISC_perspective_pval{i}<0.05;    
    for j=(i+1):N,
        map2=ISC_perspective_pval{j}<0.05;
        cor(i,j)=nnz(map1 & map2)/nnz(map1 | map2);
    end
end
cor = cor+cor';