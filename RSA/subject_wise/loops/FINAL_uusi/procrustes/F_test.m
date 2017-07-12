function [Quantile, p,F]=F_test(p1,p2)
addpath(genpath('statbx42'))

m=size(p1,3);
n=size(p2,3);
[mean_shape1 pc_shape Sx ]=tangent_pca_shape(p1);
[mean_shape2 pc_shape Sy ]=tangent_pca_shape(p2);
Sx=sum(Sx.^2);
Sy=sum(Sy.^2);
F=(Sx./Sy);

Quantile(1)=fq(0.975,m-1,n-1);
Quantile(2)=fq(0.025,m-1,n-1);
p=2*min(fp(F,m-1,n-1),1-fp(F,m-1,n-1));
