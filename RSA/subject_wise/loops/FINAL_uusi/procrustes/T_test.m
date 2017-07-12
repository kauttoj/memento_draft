function [Quantile, p, T, df]=T_test(p1,p2)
addpath(genpath('statbx42'))

n1=size(p1,3);
n2=size(p2,3);
[mean_shape1 pc_shape s1 ]=tangent_pca_shape(p1);
[mean_shape2 pc_shape s2 ]=tangent_pca_shape(p2);
s1=sum(s1.^2);
s2=sum(s2.^2);

[distance]= finddistance(mean_shape1,mean_shape2);

total_var=sqrt(s1./n1+s2./n2);
T=distance./total_var;
df=((s1./n1+s2./n2)^2)./((s1./n1)^2./(n1-1)+(s2./n2)^2./(n2-1));

Quantile(1)=tq(0.975,df);
Quantile(2)=tq(0.025,df);
p=2*(1-tp(T,df));
