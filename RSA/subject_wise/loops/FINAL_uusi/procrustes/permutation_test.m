function [dist p]=permutation_test(p1,p2, iter)
addpath(genpath('statbx42'))

n1=size(p1,3);
n2=size(p2,3);
mean_shape1=mean_shape(p1);
mean_shape2=mean_shape(p2);
[dist]= finddistance(mean_shape1,mean_shape2);
p_all=cat(3,p1,p2);
dist1=zeros(1,iter);
for i=1:iter
    index_all=randsample(n1+n2, n1+n2,'true');
    p1_test=p_all(:,:,index_all(1:n1));
    p2_test=p_all(:,:,index_all(n1+1:end));
    
    mean_shape_test1=mean_shape(p1_test);
    mean_shape_test2=mean_shape(p2_test);
    [dist1(i)]= finddistance(mean_shape_test1,mean_shape_test2);
end
p=numel(find(dist1>dist))./numel(dist1);

