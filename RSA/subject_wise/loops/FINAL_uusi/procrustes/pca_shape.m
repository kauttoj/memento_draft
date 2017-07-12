function [mean_shape pc_shape std_shape pc_projection]=pca_shape(new_pt)
% input new_pt in dimension n*3*m
% after alignment

% output:
% mean_shape in n*3 for the mean shape
% pc_shape in n*3*k the first k pcs for the group of shapes
% std_shape stardard deviation along 
m=size(new_pt,1);
d=size(new_pt,2);
n=size(new_pt,3);
pc_shape=zeros(m,d,n);
mean_shape=mean(new_pt,3);
new_pt1=reshape(new_pt,m*d,n)-repmat(reshape(mean(new_pt,3),m*d,1),1,n);
[S V D]=svd(new_pt1'*new_pt1./(n-1));
pc_shape_long=(new_pt1./sqrt(n-1))*S./repmat(sqrt(diag(V)'),m*d,1);
for i=1:n
    pc_shape(:,:,i)=reshape(pc_shape_long(:,i),m,d);
end
    
std_shape=sqrt(diag(V));
pc_projection=new_pt1'*pc_shape_long;