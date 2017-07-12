function [mean_shape pc_shape std_shape pc_projection new_pt_out]=tangent_pca_shape(new_pt)
% input new_pt in dimension n*3*m

% output:
% mean_shape in n*3 for the mean shape
% pc_shape in n*3*k the first k pcs for the group of shapes
% std_shape stardard deviation along 

m=size(new_pt,1);
d=size(new_pt,2);
n=size(new_pt,3);
pc_shape=zeros(m,d,n);
[W mean_shape new_pt1]=p2tangent(new_pt);
[S V D]=svd(new_pt1'*new_pt1./(n-1));
pc_shape_long=(new_pt1./sqrt(n-1))*S./repmat(sqrt(diag(V)'),m*d,1);
for i=1:n
    pc_shape(:,:,i)=reshape(pc_shape_long(:,i),m,d);
    new_pt_out(:,:,i)=pre_shape(new_pt(:,:,i));
end
    
std_shape=sqrt(diag(V));
pc_projection=new_pt1'*pc_shape_long;


function [W p Z]=p2tangent(pi)
% pi is m preshape sample, pi is in the dimension n *k *m
% k is the dimension of the points usually k=2,3.
% output p is the frechet mean shape
% output W is the vector in the tangent space with the length of
% ditance(p,pi(:,:,i))
p=mean_shape(pi);
for i=1:size(pi,3)
    [distance F2 temp]= finddistance(p,pi(:,:,i));
    W(:,:,i)=(distance/norm(F2-temp*p,'fro'))*(F2-temp*p);
    Z(:,i)=reshape(W(:,:,i),size(W,1)*size(W,2),1);
end