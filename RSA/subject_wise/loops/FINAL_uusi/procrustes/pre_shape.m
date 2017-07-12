% Find the Pre_shape for any arbitary shape
function [p c scale]=pre_shape(s)
% s is the shape sample in the n*k, n points and each point lies in R^k
c=mean(s);
p=s-repmat(c,size(s,1),1);% centerlize the data.
scale=norm(p,'fro');
p=p/scale;
end
