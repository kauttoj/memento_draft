% find distance for two shapes
function [distance F2 temp tran]= finddistance(dataN,dataN1)
% dataN and dataN1 N*k  two preshapes
% N is the number of points, and k is the dimension

%output
%distance is the geodestic distance between dataN and dataN1;

[dataN cN scaleN]=pre_shape(dataN);
[dataN1 cN1 scaleN1]=pre_shape(dataN1);
C=dataN'*dataN1;
[S,V,D]=svd(C);
tran=S*transpose(D);
F=S*transpose(D)*dataN1';% the final minimizer
F=F';
F2=F/norm(F,'fro');
temp=sum(diag(V));
distance = real(acos(temp));
end
