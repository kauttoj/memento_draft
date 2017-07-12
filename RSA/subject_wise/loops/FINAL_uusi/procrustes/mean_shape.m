function p=mean_shape(pi)
% pi is m preshape sample, pi is in the dimension n *k    *m
% k is the dimension of the points usually k=2,3.
% output p is the frechet mean shape
err=1;
threshold=10^(-15);
for i=1:size(pi,3)
    pi(:,:,i)=pre_shape(pi(:,:,i));
end
p=pi(:,:,1);
k=1;
while(err>threshold&&k<50)
    vsum=zeros(size(pi,1),size(pi,2));
    lambdasum=0;
    for i=1:size(pi,3)
        [distance F2 temp]=finddistance(p,pi(:,:,i));
        if temp==1;
            temp=temp-10^(-8);
        end
        vsum=vsum-(distance/sqrt(1-temp^2))*F2;
        lambdasum=lambdasum-temp*distance/sqrt(1-temp^2);
    end
    ptemp=(sign(lambdasum)/norm(vsum,'fro'))*vsum;
    err=norm(ptemp-p,'fro');
    p=ptemp;
    k=k+1;
end

%fprintf('This uses %d steps(max 49)\n', k);


    
        
        
        
