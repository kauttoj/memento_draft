function [p_new]=p2preshape(p,ii,a)
% p is the preshape n*k*m
% ii which PC
% a 
[pci W p_mean sigma]=get_PC(ii,p);
% project back to preshape space
 p_new=cos(a*sqrt(sigma))*p_mean+(sin(a*sqrt(sigma))/norm(pci,'fro'))*pci;
 distance =finddistance(p_mean,p_new);
 abs(a)*sqrt(sigma);
end