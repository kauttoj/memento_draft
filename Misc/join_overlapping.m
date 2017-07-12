function models = join_overlapping(models,min_dt)

N=length(models);
for i=1:N
   T = length(models(i).onset);
   k=1;
   while k<T
       k=k+1;
       if models(i).onset(k)-models(i).ending(k-1)<min_dt
           models(i).ending(k-1) = models(i).ending(k);
           models(i).duration(k-1) = models(i).ending(k-1)-models(i).onset(k-1);
           models(i).ending(k)=[];
           models(i).onset(k)=[];
           models(i).duration(k)=[];
           T=T-1;
           k=k-1;
       end
   end    
end


end