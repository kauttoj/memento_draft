function [res ,n]= convert_str2num(responses,choice_str)

N=length(responses);
L=length(choice_str);
res = nan(1,N);
for i=1:N
   for j=1:L
       if isempty(choice_str{j})
           if isempty(responses{i})
               res(i)=j;
           end
       elseif strfind(responses{i},choice_str{j})
           res(i)=j;
       end
   end
end

n=0;
if nnz(isnan(res))>0
%    warning('%i NaN''s found (bad data)',nnz(isnan(res)))
    n = nnz(isnan(res));
end

end