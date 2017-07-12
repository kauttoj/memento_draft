function res = save_data(data,filename)

res = 0;
try
save(filename,'data','-v7.3');
catch err   
   return;
end

res=1;

end