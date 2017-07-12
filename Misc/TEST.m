        
MAX_OFFSET = 1;
framerate = 25;

ORIGINAL_data(30/(1/framerate))=struct();
KRON_data(30/(1/framerate))=struct();

MAX_SHIFT = round(MAX_OFFSET/(1/framerate));

print_int = round(linspace(MAX_SHIFT/10,MAX_SHIFT,10));
kk=1;

fprintf(' ...running jitter analysis\n');

ind = 1:min(length(ORIGINAL_data),length(KRON_data));

z=0;
total_corr=nan(1,MAX_SHIFT*2+1);
for i=-MAX_SHIFT:MAX_SHIFT
    
    z=z+1;
    
    ind1 = ind;
    ind2 = ind;
    
    if i<0
        ind1(1:(-i))=[];
        ind2 = ind2(1:length(ind1));
    elseif i>0
        ind1((end-i+1):end)=[];
        ind2(1:i)=[];
    else
        ;
    end
    
    if length(ind1)~=length(ind2)
        error('!!!')
    end
    
    dd(z)=ind1(1)-ind2(1);
    
end


% %%
% ii=1;
% BACKWARDS = 3;
% FORWARDS = 15;
% RATEMOD = 0.9;
% 
% %%
% close all;
% DATA = play_video({full_corrdata(ii).ORIGINAL_file,full_corrdata(ii).KRON_file},...
%     [full_corrdata(ii).beginning_best_time(1)-BACKWARDS,full_corrdata(ii).beginning_best_time(1)+FORWARDS;...
%     full_corrdata(ii).beginning_best_time(2)-BACKWARDS,full_corrdata(ii).beginning_best_time(2)+FORWARDS],...
%     RATEMOD);
% 
% %%
% close all;
% DATA = play_video({full_corrdata(ii).ORIGINAL_file,full_corrdata(ii).KRON_file},...
%     [full_corrdata(ii).ending_best_time(1)-BACKWARDS,full_corrdata(ii).ending_best_time(1)+FORWARDS;...
%     full_corrdata(ii).ending_best_time(2)-BACKWARDS,full_corrdata(ii).ending_best_time(2)+FORWARDS],...
%     RATEMOD);
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% close all;
% DATA = play_video({full_corrdata(ii).ORIGINAL_file,full_corrdata(ii).KRON_file},...
%     [31*60+43.80-BACKWARDS,31*60+43.80+FORWARDS;...
%     13*60+6.68-BACKWARDS,13*60+6.68+FORWARDS],...
%     RATEMOD);