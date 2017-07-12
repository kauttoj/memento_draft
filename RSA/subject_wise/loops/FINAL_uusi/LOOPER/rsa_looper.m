
clc;
clear all;
close all;

home = pwd;
addpath(pwd);

HARD_LIMIT = 3;
delay_array = -26:2:20;
radius_array = 6;
counts = -ones(length(radius_array),length(delay_array),2);

save('memento_LOOPER_results.mat','counts','radius_array','delay_array');

k0=0;
for radius = radius_array
    k0=k0+1;
    k2=0;
    for delay = delay_array
        k2=k2+1;
        
        cd(home)
        s = sprintf('delay_%g_radius_%i',delay,radius);
        mkdir(s);
        folder = [home,filesep,s];
        
        fprintf('\n\n ---------- RUNNING BIG LOOP (%i,%i) of (%i,%i) ------------- \n\n',k0,k2,length(radius_array),length(delay_array));
        
        prepare_data_ver6_looper(folder,delay,HARD_LIMIT);
        
        cd(home)
        [counts(k0,k2,1),counts(k0,k2,2)] = memento_rsa_ver5_looper(folder,radius);
        
        cd(home)
        save('memento_LOOPER_results.mat','counts','radius_array','delay_array');
        
    end
end

cd(home)