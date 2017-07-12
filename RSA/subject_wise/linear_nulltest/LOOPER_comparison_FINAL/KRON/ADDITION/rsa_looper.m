
clc;
clear all;
close all;

home = pwd;
addpath(pwd);

HARD_LIMIT = 3;
delay_array = -30:2:-14;
radius_array = 6;
counts = -ones(length(radius_array),length(delay_array),2);

save('memento_LOOPER_results_KRON.mat','counts','radius_array','delay_array','-v7.3');

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
        all_cormaps{k0,k2} = SL_correlator(folder,radius);
        
        cd(home)
        save('memento_LOOPER_results_KRON.mat','all_cormaps','counts','radius_array','delay_array','-v7.3');
        
    end
end

cd(home)