
clc;
clear all;
close all;

home = pwd;
addpath(pwd);

duration_array = 4;
delay_array = -20:2.5:20;
radius_array = 5;
counts = -ones(length(radius_array),length(duration_array),length(delay_array),2);

save('memento_LOOPER_results.mat','counts','radius_array','duration_array','delay_array');

k0=0;
for radius = radius_array
    k0=k0+1;
    k1=0;
    for duration = duration_array
        k1=k1+1;
        k2=0;        
        for delay = delay_array
            k2=k2+1;
            
            cd(home)            
            s = sprintf('duration_%i_delay_%g_radius_%i',duration,delay,radius);
            mkdir(s);
            folder = [home,filesep,s];
            
            fprintf('\n\n ---------- RUNNING BIG LOOP (%i,%i,%i) of (%i,%i,%i) ------------- \n\n',k0,k1,k2,length(radius_array),length(duration_array),length(delay_array));
            
            prepare_data_ver3_looper(folder,delay,duration);
            
            cd(home)
            [counts(k0,k1,k2,1),counts(k0,k1,k2,2)] = memento_rsa_ver5_looper(folder,radius);
            
            cd(home)
            save('memento_LOOPER_results.mat','counts','radius_array','duration_array','delay_array');
            
        end
    end
end

cd(home)
