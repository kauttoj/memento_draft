%% demo code for voxel-wise HRF deconvolution
%% From NIFTI image (resting state fMRI data) to NIFTI image (HRF parameters).
%% Guo-Rong Wu, gronwu@gmail.com, UESTC, UGent, 2013.9.12
%% Reference: Wu, G.; Liao, W.; Stramaglia, S.; Ding, J.; Chen, H. & Marinazzo, D.. 
%% A blind deconvolution approach to recover effective connectivity brain networks 
%% from resting state fMRI data. Medical Image Analysis, 2013,17(3):365-374 .
clc;
clear all;
close all;
%% Mask file
num_voxel = 10;
nobs = 500; % number of time points
bsig = zeros(nobs,num_voxel); 

TR = 1.56; %in seconds 
thr = 1.5; % threshold, for example 1 SD.
event_lag_max = 7; % the (estimated) maximum lagged time from neural event to BOLD event, in points. 

%for isub=1:length(sub)

load rsig;

%S12_roi_data.mat;
%rsig = data{2}{130};

disp('Retrieving HRF ...');
tic
close all;
for i=1:length(rsig_vis)
    rsig = rsig_vis{i};
    [data_deconv,onset,hrf,event_lag,PARA] = wgr_deconv_canonhrf_par(rsig,thr,event_lag_max,TR,nworker);
    plot(((1:size(hrf,1))-1)*TR,hrf);
    pause(1);
end
close all;
for i=1:length(rsig_vis)
    rsig = rsig_high{i};
    [data_deconv,onset,hrf,event_lag,PARA] = wgr_deconv_canonhrf_par(rsig,thr,event_lag_max,TR,nworker);
    plot(((1:size(hrf,1))-1)*TR,hrf);
    pause(1);    
end
toc
disp('Done');

%save(fullfile(save_dir,[sub(isub).name,'_hrf.mat']),'event_lag_max','thr','TR','onset', 'hrf','event_lag', 'PARA','-v7.3');
    
    
    % Write HRF parameter
    % Height - h
    % Time to peak - p (in time units of TR seconds)
    % FWHM (at half peak) - w  

    v=spm_vol(brainmask);
    v.dt=[16,0]; 
    
    v.fname = fullfile(save_dir,'Height',[sub(isub).name,'_height.nii']);
    data = data_tmp;
    data(voxel_ind)=PARA(1,:);
    spm_write_vol(v,data);

    v.fname = fullfile(save_dir,'T2P',[sub(isub).name,'_Time2peak.nii']);
    data = data_tmp;
    data(voxel_ind)=PARA(2,:);
    spm_write_vol(v,data);

    v.fname = fullfile(save_dir,'FWHM',[sub(isub).name,'_FWHM.nii']);
    data = data_tmp;
    data(voxel_ind)=PARA(3,:);
    spm_write_vol(v,data);
    
    
    sub_save_dir = fullfile(save_dir2,sub(isub).name);
    mkdir(sub_save_dir)
    % writting back into nifti files
    for k = 1:length(imag)
        v.fname = fullfile(sub_save_dir,imag(k).name);
        data = data_tmp;
        data(voxel_ind) = data_deconv(k,:);
        spm_write_vol(v,data);
    end   
    
