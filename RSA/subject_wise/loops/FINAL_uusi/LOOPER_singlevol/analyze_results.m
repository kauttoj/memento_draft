clc;
clear variables;
close all;

%delay_array = -26:2:20;
radius = 5;

% counter = 0;
% for d=delay_array,
%        
%     s = sprintf('delay_%g_radius_%i',d,radius);
%     fprintf('Delay %s\n',s);
%     cd(s);
%     load memento_rsa_results.mat;
%     counter = counter+1;
%     th_p005fdr(counter) = Th_kauppi(2);
%     kauppi_voxels(counter)=nnz(mean_map>Th_kauppi(2));
%     
%     cd Maps
%     fil = dir('rs_subject*.nii');
%     img = zeros(91,109,91,length(fil));
%     fprintf('  %i subjects\n',length(fil));
%     for i=1:length(fil)
%        nii=load_nii(fil(i).name); 
%        img(:,:,:,i)=nii.img;
%     end
%     cd ..
%     cd ..
%     save_nii_oma(img,[s,'_corrmaps.nii']);    
% end

N_workers = create_workers(3);

output_txt=cell(1,length(delay_array));
output_files=cell(1,length(delay_array));

parfor i=1:length(delay_array),  
    d = delay_array(i);
    s = sprintf('delay_%g_radius_%i',d,radius);
    if d>-15 && d<15        
        [output_files{i},output_txt{i}] = randomise_ttest([pwd,filesep,s,'_corrmaps.nii'],[pwd,filesep,'supermask.nii'],1000,'fsl5.0-');
    end
end

