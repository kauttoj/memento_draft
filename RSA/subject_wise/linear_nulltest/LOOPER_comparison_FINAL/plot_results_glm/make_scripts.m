clc;
clear variables;
close all;

delays = -8:2:6;

for d=delays
        
    s = sprintf('delay%.1f',d);    

    cd POS   
    nii_p=load_nii(sprintf('delay_%.1f_zscoremats_tfce_corrp_tstat1.nii',d));
    nii_t=load_nii(sprintf('delay_%.1f_zscoremats_tstat1.nii',d));    
    img = nii_t.img.*(nii_p.img>0.99); 
    fprintf('voxels to show: %i\n',nnz(img>0));    
    mkdir(s);
    cd(s);
    save_nii_oma(img,'image_POS.nii');
    cd ..
    files{1}=[pwd,filesep,s,filesep,'image_POS.nii'];
    cd ..
    
    cd NEG
    nii_p=load_nii(sprintf('delay_%.1f_zscoremats_tfce_corrp_tstat1.nii',d));
    nii_t=load_nii(sprintf('delay_%.1f_zscoremats_tstat1.nii',d));    
    img = nii_t.img.*(nii_p.img>0.99); 
    fprintf('voxels to show: %i\n',nnz(img>0));    
    mkdir(s);
    cd(s);
    save_nii_oma(img,'image_NEG.nii');
    cd ..
    files{2}=[pwd,filesep,s,filesep,'image_NEG.nii'];    
    cd ..
    
    mkdir(s);
    filename = create_mricroGL_script(s,files,[pwd,filesep,s]);
    
end

% delays = -12:2:12;
% for d=delays    
%     s = sprintf('delay%.1f',d);
%     cd(s);
%     image_autocrop('all');
%     cd ..            
% end