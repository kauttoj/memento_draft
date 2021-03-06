clc;
clear variables;
close all;

delays = -12:2:12;

for d=delays
    
    s = sprintf('delay%i',d);
    nii_p=load_nii(sprintf('delay_%.1f_corrmats_tfce_corrp_tstat1.nii',d));
    nii_t=load_nii(sprintf('delay_%.1f_corrmats_tstat1.nii',d));
    
    img = nii_t.img.*(nii_p.img>0.99); 
    fprintf('voxels to show: %i\n',nnz(img>0));
    
    mkdir(s);
    cd(s);
    save_nii_oma(img,'image.nii');
    cd ..
        
    filename = create_mricroGL_script(s,[pwd,filesep,s,filesep,'image.nii'],[pwd,filesep,s]);
    
end

delays = -12:2:12;
for d=delays
    
    s = sprintf('delay%i',d);
    cd(s);
    image_autocrop('all');
    cd ..
            
end