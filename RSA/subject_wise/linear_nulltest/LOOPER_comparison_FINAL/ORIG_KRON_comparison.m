clc;
clear variables;
close all;

ORIG_data = load('ORIG/memento_LOOPER_results_ORIG.mat');

KRON_data = load('KRON/memento_LOOPER_results_KRON.mat');

assert(nnz(ORIG_data.delay_array - KRON_data.delay_array)==0);

N1 = length(ORIG_data.all_cormaps);
N2 = length(KRON_data.all_cormaps);

assert(N1 == N2);

N=N1;

% for i=1:N
%     A = ORIG_data.all_cormaps{i}; 
%     B = KRON_data.all_cormaps{i};
%     
%     n1=length(A);
%     n2=length(B);
%     
%     vals1 = zeros(1,n1);
%     vals2 = zeros(1,n2);
%     
%     mask1 = ~isnan(A{1});
%     mask2 = ~isnan(B{1});  
%     
%     mask = mask1 & mask2;
%     
%     fprintf('Running delay %.2f (%i of %i)\n',ORIG_data.delay_array(i),i,N);
%     fprintf('mask1 size %i, mask2 size %i, overlap %i\n',nnz(mask1),nnz(mask2),nnz(mask));
%     
%     ind = find(mask)';
%     N_vox = length(ind);
%     
%     tmap = nan(size(mask));
%     pmap = tmap;
%     
%     printind = round(linspace(N_vox/10,N_vox,10));    
%     kk=1;
%         
%     k=0;    
%     for j = ind
%         k=k+1;
%         if kk<11 && printind(kk)==k
%             fprintf('...%i/10 (%i below p<0.001)\n',kk,nnz(pmap<0.001));
%             kk=kk+1;
%         end
%         
%         for n=1:n1
%            vals1(n)=A{n}(j);
%         end
%         for n=1:n2
%            vals2(n)=B{n}(j);
%         end        
%         
%         [h,p,~,b]=ttest2(vals1,vals2,'var','unequal');
%         t = b.tstat;
%         
%         tmap(j)=t;
%         pmap(j)=p;
%         
%     end    
%     
%     all_tmaps{i}=tmap;
%     all_pmaps{i}=pmap;
%     
% end
% 
% for i=1:N
%     
%     ind = find(~isnan(all_pmaps{i}));
%     vals = all_pmaps{i}(ind);
%     p_cor=mafdr(vals,'BHFDR',true);
%     
%     a = all_pmaps{i};
%     a(ind)=p_cor;
%     
%     all_pmaps_cor{i} = a;
%     
%     res = zeros(size(a));
%     res(a<0.05)=1;
%     res = extentThreshold(res,30);            
%     
%     res(res>0)=all_tmaps{i}(res>0);
%     all_tmaps_cor{i} = res;
%     
%     if nnz(res<0)>10
%         save_nii_oma(-res.*(res<0),sprintf('delay_%.1f_p005_ext30_NEG.nii',ORIG_data.delay_array(i)));
%     end
%     if nnz(res>0)>10
%         save_nii_oma(res.*(res>0),sprintf('delay_%.1f_p005_ext30_POS.nii',ORIG_data.delay_array(i)));
%     end
%     
% end

%%

addpath('/m/nbe/scratch/braindata/kauttoj2/code/jannes_codes');
addpath('/triton/becs/scratch/braindata/kauttoj2/code/jannes_codes');

FSL_PREFIX = 'fsl5.0-';

funAlign_createpool(7);

fprintf('Starting randomise loop\n');

resfiles = cell(1,N);
output_txt = cell(1,N);

parfor i=1:N
    
    fprintf('Running delay %.2f (%i of %i)\n',ORIG_data.delay_array(i),i,N);
    
    A = ORIG_data.all_cormaps{i}; 
    B = KRON_data.all_cormaps{i};
    
    n1=length(A);
    n2=length(B);
    
    apu = [A;B];
    C=apu{1};
    for j=2:length(apu)        
       C = cat(4,C,apu{j});
    end
    for j=1:length(apu)
        assert(nnz(abs(squeeze(C(:,:,:,j))-apu{j})>eps)==0);
    end
    apu = 0;
    
    IMAGE = sprintf('delay_%.1f_corrmats.nii',ORIG_data.delay_array(i));
    save_nii_oma(C,IMAGE);
        
    mask = ~isnan(sum(C,4));
    
    if nnz(mask)<130000
        error('BAD MASK!')
    end
    
    MASK = sprintf('delay_%.1f_mask.nii',ORIG_data.delay_array(i));
    save_nii_oma(mask,MASK);
    
    DESIGN = zeros(size(C,4),2);
    DESIGN(1:n1,1)=1;
    DESIGN(n1+(1:n2),2)=1;
    
    CONTRAST = [1 -1];
    
    [resfiles{i},output_txt{i}] = randomise_ttest2([pwd,filesep,IMAGE],[pwd,filesep,MASK],5000,DESIGN,CONTRAST,FSL_PREFIX);
end


