function ORIG_KRON_comparison_glm()

ORIG_S = {'S7','S8', 'S9', 'S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'}; % S13 incomplete
KRON_S = {'KRON_2','KRON_3','KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13','KRON_15','KRON_16'};

ORIG_data = load('ORIG/memento_LOOPER_results_ORIG.mat');
KRON_data = load('KRON/memento_LOOPER_results_KRON.mat');

assert(nnz(ORIG_data.delay_array - KRON_data.delay_array)==0);

delay_array = ORIG_data.delay_array;

clear ORIG_data KRON_data

N_orig = length(ORIG_S);
N_kron = length(KRON_S);

radius = 6;
RUN_PARAMETRIC = 0;

FWHM=5;
FSLDIR = 'fsl5.0-';
bramilapath = '/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila';

HOME = pwd;

%funAlign_createpool(7);

for k=1:length(delay_array)
    delay = delay_array(k);
    
    fprintf('Running delay %.2f (%i of %i)\n',delay_array(k),k,length(delay_array));
    
    cd(HOME);
    cd('KRON');
    s = sprintf('delay_%g_radius_%i',delay,radius);
    cd(s);
    nii = load_nii('PATTERN_mask.nii');
    KRON_mask{k}=nii.img;
    for i=1:N_kron
        outfile = [KRON_S{i},'_averaged_patterns_smooth.nii'];
        if ~exist(outfile,'file')
            infile = [pwd,filesep,KRON_S{i},'_averaged_patterns.nii'];
            img = simple_smooth(infile,outfile,FWHM,FSLDIR,bramilapath);
        else        
            nii=load_nii(outfile);
            img = nii.img;
        end
        KRON_data{k}{i}=mean(img,4);
    end
    
    cd(HOME);
    cd('ORIG');
    s = sprintf('delay_%g_radius_%i',delay,radius);
    cd(s);
    nii = load_nii('PATTERN_mask.nii');
    ORIG_mask{k}=nii.img;
    for i=1:N_orig
        outfile = [ORIG_S{i},'_averaged_patterns_smooth.nii'];
        if ~exist(outfile,'file')
            infile = [pwd,filesep,ORIG_S{i},'_averaged_patterns.nii'];
            img = simple_smooth(infile,outfile,FWHM,FSLDIR,bramilapath);
        else        
            nii=load_nii(outfile);
            img = nii.img;
        end
        ORIG_data{k}{i}=mean(img,4);
    end
end

cd(HOME);

if RUN_PARAMETRIC
    
    for i=1:length(delay_array)
        delay = delay_array(i);
        
        A = ORIG_data{i};
        B = KRON_data{i};
        
        n1=length(A);
        n2=length(B);
        
        vals1 = zeros(1,n1);
        vals2 = zeros(1,n2);
        
        mask1 = ORIG_mask{i};
        mask2 = KRON_mask{i};
        
        mask = mask1 & mask2;
        
        fprintf('Running delay %.2f (%i of %i)\n',delay_array(i),i,length(delay_array));
        fprintf('mask1 size %i, mask2 size %i, overlap %i\n',nnz(mask1),nnz(mask2),nnz(mask));
        
        ind = find(mask)';
        N_vox = length(ind);
        
        tmap = nan(size(mask));
        pmap = tmap;
        
        printind = round(linspace(N_vox/10,N_vox,10));
        kk=1;
        
        k=0;
        for j = ind
            k=k+1;
            if kk<11 && printind(kk)==k
                fprintf('...%i/10 (%i below p<0.001)\n',kk,nnz(pmap<0.001));
                kk=kk+1;
            end
            
            for n=1:n1
                vals1(n)=A{n}(j);
            end
            for n=1:n2
                vals2(n)=B{n}(j);
            end
            
            [h,p,~,b]=ttest2(vals1,vals2,'var','unequal');
            t = b.tstat;
            
            tmap(j)=t;
            pmap(j)=p;
            
        end
        
        all_tmaps{i}=tmap;
        all_pmaps{i}=pmap;
        
    end
    
    for i=1:length(delay_array)
        delay = delay_array(i);
        
        ind = find(~isnan(all_pmaps{i}));
        vals = all_pmaps{i}(ind);
        p_cor=mafdr(vals,'BHFDR',true);
        
        a = all_pmaps{i};
        a(ind)=p_cor;
        
        all_pmaps_cor{i} = a;
        
        res = zeros(size(a));
        res(a<0.05)=1;
        res = extentThreshold(res,30);
        
        res(res>0)=all_tmaps{i}(res>0);
        all_tmaps_cor{i} = res;
        
        if nnz(res<0)>10
            save_nii_oma(-res.*(res<0),sprintf('delay_%.1f_p005_ext30_UNIVARIATE_NEG.nii',delay_array(i)));
        end
        if nnz(res>0)>10
            save_nii_oma(res.*(res>0),sprintf('delay_%.1f_p005_ext30_UNIVARIATE_POS.nii',delay_array(i)));
        end
        
    end
    
end

%%

addpath('/m/nbe/scratch/braindata/kauttoj2/code/jannes_codes');
addpath('/triton/becs/scratch/braindata/kauttoj2/code/jannes_codes');

FSL_PREFIX = FSLDIR;

funAlign_createpool(7);

fprintf('Starting randomise loop\n');

N=length(delay_array);

resfiles = cell(1,N);
output_txt = cell(1,N);

parfor i=1:N
    
    fprintf('Running delay %.2f (%i of %i)\n',delay_array(i),i,N);
    
    A = ORIG_data{i}';
    B = KRON_data{i}';
    
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
    
    IMAGE = sprintf('delay_%.1f_zscoremats.nii',delay_array(i));
    save_nii_oma(C,IMAGE);
    
    mask1 = ORIG_mask{i};
    mask2 = KRON_mask{i};
    
    mask = mask1 & mask2;
    
    if nnz(mask)<130000
        error('BAD MASK!')
    end
    
    MASK = sprintf('delay_%.1f_mask.nii',delay_array(i));
    save_nii_oma(mask,MASK);
    
    DESIGN = zeros(size(C,4),2);
    DESIGN(1:n1,1)=1;
    DESIGN(n1+(1:n2),2)=1;
    
    CONTRAST = [-1 1];
    
    [resfiles{i},output_txt{i}] = randomise_ttest2([pwd,filesep,IMAGE],[pwd,filesep,MASK],5000,DESIGN,CONTRAST,FSL_PREFIX);
end

end

function img = simple_smooth(infile,outfile,FWHM,FSLDIR,bramilapath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fprintf('...smoothing data\n');
addpath(bramilapath);
cfg.bramilapath = bramilapath;
cfg.FSLDIR=FSLDIR;
cfg.infile = infile;
cfg.smooth_FWHM = FWHM;
cfg.smooth_method='SPM';
cfg.write=0;
cfg.infile_orig=cfg.infile;
cfg = bramila_smooth(cfg);
singleSubjectVols = cfg.vol;
save_nii_oma(singleSubjectVols,outfile);
img = singleSubjectVols;

end

