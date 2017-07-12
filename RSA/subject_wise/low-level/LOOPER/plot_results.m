clc;
clear all;
close all;

addpath('/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA');

HOME = pwd;

delay_array = -20:2.5:25;

for i=1:length(delay_array)
   
    fprintf('----- Delay = %g\n',delay_array(i));
        
    str = sprintf('duration_4_delay_%g_radius_5',delay_array(i));
    cd(str);
    load memento_rsa_results;
        
    fprintf('Subject models sizes:\n');
    for k=1:length(all_models)
         fprintf('...%i: %i\n',k,size(all_models{k}.RDM,1));
    end  
    
    map = mean_map;
    map(isnan(map))=0;
    map(map<Th_kauppi(6))=0;    
    allmaps_size_k1(i)=nnz(map);  
    map = extentThreshold(map,20);
    allmaps{i}=map;
    allmaps_size(i)=nnz(map);
    
    map = mean_map;
    map(isnan(map))=0;
    map(map<Th_kauppi(9))=0;
    map = extentThreshold(map,50);
    allmaps_k50{i}=map;
    allmaps_size_k50(i)=nnz(map);
    
    fprintf('FDR thresholds are = %f (p<0.05), %f (p<0.01)\n',Th_kauppi(2),Th_kauppi(6));
    
    fprintf('Final maps size is %i\n',allmaps_size_k1(i));
    
    [glmvals{i},glmpvals{i},glmpvals_corrected{i},resultmap,glm_mask]=glm_analysis('averaged_patterns_summary.mat',pwd,'EPI_mask.nii','right');
    
    map=0*glm_mask;
    ind1 = find(glm_mask);
    ind2 = glmpvals_corrected{i}<0.01 & glmvals{i}>0;    
    map(ind1(ind2))=glmvals{i}(ind2); 
    map = extentThreshold(map,20); 
    glmallmaps{i}=map;
    glmallmaps_size(i)=nnz(map);

    map=0*glm_mask;
    ind1 = find(glm_mask);
    ind2 = glmpvals{i}<0.001 & glmvals{i}>0;
    map(ind1(ind2))=glmvals{i}(ind2);
    map = extentThreshold(map,50); 
    glmallmaps_k50{i}=map;
    glmallmaps_size_k50(i)=nnz(map);
    
    
    cd(HOME);       
    
end

img = zeros(91,109,91,length(allmaps));
cormap = zeros(length(allmaps),length(allmaps));
overlap = cormap;

glmcormap = cormap;
glmoverlap = overlap;

for i=1:length(allmaps)
    img(:,:,:,i)=allmaps{i};
    if nnz(allmaps{i})>0
        save_nii_oma(allmaps{i},sprintf('meanmap_p001fdr_k20_delay%g.nii',delay_array(i)));
    end
    for ii=(i+1):length(allmaps)
       %a=allmaps_k50{i}(supermask>0);
       %b=allmaps_k50{ii}(supermask>0);
       a=allmaps{i}(supermask>0);
       b=allmaps{ii}(supermask>0);       
       
       if nnz(a)>0 && nnz(b)>0
         cormap(i,ii)=corr2(a(:),b(:));
         overlap(i,ii)=2*nnz((allmaps{i}>0).*(allmaps{ii}>0))/(nnz(allmaps{i}>0) + nnz(allmaps{ii}>0));
       end
       
       a=glmallmaps_k50{i}(glm_mask>0);
       b=glmallmaps_k50{ii}(glm_mask>0);
       if nnz(a)>0 && nnz(b)>0
         glmcormap(i,ii)=corr2(a(:),b(:));
         glmoverlap(i,ii)=2*nnz((glmallmaps_k50{i}>0).*(glmallmaps_k50{ii}>0))/(nnz(glmallmaps_k50{i}>0) + nnz(glmallmaps_k50{ii}>0));
       end       
       
    end
end
cormap=cormap+cormap';
overlap=overlap+overlap';

glmcormap=glmcormap+glmcormap';
glmoverlap=glmoverlap+glmoverlap';

a=sum(cormap)>0;
cut = max(1,(find(a,1,'first')-1)):min(length(a),(find(a,1,'last')+1));

delay_array_cut = delay_array(cut);
cormap_cut=cormap(cut,cut);
overlap_cut=overlap(cut,cut);

plot_matrix_simple(overlap,'test',delay_array_cut,delay_array_cut)

save_nii_oma(img,'meanmaps_p001fdr_k20.nii');

save('postprosessed_results.mat','-v7.3');

%%

% [a1,a2,a3]=plotyy(delay_array,glmallmaps_size_k50,delay_array,allmaps_size_k50);
% a=get(a1,'Ylabel');
% set(a{1},'string','GLM','FontSize',16)
% set(a{2},'string','MVPA','FontSize',16)
% xlabel('Delay [s]','FontSize',16)
% set(a1(1),'FontSize',14)
% set(a1(2),'FontSize',14)

plot(delay_array,[glmallmaps_size_k50;allmaps_size_k50]/1000);
axis tight;
box on;
ylabel('#voxels (\times 10^4)','FontSize',16);
xlabel('Onset delay [s]','FontSize',16);
set(gca,'FontSize',14);
legend('GLM','MVPA')
hold on;
plot([0,0],[0,max(glmallmaps_size_k50+50)/1000],'r')

