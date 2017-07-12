clc;
close all;
clear variables;


A1 = load('KRON/memento_LOOPER_results_KRON.mat');
A2 = load('KRON/ADDITION/memento_LOOPER_results_KRON.mat');

KRON_corrmaps = [A2.all_cormaps,A1.all_cormaps];
KRON_delays = [A2.delay_array,A1.delay_array];

A1 = load('ORIG/memento_LOOPER_results_ORIG.mat');
A2 = load('ORIG/ADDITION/memento_LOOPER_results_ORIG.mat');

ORIG_corrmaps = [A2.all_cormaps,A1.all_cormaps];
ORIG_delays = [A2.delay_array,A1.delay_array];

assert(all(ORIG_delays-KRON_delays==0));
delays = KRON_delays;

nii=load_nii('ORIG_KRON_mask.nii');
mask = nii.img;

load('/m/nbe/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila/external/HO_2mm_rois.mat');

N_roi = length(rois);
T = length(delays);

N_ORIG = length(ORIG_corrmaps{1});
N_KRON = length(KRON_corrmaps{1});

ORIG_ts = nan(N_ORIG,N_roi,T);
KRON_ts = nan(N_KRON,N_roi,T);

mask = ones(size(mask));
for t=1:T
    for j=1:N_ORIG
        mask = mask.*~isnan(ORIG_corrmaps{t}{j});
    end
    for j=1:N_KRON
        mask = mask.*~isnan(KRON_corrmaps{t}{j});
    end
end

fprintf('mask size %i\n',nnz(mask));

bad = [];
for n=1:N_roi
    ind = sub2ind(size(mask),rois(n).map(:,1),rois(n).map(:,2),rois(n).map(:,3));
    ind(~mask(ind))=[];
    
    koko(n)=length(ind);
    
    if length(ind)<10
        bad(end+1)=n;
    else
        
        for t=1:T
            for j=1:N_ORIG
                vals = ORIG_corrmaps{t}{j}(ind);
                if nnz(isnan(vals))>0
                    error('!!')
                end
                ORIG_ts(j,n,t)=mean(vals);
            end
            for j=1:N_KRON
                vals = KRON_corrmaps{t}{j}(ind);
                if nnz(isnan(vals))>0
                    error('!!')
                end                
                KRON_ts(j,n,t)=mean(vals);
            end
        end
    end
end

ORIG_ts(:,bad,:)=[];
KRON_ts(:,bad,:)=[];
rois(bad)=[];
koko(bad)=[];
N_roi = length(rois);

assert(nnz(isnan(ORIG_ts))==0 && nnz(isnan(KRON_ts))==0);

pvals1=nan(N_roi,T);
tvals1=nan(N_roi,T);
pvals11=nan(N_roi,T);
tvals11=nan(N_roi,T);
pvals2=nan(N_roi,T);
tvals2=nan(N_roi,T);
mean_ORIG_ts=nan(N_roi,T);
mean_KRON_ts=nan(N_roi,T);
for t=1:T
    [~,p,~,s]=ttest2(ORIG_ts(:,:,t),KRON_ts(:,:,t),'tail','both','var','unequal');
    pvals2(:,t)=p;
    tvals2(:,t)=s.tstat;    
    
    mean_ORIG_ts(:,t) = mean(ORIG_ts(:,:,t));
    mean_KRON_ts(:,t) = mean(KRON_ts(:,:,t));
        
    [~,p,~,s]=ttest(squeeze(ORIG_ts(:,:,t)),0,'tail','both');
    pvals1(:,t)=p;
    tvals1(:,t)=s.tstat;
    
    [~,p,~,s]=ttest(squeeze(KRON_ts(:,:,t)),0,'tail','both');
    pvals11(:,t)=p;
    tvals11(:,t)=s.tstat;    
end

all_tvals1=tvals1;
all_pvals1=pvals1;

all_tvals2=tvals2;
all_pvals2=pvals2;

all_tvals11=tvals11;
all_pvals11=pvals11;


%% ORIG all different and nonzero

ind = delays>=-6 & delays<=+4;
pvals1 = all_pvals1(:,ind);
tvals1 = all_tvals1(:,ind);
[pvals_cor1,qvals1] = mafdr(pvals1(:));
pvals_cor1 = reshape(pvals_cor1,size(pvals1));
pvals2 = all_pvals2(:,ind);
tvals2 = all_tvals2(:,ind);
[pvals_cor2,qvals2] = mafdr(pvals2(:));
pvals_cor2 = reshape(pvals_cor2,size(pvals2));

SELECTED_rois =find(sum(pvals_cor2<0.01,2)>0);
SELECTED_rois1 =find(sum(pvals_cor1<0.01,2)>0 & sum(pvals_cor2<0.01,2)>0);
if ~isempty(setdiff(SELECTED_rois,SELECTED_rois1))
    error('!!!!')
end

ind1=find(ind);
ind = find(delays<=6 & delays>=-30);

pval_table = all_pvals2(SELECTED_rois,ind);%zeros(length(SELECTED_rois),lenght(ind1));
max_table = 0*pval_table;

leg = [];
lat=[];
maxtime=[];
k=0;
for i=SELECTED_rois'
   k=k+1;
      
   [~,j]=max(abs(mean_ORIG_ts(i,ind)));
   
   val = mean_ORIG_ts(i,ind(j));
   if val<0
       error('!!')
   end
   
   maxtime(k)=delays(ind(j));
   
   max_table(k,j)=1;
   
   leg{k}=rois(i).label;
   leg{k}=simplify_roitext(leg{k});
   
   leg{k}=[leg{k},' (',num2str(koko(i)),')'];
   
   if strfind(leg{k},'R')
      lat(k)=1;
   else
      lat(k)=0; 
   end
   
   fprintf('%i: peak location %+.2fs\n',k,maxtime(k));
end

data = mean_ORIG_ts(SELECTED_rois,ind);
[~,i]=sort(maxtime,'descend');
data = data./repmat(max(data')',[1,size(data,2)]);
handle = plot_matrix(data(i,:),'',delays(ind),leg(i));
%plug_annotations(handle,max_table(i,:),'\makebox[5][5]',32);
table = (pval_table(i,:)<0.001);
table1 = 0*table;
j=ismember(ind,ind1);
table1(:,find(j)) = pvals_cor2(SELECTED_rois,:)<0.01;
plug_annotations(handle,table.*(table1(i,:)==0),'\textasteriskcentered',14);
plug_annotations(handle,table1(i,:),'\textasteriskcentered\textasteriskcentered',10);

hold on;
plot(find(delays(ind)==-6)*[1,1]-0.5,[0,size(data,1)+1],'k');
plot(find(delays(ind)==6)*[1,1]-0.5,[0,size(data,1)+1],'k');
plot(find(delays(ind)==0)*[1,1],[0,size(data,1)+1],'k--');
xlabel('Onset delay [s]')

annotation('textbox',[.84 .7 .1 .1],'String','\textasteriskcentered=$p<0.001$','FitBoxToText','on','interpreter','latex','LineStyle','none','FontSize',11);
annotation('textbox',[.84 .76 .1 .1],'String','\textasteriskcentered\textasteriskcentered=$p<0.01$ (FDR)','FitBoxToText','on','interpreter','latex','LineStyle','none','FontSize',11);
annotation('textbox',[.0 .0 .1 .1],'String','','FitBoxToText','on','interpreter','latex','LineStyle','none','FontSize',11);

close all;

figure('Position',[675   536   839   548]);
i=find(maxtime>=0 & lat>0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

figure('Position',[675   536   839   548]);
i=find(maxtime>=0 & lat==0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

figure('Position',[675   536   839   548]);
i=find(maxtime<0 & lat>0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

figure('Position',[675   536   839   548]);
i=find(maxtime<0 & lat==0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

%% ORIG nonzero

ind = delays>=-6 & delays<=+4;
pvals1 = all_pvals1(:,ind);
tvals1 = all_tvals1(:,ind);
[pvals_cor1,qvals1] = mafdr(pvals1(:));
pvals_cor1 = reshape(pvals_cor1,size(pvals1));
pvals2 = all_pvals2(:,ind);
tvals2 = all_tvals2(:,ind);
[pvals_cor2,qvals2] = mafdr(pvals2(:));
pvals_cor2 = reshape(pvals_cor2,size(pvals2));
SELECTED_rois =find(sum(pvals_cor1<0.01,2)>0);

ind = find(delays<7 & delays>-21);

leg = [];
lat=[];
maxtime=[];
k=0;
for i=SELECTED_rois'
   k=k+1;
      
   [~,j]=max(mean_ORIG_ts(i,ind));
   
   maxtime(k)=delays(ind(j));
   
   leg{k}=rois(i).label;
   leg{k}=simplify_roitext(leg{k});
   
   leg{k}=[leg{k},' (',num2str(maxtime(k)),'s)'];
   
   if strfind(leg{k},'R')
      lat(k)=1;
   else
      lat(k)=0; 
   end
   
   fprintf('%i: peak location %+.2fs\n',k,maxtime(k));
end

data = mean_ORIG_ts(SELECTED_rois,ind);
[~,i]=sort(maxtime,'descend');
data = data./repmat(max(data')',[1,size(data,2)]);
plot_matrix(data(i,:),'',delays(ind),leg(i))
hold on;
plot(find(delays(ind)==-6)*[1,1]-0.5,[0,size(data,1)+1],'k');
plot(find(delays(ind)==6)*[1,1]-0.5,[0,size(data,1)+1],'k');
plot(find(delays(ind)==0)*[1,1],[0,size(data,1)+1],'k--');
xlabel('Onset delay [s]','FontSize',12)

close all;

figure('Position',[675   536   839   548]);
i=find(maxtime>=0 & lat>0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

figure('Position',[675   536   839   548]);
i=find(maxtime>=0 & lat==0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

figure('Position',[675   536   839   548]);
i=find(maxtime<0 & lat>0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

figure('Position',[675   536   839   548]);
i=find(maxtime<0 & lat==0);
%subplot(2,1,1);
plot(delays(ind),mean_ORIG_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;

%% KRON nonzero


ind = delays>=-6 & delays<=+4;
pvals11 = all_pvals11(:,ind);
tvals11 = all_tvals11(:,ind);
[pvals_cor11,qvals11] = mafdr(pvals11(:));
pvals_cor11 = reshape(pvals_cor11,size(pvals11));
SELECTED_rois =find(sum(pvals_cor11<0.01,2)>0);

ind = find(delays<7 & delays>-21);

leg = [];
lat=[];
maxtime=[];
k=0;
for i=SELECTED_rois'
   k=k+1;
      
   [~,j]=max(mean_KRON_ts(i,ind));
   
   maxtime(k)=delays(ind(j));
   
   leg{k}=rois(i).label;
   leg{k}=simplify_roitext(leg{k});
   
   leg{k}=[leg{k},' (',num2str(maxtime(k)),'s)'];
   
   if strfind(leg{k},'R')
      lat(k)=1;
   else
      lat(k)=0; 
   end
   
   fprintf('%i: peak location %+.2fs\n',k,maxtime(k));
end

data = mean_KRON_ts(SELECTED_rois,ind);
[~,i]=sort(maxtime,'descend');
data = data./repmat(max(data')',[1,size(data,2)]);
plot_matrix(data(i,:),'',delays(ind),leg(i))
hold on;
plot(find(delays(ind)==-6)*[1,1]-0.5,[0,size(data,1)+1],'k');
plot(find(delays(ind)==6)*[1,1]-0.5,[0,size(data,1)+1],'k');
plot(find(delays(ind)==0)*[1,1],[0,size(data,1)+1],'k--');
xlabel('Onset delay [s]','FontSize',12)

close all;

figure('Position',[675   536   839   548]);
i=find(maxtime>=0 & lat>0);
%subplot(2,1,1);
try
plot(delays(ind),mean_KRON_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;
catch
    
end

figure('Position',[675   536   839   548]);
i=find(maxtime>=0 & lat==0);
try
%subplot(2,1,1);
plot(delays(ind),mean_KRON_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;
catch
    
end

figure('Position',[675   536   839   548]);
i=find(maxtime<0 & lat>0);
try
%subplot(2,1,1);
plot(delays(ind),mean_KRON_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;
catch
    
end

figure('Position',[675   536   839   548]);
i=find(maxtime<0 & lat==0);
try
%subplot(2,1,1);
plot(delays(ind),mean_KRON_ts(SELECTED_rois(i),ind),'LineWidth',2);
set(gca,'FontSize',14)
xlabel('Onset delay [s]','FontSize',24,'Interpreter','Latex');
ylabel('$\hat{r}$','Interpreter','latex','FontSize',24);
a = legend(leg(i),'location','northwest','FontSize',16);
legend('boxoff')
omagrid;
catch
    
end