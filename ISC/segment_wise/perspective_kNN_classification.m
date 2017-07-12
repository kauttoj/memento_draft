function [accuracy_out,results_out,pval_out,nulldist] = perspective_kNN_classification(ISC_matrices,KK,groups,iterations)

%Simple kNN-classification for ISC matrices.
%--------------------------------------------------------------------------
% ISC_matrices = subject_pair x voxel
%
% if ~iscell(ISC_matrices)
%     a{1}=ISC_matrices;
%     ISC_matrices=a;
% end

x=size(ISC_matrices,1);
N=(.5+sqrt(.5^2+x*2));

% if length(ISC_matrices)>1
%     a=0;
%     for s=1:length(ISC_matrices)
%         a = a + ISC_matrices{s};
%     end
%     ISC_matrices{s+1}=a/length(ISC_matrices);
% end

fprintf('\n----- Starting classification analysis -----\n');

%for s=1:length(ISC_matrices)
    
    %fprintf('Starting classification (set %i of %i)\n',s,length(ISC_matrices));
    fprintf('...Computing kNN classification accuracy\n');    
    corMat = ISC_matrices;
    results=zeros(size(corMat,2),N);
    for sub=1:N
        subs=1:N;
        subs(sub)=[];
        grps=groups;
        grps(sub)=[];
        mask=false(N);
        mask(sub,:)=true;
        mask(:,sub)=true;
        
        ind = mask(triu(true(N),1));
        [~,idx]=sort(corMat(ind,:),'descend');
        for K=KK
            results(:,sub)=results(:,sub) + mode(grps(idx(1:K,:))==groups(sub))';
        end
               
    end
    results = results/length(KK);
    
    fprintf('...Computing null-distribution\n');
    voxels = randsample(size(corMat,2),iterations,true);
    nulldist = nan(1,iterations);
    indmat = triu(true(N),1);
    parfor i=1:iterations
        arr = zeros(1,N);
        groups_mixed = groups(randperm(N));
        for sub=1:N
            subs=1:N;
            subs(sub)=[];
            grps=groups_mixed;
            grps(sub)=[];
            mask=false(N);
            mask(sub,:)=true;
            mask(:,sub)=true;
            ind = mask(indmat);
            [~,idx]=sort(corMat(ind,voxels(i)),'descend');
            for K=KK
                arr(sub)=arr(sub) + mode(grps(idx(1:K))==groups_mixed(sub));            
            end
            
        end
        arr=arr/length(KK);
        nulldist(i)=mean(arr);
    end
    
    results_out=results;
    accuracy_out=sum(results,2)/N;
    pval=0*accuracy_out;
    fprintf('...Computing p-values ');
    for i=1:length(accuracy_out)
        pval(i)=1 - (nnz(accuracy_out(i)>nulldist)/iterations);
    end
    pval_out=pval;
    fprintf('(%i of %i with p<0.01)\n',nnz(pval<0.01),length(pval));
    fprintf('...all done!\n\n');
%end
