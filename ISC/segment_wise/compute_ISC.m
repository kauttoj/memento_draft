function [isc_matrix,element_indices,ISC_vals1,ISC_vals2,permdist1,permdist2] = compute_ISC(DATA1,DATA2,iter,nworker)
%COMPUTE_ISC Summary of this function goes here
%   Detailed explanation goes here

if length(size(DATA1))==2 && length(size(DATA2))==2
    a(:,1,:)=DATA1;
    DATA1=a;
    b(:,1,:)=DATA2;
    DATA2=b;
    clear a b;
end

N_sub1=size(DATA1,3);
N_sub2=size(DATA2,3);
N_sub = N_sub1+N_sub2;

T = size(DATA1,1);
N_voxel = size(DATA1,2);

if T~=size(DATA2,1) || N_voxel~=size(DATA2,2)
    error('Data size does not match!')
end

mat_ind = find(triu(ones(N_sub,N_sub),1));
isc_matrix = nan(length(mat_ind),N_voxel);

mat = zeros(N_sub,N_sub);
mat(1:N_sub1,1:N_sub1)=1;
mat(N_sub1+(1:N_sub2),N_sub1+(1:N_sub2))=2;
mat((1:N_sub1),N_sub1+(1:N_sub2))=3;
element_indices = mat(mat_ind);

element_indices1 = mat_ind(element_indices==1);
element_indices2 = mat_ind(element_indices==2);

ISC_vals1 = nan(1,N_voxel);
ISC_vals2 = nan(1,N_voxel);

print_ind = round(linspace(N_voxel/20,N_voxel-N_voxel/20,20));
k=1;

fprintf('\nStarting ISC computing\n');
for voxel = 1:N_voxel
    
    data1 = squeeze(DATA1(:,voxel,:));
    data2 = squeeze(DATA2(:,voxel,:));
    
    data = [data1,data2];
    mat = corr(data);
    
    isc_matrix(:,voxel) = mat(mat_ind);
    
    ISC_vals1(voxel) = mean(mat(element_indices1));
    ISC_vals2(voxel) = mean(mat(element_indices2));
    
    if k<21 && voxel==print_ind(k)
        fprintf('... %i/20 done\n',k);
        k=k+1;
    end
    
end

if nargin==4 && nargout==6 && iter>100
    
    voxel = randsample(N_voxel,iter,true);
    permdist1=nan(1,iter);
    permdist2=nan(1,iter);
    mat_ind1 = (triu(ones(N_sub1,N_sub1),1))>0;
    mat_ind2 = (triu(ones(N_sub2,N_sub2),1))>0;
    %create local cluster
    try
        myCluster = gcp('nocreate');
        if isempty(myCluster)
            %delete(gcp)
            myCluster = parcluster('local');
            if nargin>4 && ~isempty(nworker)
                myCluster.NumWorkers=nworker;
            end
            parpool(myCluster);
        end
        N_workers = myCluster.NumWorkers;
    catch err % old matlab?
        if ~matlabpool('size')
            if nargin>4 && ~isempty(nworker)
                eval(['matlabpool local ',num2str(nworker)]);
            else
                matlabpool local
            end
        end
        N_workers = matlabpool('size');
    end
    
    
    fprintf('\nStarting mean ISC permutations (total %i)\n',iter);
    
    try
        tic;
        parfor i = 1:iter
            
            data1 = squeeze(DATA1(:,voxel(i),:));
            data2 = squeeze(DATA2(:,voxel(i),:));
            
            shift1 = ceil(T*rand(1,N_sub1));
            shift2 = ceil(T*rand(1,N_sub2));
            
            for j=1:N_sub1
                data1(:,j)=data1([shift1(j):end,1:(shift1(j)-1)],j);
            end
            for j=1:N_sub2
                data2(:,j)=data2([shift2(j):end,1:(shift2(j)-1)],j);
            end
            
            mat1 = corr(data1);
            mat2 = corr(data2);
            
            permdist1(i) = mean(mat1(mat_ind1));
            permdist2(i) = mean(mat2(mat_ind2));
            
        end
    catch err
        warning('Failed to run parallel loop! Returning to serial...');
        tic;
        for i = 1:iter
            
            data1 = squeeze(DATA1(:,voxel(i),:));
            data2 = squeeze(DATA2(:,voxel(i),:));
            
            shift1 = ceil(T*rand(1,N_sub1));
            shift2 = ceil(T*rand(1,N_sub2));
            
            for j=1:N_sub1
                data1(:,j)=data1([shift1(j):end,1:(shift1(j)-1)],j);
            end
            for j=1:N_sub2
                data2(:,j)=data2([shift2(j):end,1:(shift2(j)-1)],j);
            end
            
            mat1 = corr(data1);
            mat2 = corr(data2);
            
            permdist1(i) = mean(mat1(mat_ind1));
            permdist2(i) = mean(mat2(mat_ind2));
            
        end
    end
    
    a=toc;
    fprintf('Finished (%.1fs)!\n',a);
    
end

fprintf('all done!\n\n');
