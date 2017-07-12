
load('memento_kron_addition.mat');

clear kron_addition_response kron_addition_response_srt;
for i=1:size(memento_kron_addition,1)
    k=0;
    for j=4:3:size(memento_kron_addition,2)
        k=k+1;
        kron_addition_response(i,k)=0;
        kron_addition_response_srt{i,k}=memento_kron_addition{i,j};        
        if ~isempty(memento_kron_addition{i,j})
            if findstr('Muistan yllä olevan kuvan toistuneen',memento_kron_addition{i,j})
                kron_addition_response(i,k)=1;
            else
               error('!!!!');
            end
        end
    end
end
fprintf('Success!\n')