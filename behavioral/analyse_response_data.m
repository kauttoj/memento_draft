clc;
clear all;
close all;

% Pian tekemää excel tiedostoa muutettu siten, että kaikki turhat rivit ja sarakkeet on
% poistettu. Lisäksi sarakkeiden järjestys on muutettu siten, että
% importance ja eureka kysymykset tulee samassa järjestyksessä (aina ensin
% importance). Lisäksi taulukko tallennettu officella (muuten ei toimi matlabissa).

% HUOM! SARAKEINDEKSIT ON HARD-KOODATTU TÄMÄN FILEN MUKAISESTI!
FILENAME = 'Memento jälkikysely ALL (Responses) 15June2015_pt_EDITED_win.xlsx';
% mitkä kuvista on key-frameja
KEYFRAME_IMAGE_IDs = 1:15;
REMOVE_IDs = [];

CHOSEN_SUBJECTS = {'MEM_5','MEM_7','MEM_8','MEM_9','MEM_10','MEM_12','MEM_13','MEM_14','MEM_15','MEM_16','MEM_17','MEM_19','MEM_21','MEM_22','MEM_23'};
%CHOSEN_SUBJECTS = {'KRON_2','KRON_3','KRON_4','KRON_5','KRON_6','KRON_7','KRON_8','KRON_9','KRON_10','KRON_12','KRON_13','KRON_15','KRON_16'};

[num,txt,raw] = xlsread(FILENAME);
last_col = length(num(1,:));

ROWS = max(size(num,1),size(txt,1));
COLS = max(size(num,2),size(txt,2));

raw = raw(1:ROWS,1:COLS);

fprintf('Table has %i rows and %i columns\n',size(txt,1),size(txt,2));

k=0;
while 1
    if isnan(num(1,end-k))
       k=k+1;
    else
       break; 
    end
end

subject_ID = txt(3:end,2);
image_desc = txt(2,5:(last_col-k));
image_ID = num(1,5:(last_col-k));
responses = txt(3:end,5:(last_col-k));
gender = txt(3:end,4);
birthyear = num(3:end,3);

% for i=1:length(subject_ID)
%     if strcmp(subject_ID{i},'MEM6')
%         subject_ID{i} = 'MEM_6';
%     elseif strcmp(subject_ID{i},'MEM7')
%         subject_ID{i} = 'MEM_7';            
%     end
% end
add_kron_15=0;
add_kron_16=0;

pois=[];
replace_table{1} = {'mem_14_1','MEM_14'};
for i=1:length(subject_ID)
    for j=1:length(replace_table)
        if strcmp(subject_ID{i},replace_table{j}{1,1})
            fprintf('Replaced %s with %s\n',subject_ID{i},replace_table{j}{1,2});
            subject_ID{i}=replace_table{j}{1,2};            
        end
    end        
    

    
end

for i=1:length(CHOSEN_SUBJECTS)
    if strcmp(CHOSEN_SUBJECTS{i},'KRON_15')
        add_kron_15 = 1;
        pois(end+1)=i;
    end    
    if strcmp(CHOSEN_SUBJECTS{i},'KRON_16')
        add_kron_16 = 1;
        pois(end+1)=i;
    end    
end
CHOSEN_SUBJECTS(pois)=[];

s = convert_str2num(subject_ID,CHOSEN_SUBJECTS);
ind = find(~isnan(s));
if length(ind)~=length(CHOSEN_SUBJECTS)
   error('Not all subject IDs were found (check data)!') 
end
subject_ID=subject_ID(ind);
responses=responses(ind,:);

% kuvan toistumiset:
repeat_str = {'','Muistan yllä olevan kuvan toistuneen'};
repeat_responses = responses(:,1:30);
repeat_image_IDs=image_ID(1:30);
if length(unique(repeat_image_IDs))<length(repeat_image_IDs)
    error('Repeating IDs!')
end
[~,ind]=sort(repeat_image_IDs,'ascend');
repeat_responses=repeat_responses(:,ind);
repeat_image_IDs=repeat_image_IDs(ind);
pois = ismember(repeat_image_IDs,REMOVE_IDs);
repeat_image_IDs(pois)=[];
repeat_responses(:,pois)=[];

% kuvan tärkeys
importance_str={'En osaa sanoa','Ei erityisen tärkeä','Oli tärkeä'};
importance_responses=responses(:,31:2:end);
importance_image_IDs=image_ID(31:2:end);
if length(unique(importance_image_IDs))<length(importance_image_IDs)
    error('Repeating IDs!')
end
[~,ind]=sort(importance_image_IDs,'ascend');
importance_responses=importance_responses(:,ind);
importance_image_IDs=importance_image_IDs(ind);
pois = ismember(importance_image_IDs,REMOVE_IDs);
importance_image_IDs(pois)=[];
importance_responses(:,pois)=[];

% oivallukset:
eureka_str = {'En osaa sanoa','Ei erityistä oivallusta','Koin merkittävän oivalluksen'};
eureka_responses=responses(:,32:2:end);
eureka_image_IDs=image_ID(32:2:end);
if length(unique(eureka_image_IDs))<length(eureka_image_IDs)
    error('Repeating IDs!')
end
[~,ind]=sort(eureka_image_IDs,'ascend');
eureka_responses=eureka_responses(:,ind);
eureka_image_IDs=eureka_image_IDs(ind);
pois = ismember(eureka_image_IDs,REMOVE_IDs);
eureka_image_IDs(pois)=[];
eureka_responses(:,pois)=[];

%----------------------------------------

if nnz(sort(repeat_image_IDs)-sort(importance_image_IDs))>0
    error('!!!')
end
if nnz(sort(eureka_image_IDs)-sort(importance_image_IDs))>0
    error('!!!')
end

if add_kron_15+add_kron_16>0
    read_kron_addition;
    if add_kron_15>0
        repeat_responses = [repeat_responses;kron_addition_response_srt(1,:)];
        subject_ID{end+1}='KRON_15';
    end
    if add_kron_16>0
        repeat_responses = [repeat_responses;kron_addition_response_srt(2,:)];
        subject_ID{end+1}='KRON_16';
    end
end

model = 0*repeat_image_IDs;
model(ismember(repeat_image_IDs,KEYFRAME_IMAGE_IDs))=2;
model(~ismember(repeat_image_IDs,KEYFRAME_IMAGE_IDs))=1;
repeat_responses_num = nan(size(repeat_responses));
for i=1:length(repeat_image_IDs)
    s = convert_str2num(repeat_responses(:,i),repeat_str);
    s(isnan(s))=[];    
    repeat_responses_num(:,i)=s;
    repeat_pval(i)=myBinomTest(nnz(s==model(i)),length(s),0.5,'right');    
end
for i=1:size(repeat_responses_num,1)
    repeat_corr_count(i)=nnz(repeat_responses_num(i,:)==model);   
    repeat_corr_count_pval_binom(i)=myBinomTest(nnz(repeat_responses_num(i,:)==model),length(model),0.5,'right'); 
%     nulvals=nan(1,20000);
%     L=length(model);
%     for j=1:length(nulvals)
%         nullmodel = model(randperm(L));
%         nulvals(j)=nnz(repeat_responses_num(i,:)==nullmodel); 
%     end
%     repeat_corr_count_pval(i)=nnz(nulvals>=repeat_corr_count(i))/length(nulvals);    
end
clear nulvals nullmodel

ind = find(ismember(repeat_image_IDs,KEYFRAME_IMAGE_IDs));
repeat_pval_keyframe = repeat_pval(ind);
%repeat_pval_corrected = mafdr(repeat_pval(ind),'BHFDR',true);
fprintf('\n---KEYFRAMES---\n');
fprintf('Correctly identified (p<0.05) frames: %s (%i of %i)\n',num2str(repeat_image_IDs(repeat_pval<0.05)),nnz(repeat_pval<0.05),length(repeat_pval));
fprintf('Correctly identified (p<0.05) keyframes: %s (%i of %i)\n',num2str(repeat_image_IDs(ind(repeat_pval_keyframe<0.05))),nnz(repeat_pval_keyframe<0.05),length(ind));
fprintf('Subjects can separate frames (p<0.05): %s (%i of %i)\n',cell2str(subject_ID(repeat_corr_count_pval_binom<0.05)),nnz(repeat_corr_count_pval_binom<0.05),length(subject_ID));  

badcount = nan(1,length(importance_image_IDs));
for i=1:length(importance_image_IDs)
    [s,badcount(i)] = convert_str2num(importance_responses(:,i),importance_str);
    s(isnan(s))=[];
    %repeat_pval(i)=myBinomTest(nnz(s==2),length(s),0.5,'right');
    importance_percentage(i) = nnz(s==3)/length(s);
    %[~,importance_pval(i)] = chi2gof(s,'Ctrs',[1:length(importance_str)],'Expected',[1:length(importance_str)]/length(importance_str));
end
ind = find(ismember(importance_image_IDs,KEYFRAME_IMAGE_IDs));
fprintf('\n---SCENE IMPORTANCE---\n');
fprintf('Important scenes (all): %s\n',num2str(importance_image_IDs(importance_percentage>0.5)));
%fprintf('Important keyframes: %s\n',num2str(importance_image_IDs(importance_percentage(ind)>0.5)));

badcount = nan(1,length(eureka_image_IDs));
for i=1:length(eureka_image_IDs)
    [s,badcount(i)] = convert_str2num(eureka_responses(:,i),eureka_str);
    s(isnan(s))=[];
    %repeat_pval(i)=myBinomTest(nnz(s==2),length(s),0.5,'right');
    eureka_percentage(i) = nnz(s==3)/length(s);
    %[~,importance_pval(i)] = chi2gof(s,'Ctrs',[1:length(importance_str)],'Expected',[1:length(importance_str)]/length(importance_str));
end
ind = ismember(eureka_image_IDs,KEYFRAME_IMAGE_IDs);
%eureka_image_IDs(eureka_percentage(ind)>0.5)
fprintf('\n---EUREKA MOMENTS---\n');
fprintf('Eureka scenes (all): %s\n',num2str(eureka_image_IDs(eureka_percentage>0.5)));


%--
fprintf('\n---KEYFRAME ARRAY---\n');

for i=1:size(repeat_responses_num,1)
   fprintf('%s;',subject_ID{i});
   for j=1:size(repeat_responses_num,2)
       if repeat_responses_num(i,j)==2
        fprintf('o;',repeat_responses_num(i,j));
       else
           fprintf(';',repeat_responses_num(i,j));
       end
   end
   fprintf('%f',repeat_corr_count_pval_binom(i));
   fprintf('\n');   
end
fprintf(';');
for j=1:size(repeat_responses_num,2)         
%    if repeat_pval(j)<0.05
%         fprintf('1');
%    else
%        fprintf('0');
%    end
   
   fprintf('%f',repeat_pval(j));
   
   if j<size(repeat_responses_num,2)
      fprintf(';');
   end
end
fprintf('\n');




