function res = parse_timing_file(FILE,TYPE)

if nargin<2
    TYPE=1;
end

fid = fopen(FILE);

data=[];
i=1;
tline = fgetl(fid);
data{1} = tline;
while ischar(tline)    
    i=i+1;
    tline = fgetl(fid);
    if isstr(tline)
        data{i}=tline;
    end
end
fclose(fid);

N=length(data);

segmentdata=struct();

k=0;
state = 0;

if TYPE==1
    
    for i=1:N
        %%fprintf('...line %i\n',i);
        str = strsplit(data{i},'\t');
        if state==0
            if length(str)==1 && ~isempty(str{1})
                k=k+1;
                segmentdata(k).ID = str2num(str{1});
                state = 1;
            end
        elseif state==1
            if ~isempty(str) && length(str)==2
                segmentdata(k).ORIGINAL_ses = str2num(str{1});
                [dt,st] = parsetime(str{2},TYPE);
                segmentdata(k).ORIGINAL_time = dt;
                segmentdata(k).ORIGINAL_time_str = st;
                if diff(segmentdata(k).ORIGINAL_time)<1
                    segmentdata(k).ID
                    segmentdata(k).ORIGINAL_time
                    error('Bad timing (too short or negative segment time)')
                end
                state = 2;
            else
                i
                error('Bad data!')
            end
        elseif state==2
            if ~isempty(str) && length(str)==2
                segmentdata(k).KRON_ses = str2num(str{1});
                [dt,st] = parsetime(str{2},TYPE);
                segmentdata(k).KRON_time = dt;
                segmentdata(k).KRON_time_str = st;
                if diff(segmentdata(k).KRON_time)<1
                    segmentdata(k).ID
                    segmentdata(k).KRON_time
                    error('Bad timing (too short or negative segment time)')
                end
                
                segmentdata(k).time_difference = abs(diff(segmentdata(k).ORIGINAL_time)-diff(segmentdata(k).KRON_time));
                
                state = 0;
                
                fprintf('Segment ID = %i:\n   ORIG. duration = %.2fs, CHRONO. duration = %.2fs (difference %.2fs)\n   ORIG.   ses = %i, timeslot %s\n   CHRONO. ses = %i, timeslot %s\n',...
                    segmentdata(k).ID,...
                    diff(segmentdata(k).ORIGINAL_time),...
                    diff(segmentdata(k).KRON_time),...
                    abs(diff(segmentdata(k).ORIGINAL_time)-diff(segmentdata(k).KRON_time)),...
                    segmentdata(k).ORIGINAL_ses,...
                    segmentdata(k).ORIGINAL_time_str,...
                    segmentdata(k).KRON_ses,...
                    segmentdata(k).KRON_time_str);
                
            else
                error('Bad data!')
            end
        end
    end
    
elseif TYPE==2
    
    for i=1:N
        %%fprintf('...line %i\n',i);
        str = strsplit(data{i},'\t');
        if state==0
            if length(str)==1 && ~isempty(str{1})
                k=k+1;
                segmentdata(k).ID = str2num(str{1});
                state = 1;
            end
        elseif state==1
            if ~isempty(str) && length(str)==2
                segmentdata(k).FIRST_ses = str2num(str{1});
                [dt,st] = parsetime(str{2},TYPE);
                segmentdata(k).FIRST_time = dt;
                segmentdata(k).FIRST_time_str = st;
                state = 2;
            else
                i
                error('Bad data!')
            end
        elseif state==2
            if ~isempty(str) && length(str)==2
                segmentdata(k).SECOND_ses = str2num(str{1});
                [dt,st] = parsetime(str{2},TYPE);
                segmentdata(k).SECOND_time = dt;
                segmentdata(k).SECOND_time_str = st;
                                
                state = 0;
                
%                 fprintf('Segment ID = %i:\n   ORIG. duration = %.2fs, CHRONO. duration = %.2fs (difference %.2fs)\n   ORIG.   ses = %i, timeslot %s\n   CHRONO. ses = %i, timeslot %s\n',...
%                     segmentdata(k).ID,...
%                     diff(segmentdata(k).ORIGINAL_time),...
%                     diff(segmentdata(k).KRON_time),...
%                     abs(diff(segmentdata(k).ORIGINAL_time)-diff(segmentdata(k).KRON_time)),...
%                     segmentdata(k).ORIGINAL_ses,...
%                     segmentdata(k).ORIGINAL_time_str,...
%                     segmentdata(k).KRON_ses,...
%                     segmentdata(k).KRON_time_str);
                
            else
                error('Bad data!')
            end
        end
    end    

elseif TYPE==3
    
    for i=1:N
        %%fprintf('...line %i\n',i);
        str = strsplit(data{i},'\t');
        if state==0
            if length(str)==1 && ~isempty(str{1})
                k=k+1;
                segmentdata(k).ID = str2num(str{1});
                state = 1;
            end
        elseif state==1
            if ~isempty(str) && length(str)==2
                segmentdata(k).FIRST_ses = str2num(str{1});
                [dt,st] = parsetime(str{2},TYPE);
                segmentdata(k).FIRST_time = dt;
                segmentdata(k).FIRST_time_str = st;
                state = 0;
            else
                i
                error('Bad data!')
            end
        end
    end
    
end

res = segmentdata;

end

function [dt,st] = parsetime(str,TYPE)
    str = strsplit(str,'-');
    if length(str)~=2
        str = strsplit(str{1},'–');
    end
    if TYPE==1 && length(str)~=2
        str
        error('Bad timestring (separator problem)')
    end
        
    str1=strtrim(str{1});           
    str0 = strsplit(str1,':');
    dt(1) = str2num(str0{1})*60*60 + str2num(str0{2})*60 + str2num(str0{3}) + str2num(str0{4})/100;
    st = str1;
    
    if length(str)>1
        str2=strtrim(str{2});
        str0 = strsplit(str2,':');
        dt(2) = str2num(str0{1})*60*60 + str2num(str0{2})*60 + str2num(str0{3}) + str2num(str0{4})/100;
        st = [str1,' -> ',str2];
    end

    
end