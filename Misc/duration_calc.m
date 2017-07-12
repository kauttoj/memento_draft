function out = duration_calc(str1,str2)
%DURATION_CALC Summary of this function goes here
%   Detailed explanation goes here

[dt1,st1] = parsetime(str1,1);
[dt2,st2] = parsetime(str2,2);

out = sec2min(dt2+diff(dt1));

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