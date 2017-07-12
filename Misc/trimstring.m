function res = trimstring(str)

    L=length(str);
    k1=1;
    
    for i=1:L
       if strcmp(str(i),' ')
            k1=k1+1;
       else
           break;
       end
    end
    
    k2=L;
    for i=L:-1:1
       if strcmp(str(i),' ')
            k2=k2-1;
       else
           break;
       end
    end
    
    res = str(k1:k2);

end