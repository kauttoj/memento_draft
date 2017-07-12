function bad_volumes = find_overlapping_volumes(nonkey_volumes,key_volumes)

    bad_volumes=[];
    for i=1:length(nonkey_volumes)
        for j=1:length(key_volumes)
            if abs(nonkey_volumes(i)-key_volumes(j))<15
                bad_volumes(end+1)=i;
                break;
            end
        end
    end           

end