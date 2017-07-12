function [split_timearray,inds] = split_per_session(timearray)

TR=1.56;

% first pulse is always ~0.18 seconds behind video start
% ALSO 
% slice-timing correction adds another delay (currently half-TR)
first_volume_movie_difference = 0.18 + TR/2;

ses1_offset=0;
ses2_offset=2125.0;
ses3_offset=4242.24;

for i=1:length(timearray)
   ind=find(timearray<ses2_offset);
   split_timearray{1} = timearray(ind) - ses1_offset - first_volume_movie_difference;
   inds{1}=ind;

   ind=find(timearray>=ses2_offset & timearray<ses3_offset);
   split_timearray{2} = timearray(ind) - ses2_offset - first_volume_movie_difference;
   inds{2}=ind;
   
   ind=find(timearray>=ses3_offset);
   split_timearray{3} = timearray(ind) - ses3_offset - first_volume_movie_difference;
   inds{3}=ind;
end

end