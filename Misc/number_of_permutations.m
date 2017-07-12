function res = number_of_permutations(N)
% Given two item sets (sizes N), this functions returns the total number of
% unique permutations when these sets are mixed

% i.e., N=4, original set is [1,1,1,1,2,2,2,2]
% two permutations are [1,2,1,1,2,2,1,2], [1,1,2,2,1,2,2,1]
% one permutation is always the original set (unchanged)

res = (2.^(2*N).*gamma(1/2+N))./(sqrt(pi)*gamma(1+N));


end

