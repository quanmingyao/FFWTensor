function [ C ] = spmulti( Pi, Pj, A, B )
%SPMULTI  efficiently compute P(AB), here P is a sparse
%  sampling operator
%   Detailed explanation goes here

[n, ~] = size(A);
[~, m] = size(B);




[i, j, ~] = find(P);

[i_ind, ind_tmp] = sort(i);
j_ind = j(ind_tmp);
val = length(i);


% %iterate by row (i)
% i2 = unique(i);
% ind = 0;
% for t = i2'
%     
%     %find all position with row t
%     i_ind_tmp = (i==t);
%     sum_t = sum(i_ind_tmp);   
%     j_ind_tmp = j(i_ind_tmp);
%     
%     
%     val(ind+1: ind+ sum_t) = A(t, :) * B(:, j_ind_tmp);
%     ind = ind+sum_t;
% end

C = sparse(i_ind, j_ind, val, n, m);
end

