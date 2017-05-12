function [ auc ] = calAUC( X, data_test, k)
%CALAUC calculate AUC 
%   X, data_test are vectors of the same length, containing nnz elements

exist_link = find(data_test > 0);
nonexist_link = find(data_test < 0);

l_exist = length(exist_link);
l_nonexist = length(nonexist_link);

ind = randi(l_exist, k, 1);
exist_ind = exist_link(ind);

ind = randi(l_nonexist, k, 1);
nonexist_ind = nonexist_link(ind);

auc = (sum(X(exist_ind) > X(nonexist_ind)) + 1/2 * sum(X(exist_ind) == X(nonexist_ind)))/k;

end

