function [ data_train, P_train, data_test, P_test, data_mean, data_std ] = tensplitspdata(data, ratio)
%TENSPLITSPDATA normalizes the data to zero mean, unit std, and 
% splits the sparse data into a trainning part and testing part by randomly
% selection
%   data: the sparse tensor data
%   ratio: the ratio of training data
%   outputs are all sptensor

data_mean = mean(data.vals);
data_use = elemfun(data, @(x) x - data_mean);

data_std = std(data_use.vals);
data_use = data_use/data_std;

data_nnz = nnz(data_use);
train_nnz = round(data_nnz * ratio);


ind = randperm(data_nnz);
train_ind = ind(1:train_nnz);
test_ind = ind(train_nnz+1:end);

train_subs = data_use.subs(train_ind, :);
train_vals = data_use.vals(train_ind, :);

test_subs = data_use.subs(test_ind, :);
test_vals = data_use.vals(test_ind, :);

P_train = sptensor(train_subs, ones(size(train_subs,1), 1), size(data));
P_test = sptensor(test_subs, ones(size(test_subs,1), 1), size(data));

data_train = sptensor(train_subs, train_vals, size(data));
data_test = sptensor(test_subs, test_vals, size(data));
end
