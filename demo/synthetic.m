
function [ X, X_com ] = synthetic( ndims, r, data_type, work_dims)
%SYNTHETIC Generate low rank synthetic data (boolean)
%   This function generate a synthetic data set which 
%   is a tensor of size ndims with low rank r
%   data_type represents the low rank type:
%       'MIX': a mixture of tensors different mode rank
%       'CORE': a tensor has mode rank r
%   X is the original data with components in X_com (only exists when
%   data_type = 'MIX'

nmodes = length(ndims);

X = tensor(zeros(ndims));
X_com = cell(nmodes, 1);

if exist('work_dims', 'var')
    dims = work_dims;
else
    dims = 1:nmodes;
end

switch data_type
    case 'MIX'
     for d = dims
        n = ndims(d);
        m = prod(ndims)/n;
        A = rand(n, r(d));
        A = A - mean(A(:));
        B = rand(r(d), m);
        B = B - mean(B(:));
        rdims = d;
        cdims = 1:nmodes;
        cdims = cdims(cdims ~= d);
        C = A*B;
        C = C - mean(mean(C));
        X_com{d} = tensor(tenmat(C, rdims, cdims, ndims));
        X = X + X_com{d};
     end
    case 'CORE'
        core_t = tensor(rand(r));
        core_t = core_t - mean(core_t(:));
        A = cell(nmodes,1);
        for d = 1:nmodes
            A{d} = rand(ndims(d), r(d));
            A{d} = A{d} - mean(A{d}(:));
        end
        X = ttm(core_t, A, 1:nmodes);
        X = X - mean(X(:));
    otherwise
        error('Unexpected data type');
end
 
X = X - mean(X(:));
X = X/std(X(:));

X = double(X);
X(X>0) = 1;
X(X<=0) = -1;
X = tensor(X);
end
