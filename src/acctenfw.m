function [  out, obj, info ] = acctenfw( data_train, tau, data_test, test_P, opts )
%ACCTENFW Summary of this function goes here
% INPUT:
%   data_train: training data (sptensor)
%   data_test: testing data (sptensor)
%   tau: tau in the model (nuclear norm contraint)
%   Xout: a structure contains U, Sigma(S), V (solution)
%   opts:
%       maxIter: max iteration number
%       ex_dims: exclude unfolding dimension
%       verbose: control display information
%       err_type: err_type ['RMSE', 'AUC'] 
%       maxK: maximum bases nubmer (pre-allocate space) 
%       shrinkIter: maximum basis number to perform bs 
%       rho_s:  rho in bs
%       accelerated: not used, set 0
%       compute_test_per: whether to compute test performance
%       sample size to estimate AUC

if isempty(opts), opts=struct('a',1) ; end

function out = set_opts( opts, field, default )
    if ~isfield( opts, field )
        opts.(field)    = default;
    end
    out = opts.(field);
end

shrinkIter = set_opts(opts, 'shrinkIter', 400);
maxIter = set_opts(opts, 'maxIter', 200);
maxK = set_opts(opts, 'maxK', 400 * ones(ndims(data_train), 1));
ex_dims = set_opts(opts, 'ex_dims', 3);
l_search = set_opts(opts, 'l_search', 1);
compute_test_per = set_opts(opts, 'compute_test_per', 1);
accelerated = set_opts(opts, 'accelerated', 0);
post_processed = set_opts(opts, 'post_processed', 1);
verbose = set_opts(opts, 'verbose', 1);
AUC_size = set_opts(opts, 'AUCsize', 1e6);
err_type = set_opts(opts, 'err_type', 1);
rho_s = set_opts(opts, 'rho', 0.4);

A.vals = data_train.vals;
A.nmodes = ndims(data_train);
A.size = size(data_train);
A.subs = data_train.subs;

D = A.nmodes;
Asize = A.size;
nnzA = length(A.subs);

fac_size = sqrt(Asize);
dims = 1:D;
dims = setdiff(dims, ex_dims);

nmodes = A.nmodes;
rho = rho_s * ones(nmodes, 1);

data_test = data_test.vals;

Sigma = cell(nmodes,1);
YSold = cell(nmodes, 1);
YSnew = cell(nmodes, 1);

U = cell(nmodes, 1);
V = cell(nmodes, 1);

Pi = cell(nmodes, 1);
Pj = cell(nmodes, 1);
sppattern = cell(nmodes, 1);
sppattern_index = cell(nmodes, 1);

for d = dims
    Sigma{d} = zeros(maxK(d), 1);
    U{d} = zeros(Asize(d), maxK(d));
    V{d} = zeros(prod(Asize)/Asize(d), maxK(d));
    
    YSold{d} = zeros(maxK(d), 1);
    
    rdim = d;
    cdims = 1:nmodes;
    cdims = cdims(cdims ~= rdim);
    
    transP = spfold(A.subs, rdim, cdims, Asize);
    Pi{d} = transP(:,1);
    Pj{d} = transP(:,2);
    
    transP_test = spfold(test_P.subs, rdim, cdims, Asize);
    Pi_test{d} = transP_test(:, 1);
    Pj_test{d} = transP_test(:, 2);
    
    sppattern{d} = sparse(Pi{d}, Pj{d}, ones(size(Pi{d})), Asize(d), prod(Asize)/Asize(d));
    sppattern_index{d} = sparse_order(Pi{d}, Pj{d}, [Asize(d), prod(Asize)/Asize(d)]);
end

k = zeros(nmodes, 1);
kr = zeros(nmodes, maxIter);

dualgap = zeros(maxIter, 1);
obj = zeros(maxIter, 1);

iter_time = zeros(maxIter, 1);
test_per = zeros(maxIter, 1);

X2.vals = zeros(nnzA,1);
X2.subs = A.subs;

X.vals = zeros(nnzA,1);
X.subs = A.subs;

Ynew = X;
Yold = X;

shrinkage = 0;

lambda_old = 1;

disp('acctenfw start running');
for t = 1:maxIter
    
    tic;
    for d = dims
        if (sum(k) > shrinkIter)
            shrinkage = 1;
        end
    end

    gamma = 2/(t+2);  %stepsize
    
    X2.vals = A.vals - X.vals;
    
    X2.nmodes = nmodes;
    X2.size = Asize;
    [j, ~, u, v] = fwsubproblem(X2, dims, sppattern, sppattern_index, fac_size);

    newcomp =  tau * fac_size(j)*  spmultic(Pi{j}, Pj{j}, u, v');    
    
    if l_search
        gamma = linesearch(X, newcomp, X2);
    end
    if gamma == 0
        gamma = 1e-6;
    end
    
    k(j) = k(j) + 1;
    
    kr(:,t) = k;
    
    if accelerated
        lambda_new = (1+sqrt(1+4*lambda_old^2))/2;
        gs = (lambda_old-1)/lambda_new;
        
        Ynew.vals  = (1-gamma)*X.vals + gamma * (newcomp);
        for d = dims
            YSnew{d} = (1-gamma) * Sigma{d}; 
        end
        YSnew{j}( k(j)) = gamma * tau;
        
        X.vals =  Ynew.vals + gs *(Ynew.vals -  Yold.vals);
        for d = dims
            Sigma{d} = YSnew{d} + gs *(YSnew{d} - YSold{d});
        end
        
        Yold = Ynew;
        % add new basis
        for d = dims
            YSold{d} = YSnew{d};
        end
        
        lambda_old = lambda_new;
        if norm((X.vals - A.vals), 'fro') > norm((Ynew.vals-A.vals), 'fro')
            disp('not accelerate');
            accelerated = 0;
            X.vals = Ynew.vals;           
            for d = dims
                Sigma{d} = YSnew{d};
            end    
        end               
    else
        X.vals  = (1-gamma)*X.vals + gamma * (newcomp);
        % add new basis
        for d = dims
            Sigma{d} = (1-gamma) * Sigma{d};
        end
        Sigma{j}(k(j)) = gamma *fac_size(j) *  tau;
    end
      
    U{j}(:, k(j)) = u;
    V{j}(:, k(j)) = v;     
    
    dualgap(t) = ( X.vals - newcomp)' * (-X2.vals) ;
    obj(t) = norm(X2.vals)^2; 
    
    if shrinkage  
        disp(t);
        disp('shrink');
        [U, Sigma, V, k] = bshrinkage(A, U, (Sigma), V, tau, dims, nmodes, Asize, Pi, Pj, rho, k );

        X.vals = zeros(size(X.vals));
        for d = dims
            X.vals = X.vals + spmultic(Pi{d}, Pj{d}, U{d}(:,1:k(d)) * diag(Sigma{d}(1:k(d))), V{d}(:,1:k(d))');            
            
            YSold{d} = Sigma{d};
        end
        shrinkage = 0;
    end
    iter_time(t) = toc;
    
    if compute_test_per        
        
        Xout = zeros(nnz(test_P), 1);
        for d = dims            
            tmp = spmultic(Pi_test{d}, Pj_test{d}, U{d}(:,1:k(d)) * diag(Sigma{d}(1:k(d))), V{d}(:,1:k(d))');
            Xout = Xout + tmp;
        end

        switch err_type
            case 'RMSE'
                test_per(t) = sqrt(mean((Xout - data_test).^2));
                
            case 'AUC'                
                test_per(t) = calAUC(Xout, data_test, AUC_size);
            otherwise
                error('wrong performance evaluation method');
        end
        
    end
    if verbose
        fprintf('Iteration: %d, gamma: %.2e, obtained mode: %d, #bases: %d, \n',  t, gamma, j, sum(kr(:, t)));
        fprintf('    obj: %.2e, test %s: %.2e. , running time: %.2f \n', obj(t), err_type, test_per(t), sum(iter_time));
    end
    
    if sum(iter_time) > 1e3
        break;
    end
end

Sigma_p = Sigma;
Xout = zeros(nnz(test_P), 1);
for d = dims
    tmp = spmultic(Pi_test{d}, Pj_test{d}, U{d}(:,1:k(d)) * diag(Sigma_p{d}(1:k(d))), V{d}(:,1:k(d))');
    Xout = Xout + tmp;
end
switch err_type
    case 'RMSE'
        f_result = sqrt(mean((Xout - data_test).^2));
        
    case 'AUC'
        f_result = calAUC(Xout, data_test, AUC_size);
    otherwise
        error('wrong performance evaluation method');
end

out.U = U;
out.V = V;
out.S = Sigma;

iter_time = cumsum(iter_time);

info.f_result = f_result;

info.itime = iter_time;
info.kr = kr;
info.test_p = test_per;
info.iter = t;
info.dualgap = dualgap;
end
