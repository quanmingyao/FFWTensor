function [ Sigma_p ] = post_process( U, V, k, Pi, Pj, dims, A )
%POST_PROCESS do post process 
%   Detailed explanation goes here
n_basis = sum(k);
W = zeros(length(Pi{1}), n_basis);

% obtain W
i = 1;
for d = dims
    for s = 1:k(d)
        i = i + 1;
        W(:, i) = spmultic(Pi{d}, Pj{d}, U{d}(:,1:k(d)) , V{d}(:,1:k(d))');
    end
end

B = W' * W;
a = W' * A.vals;

Sigma_p = diag(a/B);

end

