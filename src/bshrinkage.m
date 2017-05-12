function [Ul, Sigmal, Vl, k] = bshrinkage(A, Ul, Sigmal, Vl, tau, dims, nmodes, Asize, Pi, Pj, rho, k )
%SHRINKAGE the shrinkage algorithm

% A: observation
% Ul Sigmal Vl:  basis
% tau: tau in the model
% taus: nuclear norm of each component
% dims: working dimensions
% nmodes: mode number
% Asize: size of the tensor
% Pi, Pj:  indices of nonzero entries in each mode unfolding
% rho: proximal step size
% k: contains basis numbers

PUSV = cell(nmodes, 1);
PUSVall.vals = zeros(size(A.vals));

taus = zeros(nmodes, 1);
U = cell(nmodes, 1);
V = cell(nmodes, 1);
S = cell(nmodes, 1);
eps = 1;

for d = dims
    U{d} = Ul{d}(:, 1:k(d));
    V{d} = Vl{d}(:, 1:k(d));
    S{d} = diag(Sigmal{d}( 1:k(d)));
    
    US = U{d} * S{d};
    
    PUSV{d}.vals = spmultic(Pi{d}, Pj{d}, US, V{d}'); 
    PUSVall.vals = PUSVall.vals + PUSV{d}.vals;
    
    [U{d}, RU] = qr(U{d}, 0);
    [V{d}, RV] = qr(V{d}, 0);
    S{d} = RU * S{d} * RV';
    
    [~, s, ~] = svd(S{d});
    taus(d) = sum(diag(s));
end

PUSVall.vals =  PUSVall.vals - A.vals;

for d = dims
   
    if(k(d) < 1)
        continue;
    end
    B = sparse(Pi{d}, Pj{d}, PUSVall.vals, Asize(d), prod(Asize)/Asize(d));
    B =  B * V{d};
    B = U{d}' * B;
    
    C = S{d} - rho(d) * B;

    if ~isempty(C)
        [Us, Ss, Vs, t] = topsumsvd(C, taus(d)*eps);
    end    
    
    tmp = U{d} * Us;
    Ul{d} = zeros(size(Ul{d}));
    Ul{d}(:, 1:t)  = tmp;
    
    tmp = V{d} * Vs;
    Vl{d} = zeros(size(Vl{d}));
    Vl{d}(:, 1:t) = tmp;
    
    Sigmal{d} = zeros(size(Sigmal{d}));
    Sigmal{d}(1:t) = diag(Ss);    
    
    k(d)  = t;
    
    PUSVall.vals = PUSVall.vals - PUSV{d}.vals;
    PUSV{d}.vals = spmultic(Pi{d}, Pj{d}, Ul{d}(:,1:t)*Ss, Vl{d}(:,1:t)'); 
    PUSVall.vals = PUSVall.vals + PUSV{d}.vals;
end
end
