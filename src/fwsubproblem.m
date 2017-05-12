function [ j, sigma, u, v ] = fwsubproblem( X, dims, sppattern, sppattern_index, fac_size )
%FWSUBPROBLEM solve the subproblem of frank-wolfe
%   Detailed explanation goes here

D = length(dims);

ul = cell(D,1);
vl = cell(D,1);
sigmal = zeros(D,1);

X2 = cell(D,1);
for d = dims

rdim = d;
cdims = 1:X.nmodes;
cdims = cdims(cdims ~=rdim);

X2{d} = setsparse(sppattern{d}, X.vals(sppattern_index{d}));

options.maxIter = 3;
[sigmal(d), ul{d}, vl{d}] = powermethod(X2{d},options);

sigmal(d) = sigmal(d) * fac_size(d);
end
[~, j] = max(sigmal);
u = ul{j};

options.u = u;
options.maxIter = 8;  %for large tensor
[sigma, u,v] = powermethod(X2{j}, options);

end

