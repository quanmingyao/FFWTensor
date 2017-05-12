function [ sigma, u, v ] = powermethod( X, options )
%POWERMETHOD Find the largest singular value/vector
%   Detailed explanation goes here

if ~isfield(options,'maxIter')
    maxIter = 100;
else
    maxIter = options.maxIter;
end

[n, m] = size(X);

z = rand(m, 1);

y = X * z;
y = y/norm(y,2);

if isfield(options,'u')
    y = options.u;
end

for t = 1:maxIter
    tmp = X' * y;
    y = X * tmp;
    normy = norm(y,2);
    if normy ~= 0
        y = y/normy;      
    end
end
    
    b = X' * y;
    sigma = norm(b,2);
    if sigma ~= 0;
        v = b/sigma;
    else
        v = b;
    end
    u = y;
    
    
end

