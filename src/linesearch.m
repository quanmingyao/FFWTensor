function [ gamma ] = linesearch( X, newcomp, X2 )
%LINESEARCH find the best gamma
%   X: X^t
%   newcomp: new component
%   X2: A.vals - X.vals

cpx = newcomp - X.vals;
center = -(-X2.vals)' * cpx /(norm(cpx)^2);

if center < 0
    gamma = 0;
else
    if center > 1
    gamma = 1;
    else
        gamma = center;
    end
end
end
