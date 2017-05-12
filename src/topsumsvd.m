function [U, S, V, t]  = topsumsvd(A, tau)
%TOPSUMSVD Summary of this function goes here
%   Detailed explanation goes here

[U, S, V] = svd(A, 0);
s = diag(S);

sum_s = sum(s);
num_s = sum(s>0);
s = s(1:num_s);

while(sum_s > tau)    
    s = s - (sum_s - tau)/num_s;    
    num_s = sum(s>0);
    s = s(1:num_s);
    sum_s = sum(s); 
end

U = U(:,1:num_s);
V = V(:,1:num_s);
S = diag(s);
t = num_s;

end

