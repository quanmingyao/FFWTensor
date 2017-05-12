function [ o_subs ] = spunfold( i_subs, rdim, cdims, Asize )
%SPUNFOLD Summary of this function goes here
%   Detailed explanation goes here
[N,~] = size(i_subs);

Asize = Asize(cdims);

k = [1 cumprod(Asize(1:end-1))];
n = length(Asize);

subs = i_subs(:,2);
subs = subs - 1;

o_subs = zeros(N, numel(cdims)+1);

for i = n : -1 : 1
    o_subs(:,cdims(i)) = floor(subs / k(i)) + 1;
    subs = rem(subs,k(i));

end
o_subs(:,rdim)  = i_subs(:,1);
