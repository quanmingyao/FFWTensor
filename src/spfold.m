function [ o_subs, Bsize ] = spfold( subs, rdim, cdims, Asize_i )
%SPFOLD obtains nonzero subs of the unfolding tensor (map to matrix)
%   subs is the nozero subs of the tensor
%   rdim: the dim mapped to row (unfolding mode)
%   cdims: the remaining dims, mapped to column
%   Asize_i: size of the original tensor

% o_subs: mapped subs
% Bsize: mapped size

%obtain the remaining size
Asize = Asize_i(cdims);

o_subs = zeros(length(subs(:,1)), 2);
o_subs(:,1) = subs(:,rdim);
mult = [1 cumprod(Asize(1:end-1))];
o_subs(:,2) = (subs(:,cdims) - 1) * mult' + 1;

Bsize = [Asize_i(rdim), prod(Asize)];
end

