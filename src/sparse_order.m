function [ idx ] = sparse_order( i,j, Asize )
%SPARSE_ORDER obtain the storing order of a sparse matrix
%   i: row idx
%   j: column idx
%   Asize: size of the sparse matrix

nnz = length(i);
spm = sparse(i,j,1:nnz,Asize(1), Asize(2));

idx = obtain_sp_order(spm);

end

