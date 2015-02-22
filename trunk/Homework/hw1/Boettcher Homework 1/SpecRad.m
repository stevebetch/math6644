function [ G,eigVals ] = SpecRad( A )
%SPECRAD This function calculates the spectral radius of matrix A for hw1.
%   Detailed explanation goes here

M=diag(diag(A));
N=-(tril(A)-diag(diag(A)))-(triu(A)-diag(diag(A)));

G=M^-1*N;

eigVals=eigs(G);


end

