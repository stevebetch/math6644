%Creates the A matrix for problem 5!
function [A,b]=SetupProb2a(n,p)
% This function sets up the A, b matrices for part 2 of the project. This
% function creates the A matrix based on the function a_k=|k+1|^-p.

x=1:n;
ak=abs(x+1).^-p;
A=toeplitz(ak); % generates a symmetric toeplitz matrix based on the row given.
b= randn(n,1);




end