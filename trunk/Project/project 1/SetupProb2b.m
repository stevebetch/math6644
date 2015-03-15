%Creates the A matrix for problem 5!
function [A,b,u]=SetupProb5b(n)
% This function sets up the A, b matrices for part 2 of the project. This
% function creates the A matrix based on the function a_k=|k+1|^-p.

theta=-pi:(2*pi)/(n-1):pi;
vals=fft(theta);
a=real(vals);
b=imag(vals);
A=toeplitz([a,b]); % generates a symmetric toeplitz matrix based on the row given.
b= randperm(100,n)';
b(end+1:end+n,1)=0;




end