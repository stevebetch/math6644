%Creates the A matrix for problem 5!
function [A,b]=SetupProb2b(n,p)
% This function sets up the A, b matrices for part 2 of the project. This
% function creates the A matrix based on the function a_k=|k+1|^-p.

theta=-pi:(2*pi)/(n-1):pi;
if(max(size(theta))<n)
    theta=-pi:(2*pi)/(n):pi;
end
vals=fft(theta);
a=real(vals);
b=imag(vals);
A=toeplitz(vals); % generates a symmetric toeplitz matrix based on the row given.
b= randn(length(A),1);




end