%Creates the A matrix for problem 5!
function [A,b,u]=SetupProb2b(n,p)
% This function sets up the A, b matrices for part 2 of the project. This
% function creates the A matrix based on the function a_k=|k+1|^-p.

theta=-pi+pi/10000:(2*pi+pi/10000)/(n):pi;
if(max(size(theta))<n)
    theta=-pi+pi/10000:(2*pi+pi/10000)/(n+1):pi;
end
vals=fft(theta);
a=real(vals);
b=imag(vals);
A=toeplitz([a]); % generates a symmetric toeplitz matrix based on the row given.
b= randn(n,1);




end