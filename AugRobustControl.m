function [A,DA,B,Q,H] = AugRobustControl(A0,DA0,B0,Q0,H0,delay,delta)

% [A,DA,B,Q,H] = AugRobustControl(A0,DA0,B0,Q0,H0,delay,delta)
%
% Augment the system matrices to take the feedback delay into account. 
%   [A0, DA0, B0, Q0, H0] are the state space representation matrices 
%   without delay. The output matrices include the delay. 
%
%   DELAY: hard time shift in the closed loop system in ms. 
%   DELTA: discretization step. 
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances. 
%   Crevecoeur F., Scott S. H., Cluff T. 
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

h = floor(delay/delta); %Feedback delay in number of sample times

n = size(A0,1);
m = size(B0,2);
t = size(Q0,3);
p = size(H0,1);

A = zeros((h+1)*n,(h+1)*n);
DA = zeros((h+1)*n,(h+1)*n);
B = zeros((h+1)*n,m);
Q = zeros((h+1)*n,(h+1)*n,t);
H = zeros(p,(h+1)*n);

A(1:n,1:n) = A0;
A(n+1:end,1:end-n) = eye(h*n);
DA(1:n,1:n) = DA0;
B(1:n,:) = B0;
H(:,end-n+1:end) = H0;

% Adding h times the constraint Q1:
Qaug = zeros(n,n,t+h);
for i = 1:h
    Qaug(:,:,i) = Q0(:,:,1);
end
for i = 1:t
    Qaug(:,:,i+h) = Q0(:,:,i);
end

%Filling the diagonal Q matrices
for time = 1:t
    for ii = 0:h
        Q(ii*n+1:(ii+1)*n,ii*n+1:(ii+1)*n,time) = Qaug(:,:,time+h-ii)/(h+1);
    end
end





