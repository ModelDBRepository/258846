function [L,K,newCost,Sx,Se,s,JC] = extLQG(A,B,C,D,H,Q,R,cxi,comega,ceta,ce,cxhat,xInit)

% [L,K,newCost,Sx,Se,s,sall,JC] = extLQG(A,B,C,D,H,Q,R,cxi,comega,ceta,ce,cxhat,cxe,xInit)
% 
% Calculates a series of feedback gains and Kalman gains for extended LQG
% control. The routine follows (Todorov, 2005, Neural Comp, 17, 1084-1108).
%   [A, B, C, D, H]     matrices for the state pace representation
%   [Q, R]              cost matrices
%   CXI, COMEGA, CETA   motor, sensory and internal noise
%                       covariance matrices
%   CE, CXHAT           initial covariance matrices for the error
%                       and the estimation
%   XINIT               initial state vector 
%
% The output matrices are:
%   L                   time series of optimal feedback gains
%   K                   time series of non adaptive Kalmen gains
%   NEWCOST             expected cost after stopping the iterations
%   SX, SE              series of cost matrices determined by the backwards
%                       recurrences for the state (SX) and error (SE) terms
%   S                   scalar component of the total expected cost
%   JC                  covariance matrix
%
%
%   Uses: > computeOFC
%         > computeExtKalman
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances. 
%   Crevecoeur F., Scott S. H., Cluff T. 
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019


n = size(A,1);
m = size(B,2);
p = size(H,1);
c = size(C,3);
step = size(R,3);

K = zeros(n,p,step);

tol = 10^-14;
current = 10^6;
itmax = 100;
count = 0;
found = false;


% The optimal control and Kalman gains are calculated iteratively (no more
% than itmax times if it does not converge)
while ~found && count < itmax
    
    [L,Sx,Se,s] = computeOFC(A,B,C,D,H,Q,R,K,cxi,comega,ceta);
    [K,JC] = computeExtKalman(A,B,C,D,H,cxi,comega,ceta,L,cxhat,ce);

    % Expected cost
    newCost = xInit'*Sx*xInit + trace((Sx+Se)*ce)+ s;
    
    % Relative improvement
    dCost = abs(current - newCost)/newCost;
    current = newCost;

    % Check: of the relative improvement is small, the solution is found
    if dCost > tol
        found = false;
    else
        found = true;
    end
    
    count = count + 1;
    
end

fprintf('Number of iterations (ext LQG): %d\n', count); % Number of iterations 
