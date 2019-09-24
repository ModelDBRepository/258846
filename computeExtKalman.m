function [K,JC] = computeExtKalman(A,B,C,D,H,oXi,oOmega,oEta,L,cxhat,ce)

%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances. 
%   Crevecoeur F., Scott S. H., Cluff T. 
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

n = size(A,1);
k = size(H,1);
d = size(D,3);
c = size(C,3);
step = size(L,3);
sigmaE = ce;
sigmaX = cxhat;
sigmaEX = zeros(n);

K = zeros(n,k,step);
JC = zeros(step,3);

for i = 1:step
    
    statedn = 0;
    sTemp = (sigmaE+sigmaX+sigmaEX+sigmaEX');
    
    for j = 1:d
        
        statedn = statedn + D(:,:,j)*sTemp*D(:,:,j)';
        
    end
    
    sdn = 0;
    
    for j = 1:c
        
        sdn = sdn + C(:,:,j)*L(:,:,i)*sigmaX*L(:,:,i)'*C(:,:,j)';
        
    end 
    
    K(:,:,i) = A*sigmaE*H'/(H*sigmaE*H'+oOmega+statedn);

    sigmaETemp = sigmaE;
    sigmaE = oXi+oEta+(A-K(:,:,i)*H)*sigmaE*A'+sdn;
    term = (A-B*L(:,:,i))*sigmaEX*H'*K(:,:,i)';
    sigmaX = oEta + K(:,:,i)*H*sigmaETemp*A'+(A-B*L(:,:,i))*sigmaX*(A-B*L(:,:,i))'...
        + term + term';
    sigmaEX = (A-B*L(:,:,i))*sigmaEX*(A-K(:,:,i)*H)'-oEta;
    
    JC(i,:) = [sigmaE(1,1),sTemp(1,1),0];
    
end


    

