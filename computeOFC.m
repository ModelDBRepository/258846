function [L,Sx,Se,s] = computeOFC(A,B,C,D,H,Q,R,K,oXi,oOmega,oEta)

%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances. 
%   Crevecoeur F., Scott S. H., Cluff T. 
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

n = size(A,1);
m = size(B,2);
c = size(C,3);
d = size(D,3);
step = size(R,3);
L = zeros(m,n,step);

currSx = Q(:,:,end);
currSe = 0;
currs = 0;

for i = step:-1:1
    
    sdn = 0;
    
    for j = 1:c
        
        sdn = sdn + C(:,:,j)'*(currSx + currSe)*C(:,:,j);
        
    end
    
    statedn = 0;
    
    for j = 1:d
        
        statedn = statedn + D(:,:,j)'*K(:,:,i)'*currSe*K(:,:,i)*D(:,:,j);
        
    end

    L(:,:,i) = (R(:,:,step) + B'*currSx*B + sdn)\(B'*currSx*A);
    currSxTemp = currSx;
    currSx = Q(:,:,i) + A'*currSx*(A-B*L(:,:,i)) + statedn;
    currSeTemp = currSe;
    currSe = A'*currSxTemp*B*L(:,:,i)+...
        (A-K(:,:,i)*H)'*currSeTemp*(A-K(:,:,i)*H);
    currs = trace(currSxTemp*oXi+currSeTemp*(...
        oXi+oEta+K(:,:,i)*oOmega*K(:,:,i)'))+currs;
    
end

Sx = currSx;
Se = currSe;
s = currs;