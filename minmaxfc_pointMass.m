function sout = minmaxfc_pointMass(xinit,xfinal,simdata)

% sout = minmaxfc_pointMass(xinit,xfinal,simdata)
%
% Calculates a trajectory with initial condition, final target and
% parameters defined in the input structure simdata.
%   XINIT: Initial State
%   XFINAL: Target State
%   Input data structure must contain:
%   SIMDATA.delta       Discretization step
%          .delay       Hard temporal delay in the closed loop system
%          .pert        1x2 vector with step force magnitude along x and y axes  
%          .time        Time horizon 
%          .gamma       1x2 with Parameter for optimal disturbance
%                       rejection level. The second entry (1 or 0) indicates 
%                       whether the routine should optimize this value.
%          .nStep       Number of time steps
%          .noise       1x2 vector of scaling parameters for noise matrices
%                       (default: [1 1])
%          .ralpha      matrix with one row per state variable and one column per time
%                       step with the cost of the corresponding state and time
%          .nsimu       number of simulation runs.
%           
%           
%   SOUT: output data structure with the following fields:
%   sout.L              Series of optimal robust control gains
%       .C              Series of optimal LQG control gains   
%       .x              State - Robust control
%       .xest           State Estimate - Robust control   
%       .z              State - LQG
%       .zest           State Estimate - LQG
%       .u              Series of Control Vector - Robust control
%       .v              Series of Control Vector - LQG
%       .minlambda      Minimum eigen value (optimized or used). Must be > 0.
%       .cost           1x2 vector with movement cost (1: Robust, 2: LQG)   
%       .gammaopt       Optimal or used gamma parameter    
%
%
%
%   Uses: > AugRobustControl
%         > extLQG
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances. 
%   Crevecoeur F., Scott S. H., Cluff T. 
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

delta = simdata.delta;
delay = simdata.delay;
gamma = simdata.gamma(1);
ralpha = simdata.ralpha;
nStep = simdata.nStep;
statedim = size(xinit,1);

% Mapping all final targets on 0 to ensure positive definiteness of the
% cost matrices
xinit = xinit-xfinal;

% System Matrices
k = .1;
tau = .066;
A = [0 0 1 0 0 0 0 0;0 0 0 1 0 0 0 0;0 0 -k 0 1 0 1 0;0 0 -0 -k 0 1 0 1;...
    0 0 0 0 -1/tau 0 0 0;0 0 0 0 0 -1/tau 0 0;0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0];

B = [0 0;0 0;0 0;0 0;1/tau 0;0 1/tau;0 0;0 0];

Aest = A;
DA = (A-Aest)*delta; % Used when there is a model error
A = eye(size(A))+delta*A;
Aest = eye(size(Aest))+delta*Aest;
B = delta*B;

% Observability Matrix
H = eye(size(A,1));
E = eye(8,1)';          %See Basar and Tamer, pp. 171

% Definition of the cost matrices:
Q = zeros(size(A,1),size(A,2),nStep);
M = Q;
TM = Q;
Id = eye(statedim);

%Filling in the cost matrices
for j = 1:nStep
    for i = 1:statedim
        
        Q(:,:,j) = Q(:,:,j) + ralpha(i,j)*Id(:,i)*Id(:,i)';
        
    end
end

% Augment the System to Include the Feedback Delay
A0 = A;
DA0 = DA;
Aest0 = Aest;
B0 = B;
Q0 = Q;
H0 = H;
[A,DA,B,Q,H] = AugRobustControl(A0,DA0,B0,Q0,H0,delay,delta);
[Aest,~,~,~,~] = AugRobustControl(Aest0,DA0,B0,Q0,H0,delay,delta);


%Signal Dependent Noise
nc = size(B,2);
Csdn = zeros(size(B,1),nc,nc);
for i = 1:nc
    
    Csdn(:,i,i) = .1*B(:,i);
    
end

M = Q;
TM = Q;
D = zeros(size(A));
D(1:8,1:8) = eye(8);

%--------------------------------------------------------------------------
% Implementing the backwards recursions
M(:,:,end) = Q(:,:,end);
SLQG = M;
L = zeros(size(B,2),size(A,1),nStep-1);  % Optimal Minimax Gains
C = L;                                   % Optimal LQG Gains
Lambda = zeros(size(A,1),size(A,2),nStep-1);

% Optimization of gamma
minlambda = zeros(nStep-1,1);
gammaK = 0.5;
reduceStep = 1;
positive = false;
relGamma = 1;

% Does the routine have to compute gamma opt?
if simdata.gamma(2)
    
    while (relGamma > .001 || ~positive)
        
        for k = nStep-1:-1:1
            
            % Minimax Feedback Control
            TM(:,:,k) = gamma^2*eye(size(A))-D'*M(:,:,k+1)*D;
            minlambda(k) = min(eig(TM(:,:,k)));
            
            Lambda(:,:,k) = eye(size(Aest))+(B*B'-gamma^-2*(D*D'))*M(:,:,k+1);
            M(:,:,k) = Q(:,:,k)+Aest'*(M(:,:,k+1)^-1+B*B'-gamma^-2*D*D')^-1*Aest;
            L(:,:,k) = B'*M(:,:,k+1)*Lambda(:,:,k)^-1*Aest;

        end
        
        oldGamma = gamma;
        
        if min(real(minlambda)) >= 0
            
            gamma = (1-gammaK)*gamma;
            relGamma = (oldGamma-gamma)/oldGamma;
            positive = true;
            
        elseif min(real(minlambda)) < 0
            
            gamma = (1-gammaK)^-1*gamma;
            reduceStep = reduceStep + 0.5;
            relGamma = -(oldGamma-gamma)/oldGamma;
            gammaK = gammaK^reduceStep;
            positive = false;
            
        end
        
    end
    
    gamma = oldGamma;
    
elseif ~simdata.gamma(2)
    
    for k = nStep-1:-1:1
        
        % Minimax Feedback Control
        TM(:,:,k) = gamma^2*eye(size(A))-D'*M(:,:,k+1)*D;
        minlambda(k) = min(eig(TM(:,:,k)));
        
        Lambda(:,:,k) = eye(size(Aest))+(B*B'-gamma^-2*(D*D'))*M(:,:,k+1);
        
        M(:,:,k) = Q(:,:,k)+Aest'*(M(:,:,k+1)^-1+B*B'-gamma^-2*D*D')^-1*Aest;
        L(:,:,k) = B'*M(:,:,k+1)*Lambda(:,:,k)^-1*Aest;

    end
   
end

%--------------------------------------------------------------------------
statedim = size(A,1);

%Forward Simulation of the System Trajectory
h = max(0,floor(delay/delta))+1;
currentX = kron(ones(h,1),xinit);
currentXEst = currentX;
x = zeros(nStep,statedim);
xest = x;

x(1,:) = currentX(1:statedim)';
xest(1,:) = currentX(1:statedim)';
u = zeros(nStep-1,size(B,2)); % size(B,2) is the control dimension
w = zeros(size(currentX,1),1);
% isForceFieldON = simdata.forcefield; 

% Parallel Simulation for LQG control
currentZ = currentX;
currentZEst = currentZ;
z = x;
zest = z;
v = u;
Oxi = 0.001*B*B';
Oxi(7:8,7:8) = Oxi(5:6,5:6);
Omega = eye(8)*Oxi(5,5)*simdata.noise(2);

%Parameters for State Estimation
Sigma = zeros(statedim,statedim,nStep);
Sigma(:,:,1) = eye(statedim)*10^-2;
SigmaK = Sigma;

pertX = false;
pertZ = false;

%--------------------------------------------------------------------------
% Extended LQG 
RLQG = zeros(2,2,nStep-1);
for i = 1:nStep-1
    RLQG(:,:,i) = eye(2);
end
Cstate = eye(statedim)*10^-2;

[C,Ke,~,~,~,~,~] = extLQG(Aest,B,Csdn,0*H,H,Q,RLQG,Oxi,Omega,0*A,Cstate,Cstate,currentZ);

%--------------------------------------------------------------------------

% Compute the total cost
cost = zeros(1,2);

for i = 1:nStep-1
    
    if i == 14 % Time step correspondng to ~1/3 of the reach path
        
        currentX(7:8) = simdata.pert;
        currentZ(7:8) = simdata.pert;
        
    end
    
    
    sensoryNoise = mvnrnd(zeros(size(Omega,1),1),Omega)';
    motorNoise = mvnrnd(zeros(size(Oxi,1),1),Oxi)';
    motorNoise(7:8) = zeros(2,1);
    
    %MINMAX HINFTY CONTROL ------------------------------------------------
    %Riccati Equation for the State Estimator
    Sigma(:,:,i+1) = Aest*(Sigma(:,:,i)^-1+H'*(E*E')^-1*H-gamma^-2*Q(:,:,i))^-1*Aest'+D*D';
    
    %Feedback Eequation
    yx = H*currentX + sensoryNoise;
    
    %Minmax Simulation with State Estimator
    u(i,:) = -B'*(M(:,:,i+1)^-1+B*B'-gamma^-2*(D*D'))^-1*Aest*...   %Control
        (eye(statedim)-gamma^-2*Sigma(:,:,i)*M(:,:,k))^-1*currentXEst;
    
    % Cost Hinf
    cost(1) = cost(1) + currentX'*Q(:,:,i)*currentX + u(i,:)*u(i,:)';
    
    %Signal Dependent Noise - Robust Control
    sdn = 0;
    
    for isdn = 1:nc
        sdn = sdn + normrnd(0,1)*Csdn(:,:,isdn)*u(i,:)';
    end
    
    %     u(i,:) = u(i,:)*(1-pertX)
    currentXEst = Aest*currentXEst + B*u(i,:)'+...
        Aest*(Sigma(:,:,i)^-1+H'*(E*E')^-1*H-gamma^-2*Q(:,:,i))^-1*(gamma^-2*Q(:,:,i)*currentXEst+H'*(E*E')^-1*(yx-H*currentXEst));
    
    % Minmax Simulation
    wx = DA*currentX; % Non zero if there is a model error. 
    currentX = Aest*currentX + B*u(i,:)'+D*wx + motorNoise + sdn;
    x(i+1,:) = currentX(1:statedim)';
    xest(i+1,:) = currentXEst(1:statedim)';
    
    %LQG CONTROL ----------------------------------------------------------
    yz = H*currentZ + sensoryNoise;
    v(i,:) = (-C(:,:,i)*currentZEst)';
    K = Ke(:,:,i);
    currentZEst = Aest*currentZEst + B*v(i,:)' + K*(yz-H*currentZEst);
    
    % Cost LQG
    cost(2) = cost(2) + currentZ'*Q(:,:,i)*currentZ + v(i,:)*v(i,:)';
    
    %Signal Dependent Noise - LQG
    sdn = 0;
    for isdn = 1:nc
        sdn = sdn + normrnd(0,1)*Csdn(:,:,isdn)*v(i,:)';
    end
    
    wz = DA*currentZ;
    currentZ = Aest*currentZ + B*v(i,:)' + D*wz + motorNoise + sdn;
    z(i+1,:) = currentZ(1:statedim)';
    zest(i+1,:) = currentZEst(1:statedim)';
    
end

% Add the final cost
cost(1) = cost(1) + currentX'*Q(:,:,end)*currentX;
cost(2) = cost(2) + currentZ'*Q(:,:,end)*currentZ;

% Output structure
sout.L = L;                      % Series of optimal robust control gains
sout.C = C;                      % Series of optimal LQG control gains   
sout.x = x(:,1:8);               % State - Robust control
sout.xest = xest(:,1:8);         % State Estimate - Robust control   
sout.z = z(:,1:8);               % State - LQG
sout.zest = zest(:,1:8);         % State Estimate - LQG
sout.u = u;                      % Series of Control Vector - Robust design
sout.v = v;                      % Series of Control Vector - LQG
sout.minlambda = minlambda;      % Min eigen value (optimized or used). Must be >0.
sout.cost = cost;                % Movement cost (1: Robust, 2: LQG)   
sout.gammaopt = gamma;           % Optimal or used gamma parameter    


