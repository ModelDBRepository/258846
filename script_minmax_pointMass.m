%   scritp_minmax_PointMass
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances. 
%   Crevecoeur F., Scott S. H., Cluff T. 
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019




% Define the data structure simdata, with parameters for the simulations:

simdata.delta = .01;        % Discretization step: 10ms
simdata.delay = .05;        % feedback loop delay, 5 time steps
simdata.pert = [0 0]';      % Perturbation magnitude, x and y coordinates, in N
simdata.time = 0.6;         % Reach time
simdata.gamma = [50000 1];  % First parameter is initial guess,
                            % Second parameter indicates whether it must be optimized
simdata.nStep = 61;         % Number of time steps corresponding to reach time (600ms), plus terminal step
% simdata.forcefield = 0;   Stay tuned
simdata.noise= [1 1];       % Sensory and motor noise, standard values.

% Populates the matrix runningalpha with the cost values:
runningalpha = zeros(8,simdata.nStep); 
for i = 1:simdata.nStep
    
    fact = min(1,(i*simdata.delta/simdata.time))^6;
    runningalpha(:,i) = [fact*10^6 fact*10^6 fact*10^5 fact*10^5 1 1 1 1]';
    
end
simdata.ralpha = runningalpha;

% Compute the optimal gamma
test = minmaxfc_pointMass([0 0 0 0 0 0 0 0]',[0 .15 0 0 0 0 0 0]',simdata);
simdata.gamma = [test.gammaopt, 0];


% Iterations
simdata.nsimu = 5; % Number of simulation runs.

costLQG = zeros(simdata.nsimu,1); % Extact movement costs
costHoo = zeros(simdata.nsimu,1);

maxLQG = zeros(simdata.nsimu,1);
maxHoo = zeros(simdata.nsimu,1);

avControlHoo = 0;
avControlLQG = 0;
averagePlot = 0;

% Normalization factor for the control variables
if simdata.pert(1) == 0
    normc = 1;
else
    normc = simdata.pert(1);
end

% Iterations
for i = 1:simdata.nsimu
    
    % Run the minmax control simulation
    test = minmaxfc_pointMass([0 0 0 0 0 0 0 0]',[0 .15 0 0 0 0 0 0]',simdata);
    
    ns = size(test.L,3);
    robustGain = zeros(1,ns);
    LQGGain = zeros(1,ns);
    
    for k = 1:ns
        robustGain(k) = norm(test.L(2,2,k));
        LQGGain(k) = norm(test.C(2,2,k));
    end
    
    % Puts hold on to add simulation to the figure
    subplot(221)
    plot(test.x(:,1),test.x(:,2),'r'), hold on, axis square;
    plot(test.z(:,1),test.z(:,2),'b')
    axis([-.1 .1 -.17 .03])
    
    % Average traces for plot
    averagePlot = averagePlot + [test.x(:,2)'+.15;test.z(:,2)'+.15;test.x(:,1)';test.z(:,1)';...
        test.x(:,4)'+.15;test.z(:,4)'+.15;test.x(:,3)';test.z(:,3)']/simdata.nsimu;
    
    % Sensitivity
    costLQG(i) = log10(test.cost(2));
    costHoo(i) = log10(test.cost(1));
    
    maxLQG(i) = max(abs(test.z(:,1)));
    maxHoo(i) = max(abs(test.x(:,1)));
    
    %Average traces for control
    % Control response is normalized to the perturbation amplitude, when
    % there is no perturbation (load = 0) the raw values are used. 
    
    avControlHoo = avControlHoo + test.u(:,1)/(abs(normc)*simdata.nsimu); 
    avControlLQG = avControlLQG + test.v(:,1)/(abs(normc)*simdata.nsimu);
    
end

subplot(221)
xlabel('x [cm]','FontSize',12);
ylabel('y [cm]','FontSize',12);

% Average forward velocity and lateral velocity
subplot(222)
plot(averagePlot(5,:),'r'), hold on
plot(averagePlot(6,:),'b');
plot(averagePlot(7,:),'r:'), hold on
plot(averagePlot(8,:),'b:');
axis square
xlabel('Time Steps','FontSize',12);
ylabel('Forward and Lateral Vel. [m/s]','FontSize',12);

% Average control responses
subplot(223)
plot(avControlHoo,'r'), hold on
plot(avControlLQG,'b');
axis square
xlabel('Time Steps','FontSize',12);
ylabel('\Delta Control [a.u.]','FontSize',12);
legend('Robust','LQG')

% Sensitivity 
subplot(224)
plot(mean(costLQG),mean(maxLQG),'bo','MarkerSize',8,'MarkerFaceColor','b'), hold on;
plot(mean(costHoo),mean(maxHoo),'ro','MarkerSize',8,'MarkerFaceColor','r');
axis square
xlabel('Movement cost (log)');
ylabel('Max lateral displacement');


