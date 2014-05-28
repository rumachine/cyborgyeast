% Initialize the system
function [opendata closedata] = CTL_insilico(narg)

load('OpenLoop7')

% How long is the experiment
Experiment_Length = 420;
dt = 5;

% Parameter set 1
Par = [44.79 0.9587 23.09 105.53 16.55];

% Transform Parameters from half-lives to rates
Par(1) = log(2)/Par(1); % signal degredation rate
Par(2) = Par(2); % mRNA production rate (Exp,1)
Par(3) = log(2)/Par(3); % mRNA degredation rate
Par(4) = log(2)/Par(4); % protein degredation rate
Par(5) = log(2)/Par(5); % protein maturation rate

% Continuous System Dynamics
Acont = [[-Par(1) 0 0 0];[Par(2) -Par(3) 0 0];[0 Par(4)+Par(5) -Par(4)-Par(5) 0]; [0 0 Par(4) -Par(4)]];
Bcont = [0;1;0;0];

% Discrete Time Dynamics (square is for 10 min sampling)
Acheat = expm([Acont Bcont; zeros(1,5)]*dt);
A1 = Acheat(1:4,1:4);
B1 = Acheat(1:4,5)*Par(3);

% Parameter set 2 (Parameter change)
Par = Par;%+ randn(1,7)*sqrtm(covp);

% Continuous System Dynamics
Acont = [[-Par(1) 0 0 0];[Par(2) -Par(3) 0 0];[0 Par(4)+Par(5) -Par(4)-Par(5) 0]; [0 0 Par(4) -Par(4)]];
Bcont = [0;1;0;0];

% Discrete Time Dynamics (square is for 10 min sampling)
Acheat = expm([Acont Bcont; zeros(1,5)]*dt);
A2 = Acheat(1:4,1:4);
B2 = Acheat(1:4,5)*Par(3);

% Initialize State (account for delay)
x{1} = [0;1;1;1];

% Troubleshoot
Amodel = A1;
Bmodel = B1;
Atrue = A2;
Btrue = B2;

% Output matrix
C = [zeros(1,3), 1];
rank(obsv(Amodel,C));

% Target set points
xtarget{1} = 7; % Can Define multiple target values for the system

% Model Sampling time is nominally 5 minutes
step = dt; % Sampling time in minutes
time(1) = 0; % Initialize the time aspect
ccycle = 3; % Control possible implemented every X cycles
ecycle = 6; % estimate is updated every X cycles (delayed by 1 cycle)
mcount = 1; % measurement count
mtime(mcount) = 0;
meandata(1:length(xtarget),1) = 1;
opendata(1:length(Open(:,1)),1) = 1;
Iter = Experiment_Length/step;

% Control Objective Variables
N = 6; % Control horizon in terms of control cycles
%W = eye(N*ccycle+1);
W = diag(1:(N*ccycle+1));
Obj = @(state,target) (state-target)*W*(state-target)'; % Euclidean distance as cost
[q1 q2 q3 q4 q5 q6] = ndgrid(0:2,0:2,0:2,0:2,0:2,0:2); % Enumeration is fast enough
%[q1 q2 q3 q4] = ndgrid(0:4,0:4,0:4,0:4);
Cont = [q1(:),q2(:),q3(:),q4(:),q5(:),q6(:)];
%Cont = [q1(:),q2(:),q3(:),q4(:)];

% State Estimation Variables
P{1} = zeros(length(Amodel)); % initial covariance matrix
R = 0.05; % measurement covariance
Q = diag([0.1 10 10 0.1]); % system process noise

% Initialize State Estimation
x_est{1} = [0;1;1;1];

%%%%%%%%%%%% End Initialization %%%%%%%%%%%%
for sys = 1:length(xtarget)
            control{sys} = zeros(1,Experiment_Length/dt);
end
% Begin Loop
for i = 1:Iter
    
    time(i+1) = i*step;
    
    %%%%%%%%% Control Step %%%%%%%%%%%%
    if mod(time(i), step*ccycle) == 0
        %pause(0.01)
        % Given the current state of the system, what is the open loop
        % control process over the next 15 (30) minutes?
        % control variables are 0, 1, 2
        % 0 - do nothing, 1 - 660 nm, 2 - 730 nm
        for sys = 1:length(xtarget)
            tcost = inf;
            for k = 1:length(Cont(:,1))
                clear xtemp ctemp
                xtemp(:,1) = x_est{sys}(:,i);
                for j = 1:N
                    % set signal state according to control sequence tested
                    if Cont(k,j) == 1 %|| Cont(k,j)==3 || Cont(k,j)==4
                        xtemp(1,(j-1)*ccycle+1) = 1;
                    elseif Cont(k,j) == 2
                        xtemp(1,(j-1)*ccycle+1) = 0;
                    end
                    % calculate system state according to signal states
                    % defined by the control
                    for n = 1:ccycle
                        xtemp(:,(j-1)*ccycle+n+1) = Amodel*xtemp(:,(j-1)*ccycle+n) + Bmodel;
                    end
                end
                ctemp = Obj(xtemp(4,:),xtarget{sys});
                if ctemp < tcost
                    bestq = k;
                    tcost = ctemp;
                end
            end
            control{sys}(i) = Cont(bestq,1);
            if control{sys}(i) == 0
            elseif control{sys}(i) == 1
            elseif control{sys}(i) == 2
                if x_est{sys}(1,i-ccycle) == 0
                    control{sys}(i) = 0;
                end
            end
        end            
    end

    
    %%%%%%%%% Measurement Step %%%%%%%%%
    % First check if we need to update the estimate with a new available
    % measurement
    meas_up = 0;
    if i > 1 && mod(time(i),step*ecycle) == 0
        clear thisdata thisdata2
        clc
        cont_on = 0;
        while cont_on == 0
            % Obtain the measurement from cytometer data
                     for sys = 1:length(xtarget)
                        thisdata(sys,1) = C*x{sys}(:,i);% + sqrt(R)*randn;
                     end

                break
           
        end
        
        % Obtain measurement data from cytometer
        meandata(:,mcount+1) = thisdata;
        mtime(mcount+1) = time(i);
        
        % note that an update needs to be made
        meas_up = 1;
        
        % progress the measurement count
        mcount = mcount + 1;
        
    end
    
    % Progress the state estimate according to the control estimate
    for sys = 1:length(xtarget)
        
        % Update the state of the system for the step before (incorporate
        % control input into the state update)
        P{sys} = Amodel*P{sys}*Amodel' + Q;
        S = C*P{sys}*C' + R;
        K = meas_up*P{sys}*C'*(S^-1);
        P{sys} = P{sys} - (K*C)*P{sys};
        
        % Update the current state estimate
        if meas_up == 1
            x_est{sys}(:,i) = Amodel*x_est{sys}(:,i-1) + Bmodel + K*(meandata(sys,mcount) - C*x_est{sys}(:,i));
        end
        
    end
    
    % Progress true and approximate state according to new control input
    for sys = 1:length(xtarget)
        U{sys} = control{sys}(i);
        % update state estimate
        if U{sys} == 1
                x{sys}(1,i) = 1;
                x_est{sys}(1,i) = 1;
        elseif U{sys} == 2
            x{sys}(1,i) = 0;
            x_est{sys}(1,i) = 0;
        end
        x{sys}(:,i+1) = Atrue*x{sys}(:,i) + Btrue;
        x_est{sys}(:,i+1) = Amodel*x_est{sys}(:,i) + Bmodel;
        
    end
    
    
end
for op = 1:length(Open(:,1))
    x_open{op} = x{1}(:,1);
    for i = 1:Iter
        if Open(op,i) == 1
                x_open{op}(1,i) = 1;
        elseif Open(op,i) == 2
            x_open{op}(1,i) = 0;
        end
        x_open{op}(:,i+1) = Atrue*x_open{op}(:,i) + Btrue;
    end
end
for sys = 1:length(xtarget)
    figure(1)
    subplot(4,1,1)
    hold on
    plot(time(1:i+1), ones(size(time(1:i+1)))*xtarget{sys},'b--')
    plot(time(1:i+1), x_open{op}(4,:),'r')
    plot(time(1:6:length(time)), x{sys}(4,1:6:length(time)),'b')
    plot(time(1:6:length(time)), x_est{sys}(4,1:6:length(time)),'c--')
    subplot(4,1,2)
    hold on
    plot(time(1:6:length(time)), x{sys}(3,1:6:length(time)),'b')
    plot(time(1:6:length(time)), x_est{sys}(3,1:6:length(time)),'c--')
    subplot(4,1,3)
    hold on
    plot(time(1:6:length(time)), x{sys}(2,1:6:length(time)),'b')
    plot(time(1:6:length(time)), x_est{sys}(2,1:6:length(time)),'c--')
    subplot(4,1,4)
    hold on
    plot(time(1:6:length(time)), x{sys}(1,1:6:length(time)),'b')
    plot(time(1:6:length(time)), x_est{sys}(1,1:6:length(time)),'c--')
    figure(2)
    hold on
    plot(time(1:6:length(time)), x{sys}(4,1:6:length(time)),'x')
    %plot(time(1:6:length(time)), x_est{sys}(4,1:6:length(time)),'c--')
    plot(time(1:i+1), ones(size(time(1:i+1)))*xtarget{sys},'k--')
    plot(time(1:i+1), x_open{op}(4,:),'b')
end
opendata = x_open{sys}(4,1:6:length(time));
closedata = x{sys}(4,1:6:length(time));
figure
stem(0:5:415,control{1})
% pause
% Open = control{1};
% save newopen2 Open