% Initialize the system
clear all

% How long is the experiment
Experiment_Length = 420;

% Parameter set
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

% Initialize State (account for delay)
x{1} = [0;1;1;1];
x{2} = [0;1;1;1];

% Discrete Time Dynamics (square is for 10 min sampling)
dt = 5;
Acheat = expm([Acont Bcont; zeros(1,5)]*dt);
A1 = Acheat(1:4,1:4);
B1 = Acheat(1:4,5)*Par(3);

% Troubleshoot
Amodel = A1;
Bmodel = B1;

% Output matrix
C = [zeros(1,3), 1];
rank(obsv(Amodel,C));

% Target set points
xtarget{1} = 4; % Can Define multiple target values for the system
xtarget{2} = 4;

% Model Sampling time is nominally 5 minutes
step = dt; % Sampling time in minutes
time(1) = 0; % Initialize the time aspect
ccycle = 3; % Control possible implemented every X cycles
ecycle = 6; % estimate is updated every X cycles (delayed by 1 cycle)
mcount = 0; % measurement count
meandata(1:length(xtarget),1) = 1;
Iter = Experiment_Length/step;

% Control Objective Variables
N = 6; % Control horizon in terms of control cycles
W{1} = diag(1:(N*ccycle+1));
W{2} = diag(1:(N*ccycle+1));
Obj = @(state,target,X) (state-target)*X*(state-target)'; % Euclidean distance as cost
[q1 q2 q3 q4 q5 q6] = ndgrid(0:2,0:2,0:2,0:2,0:2,0:2);
Cont = [q1(:),q2(:),q3(:),q4(:),q5(:),q6(:)];

% State Estimation Variables
P{1} = zeros(length(Amodel)); % initial covariance matrix
P{2} = zeros(length(Amodel)); % initial covariance matrix
R{1} = 0.05; % measurement covariance
R{2} = 0.05; % measurement covariance
Q{1} = diag([0.1 10 10 0.1]); % system process noise
Q{2} = diag([0.1 10 10 0.1]); % system process noise

% Initialize State Estimation
x_est{1} = [0;1;1;1];
x_est{2} = [0;1;1;1];

%%%%%%%%%%%% End Initialization %%%%%%%%%%%%
for sys = 1:length(xtarget)
    control{sys} = zeros(1,Experiment_Length/dt);
end
% Begin Loop
for i = 1:Iter
    
    time(i+1) = i*step;
    
    %%%%%%%%% Control Step %%%%%%%%%%%%
    if mod(time(i), step*ccycle) == 0
        pause(0.01)
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
                ctemp = Obj(xtemp(4,:),xtarget{sys},W{sys});
                if ctemp < tcost
                    bestq = k;
                    tcost = ctemp;
                end
            end
            control{sys}(i) = Cont(bestq,1);
            if control{sys}(i) == 0
                %disp(sprintf('Sample %d:  At %d minutes do nothing',sys,time(i)))
            elseif control{sys}(i) == 1
                disp(sprintf('Sample %d:  At %d minutes apply R for one minute',sys,time(i)))
            elseif control{sys}(i) == 2
                clear tempO
                tempO = find(control{sys}(1:i-1)~=0);
                if control{sys}(tempO(end))==2
                    %disp(sprintf('Sample %d:  At %d minutes do nothing',sys,time(i)))
                    control{sys}(i) = 2;
                else
                    disp(sprintf('Sample %d:  At %d minutes FR for one minute',sys,time(i)))
                end
            end
        end
    end
    
    %%%%%%%%% Measurement Step %%%%%%%%%
    % First check if we need to update the estimate with a new available
    % measurement
    meas_up = 0;
    if mod(time(i),step*ecycle) == 0
        
        % progress the measurement count
        mcount = mcount + 1;
        
        clear thisdata thisdata2
        % Display press a button
        disp('Press any key when FCS files have been uploaded')
        pause
        close all
        cont_on = 0;
        while cont_on == 0
            % Obtain the measurement from cytometer data
            thisdata2 = RunNowClosed(mcount,length(xtarget));
            
            if isempty(thisdata2)
                disp('The new measurements are not in the folder')
                disp('Please upload the new measurements and press any key')
                pause
            else
                thisdata = thisdata2(1:length(xtarget));
                clc
                break
            end
            
        end
        
        % Obtain measurement data from cytometer
        meandata(:,mcount) = thisdata;
        mtime(mcount) = time(i);
        
        % note that an update needs to be made
        meas_up = 1;
        
    end
    close all
    % Progress the state estimate according to the control estimate
    for sys = 1:length(xtarget)
        
        % Update the state of the system for the step before (incorporate
        % control input into the state update)
        P{sys} = Amodel*P{sys}*Amodel' + Q{sys};
        S = C*P{sys}*C' + R{sys};
        K = meas_up*P{sys}*C'*(S^-1);
        P{sys} = P{sys} - (K*C)*P{sys};
        
        % Update the current state estimate
        if meas_up == 1
            if i == 1
                x_est{sys}(:,i) = x_est{sys}(:,1) + K*(meandata(sys,mcount) - C*x_est{sys}(:,i));
            else
                x_est{sys}(:,i) = Amodel*x_est{sys}(:,i-1) + Bmodel + K*(meandata(sys,mcount) - C*x_est{sys}(:,i));
            end
        end
        
    end
    
    % Progress true and approximate state according to new control input
    %close all
    ddd = [{'-b*'}, {'-b*'}, {'-b*'},{'-r*'}];
    for sys = 1:length(xtarget)
        U{sys} = control{sys}(i);
        %pause
        % update state estimate
        if U{sys} == 1
                x_est{sys}(1,i) = 1;
        elseif U{sys} == 2
            x_est{sys}(1,i) = 0;
        end
        x_est{sys}(:,i+1) = Amodel*x_est{sys}(:,i) + Bmodel;
        
        
        figure(1)
        for k = 1:4
            subplot(4,1,k)
            hold on
            if k == 4
                plot(mtime, meandata(sys,:),'-r*')
                plot(time(1:i+1), ones(size(time(1:i+1)))*xtarget{sys},'k--')
                plot(time(1:i+1), x_est{sys}(k,:),'b')
            else
                plot(time(1:i+1), x_est{sys}(k,:),'b')
            end
        end
    end
    
    
end

cont_on = 0;
while cont_on == 0
    % Obtain the measurement from cytometer data
    % thisdata2 = RunNowClosed(mcount+1,length(xtarget)+length(Open(:,1)));
    thisdata2 = RunNowClosed(mcount+1,length(xtarget));
    
    if isempty(thisdata2)
        disp('The new measurements are not in the folder')
        disp('Please upload the new measurements and press any key')
        pause
    else
        thisdata = thisdata2(1:length(xtarget));
        clc
        break
    end
    
end

% Obtain measurement data from cytometer
meandata(:,mcount+1) = thisdata;
mtime(mcount+1) = time(i+1);

% note that an update needs to be made
meas_up = 1;

% progress the measurement count
mcount = mcount + 1;
ddd = [{'-b*'}, {'-b*'}, {'-b*'},{'-r*'}];
for sys = 1:length(xtarget)
    U{sys} = control{sys}(i);
    %pause
    % update state estimate
    if U{sys} == 1
            x_est{sys}(1,i) = 1;
    elseif U{sys} == 2
        x_est{sys}(1,i) = 0;
    end
    x_est{sys}(:,i+1) = Amodel*x_est{sys}(:,i) + Bmodel;
    
    
    figure(1)
    for k = 1:4
        subplot(4,1,k)
        hold on
        if k == 4
            plot(mtime, meandata(sys,:),'-r*')
            plot(time(1:i+1), ones(size(time(1:i+1)))*xtarget{sys},'k--')
            plot(time(1:i+1), x_est{sys}(k,:),'b')
        else
            plot(time(1:i+1), x_est{sys}(k,:),'b')
        end
    end
end