function [Data Measurements] = LightSim(Cont, Par, IJ,finit)

% Transform Parameters from half-lives to rates
Par(1) = log(2)/Par(1); % signal degredation rate
Par(2) = Par(2); % mRNA production rate (Exp,1)
Par(3) = log(2)/Par(3); % mRNA degredation rate
Par(4) = log(2)/Par(4); % protein degredation rate
Par(5) = log(2)/Par(5); % protein maturation rate

% Initial conditions
x = [0;1;1;1];

Acont = [[-Par(1) 0 0 0];[Par(2) -Par(3) 0 0];[0 Par(4)+Par(5) -Par(4)-Par(5) 0]; [0 0 Par(4) -Par(4)]];

% Discrete Time Dynamics (square is for 10 min sampling)
dt = 5; % minutes
Bcont = [0;1;0;0];
Acheat = expm([Acont Bcont; zeros(1,5)]*dt);
A = Acheat(1:4,1:4);
B = Acheat(1:4,5)*Par(3);

% Simulate
for i = 1:length(Cont)
    if Cont(i) == 1
            x(1,i) = 1;
    elseif Cont(i) == 2
        x(1,i) = 0;
    end
    x(:,i+1) = A*x(:,i) + B;
end

Data = x(1:4,:);

Measurements = [x(4,1:6:end)];
