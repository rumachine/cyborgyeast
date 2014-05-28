% Initialize the ABC parameters via ABCstruc structure
clear all
favg = 0.018;

load March23A
realtime1 = [mediandata];
shine1 = [control{1};control{2};control{3};control{4};control{5}];
clear control Open

load March23S
realtime2 = [mediandata];
shine2 = [control{1};control{2};control{3}];
clear control Open

for i = 1:5
    Fit{i} = [realtime1(i,:)/favg];
    control{i} = shine1(i,:);
end

for i = 6:8
   Fit{i} = [realtime2(i-5,:)/favg];
   control{i} = shine2(i-5,:);
end 

ABCstruc.Fit = Fit;
ABCstruc.control = control;

% Gaussian Noise Attributed to each Measured state
ABCstruc.MeasNoise = [0.01];
ABCstruc.SimNoise = ABCstruc.MeasNoise;

% Number of Models
ABCstruc.Models = 1;

ABCstruc.FreeParameters{1} = [1 2 3 4 5];

% Set Upper/Lower Parameter values for each model (cell array)
ABCstruc.DownParameters{1} = [10 0.5 15 100 10];
ABCstruc.UpParameters{1} = [60 5 30 140 20];

% ABC Algorithm Specifications

% Epsilon Schedule
ABCstruc.Iters = 10;
ABCstruc.Epsilon = 8;
ABCstruc.EpsAccept = 0.7;
ABCstruc.Resample = 20;
ABCstruc.Breakpoint = 0.1;

% Number of Particles
ABCstruc.Particles = 100;

% Model Perterbation Kernel (Finite State Discrete Time Markov Chain -
% Defined by P matrix)
ABCstruc.KernelModel = [1];

% Parameter Perterbation Kernel
ABCstruc.Kernel = @(x,y,sigma) (1/sqrt(2*pi*sigma^2))*exp(-(x-y)^2/2/sigma^2);
spread = 0.3;
for i = 1:1
    ABCstruc.KernelParameters{i} = spread*(ABCstruc.UpParameters{i}-ABCstruc.DownParameters{i});
end
