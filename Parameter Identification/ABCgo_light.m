function [output] = ABCgo_light(ABCstruc)
%ABCGO Performs the ABC SMC algorithm of Toni, et. al. for parameter
%identification and model selection

% Start ABC SMC Algorithm
resample = 0;
Iters = ABCstruc.Iters;

for i = 1:Iters
    
    % Set epsilon
    if i == 1
        epsilon(i) = ABCstruc.Epsilon;
    elseif resample == 1;
        epsilon(i) = epsilon(i-1);
    else
        epsilon(i) = epsilon_new;
    end
    
    % Execution differs depending on the population indicator i
    if i == 1
        
        for j = 1:ABCstruc.Particles
            
            if mod(j,10) == 0
                j
            end
            
            notaccepted = 1;
            while notaccepted
                
                clear m p Data meas
                
                % Sample model from uniform prior
                m = ceil(rand*ABCstruc.Models);
                
                % Sample parameters from uniform prior
                if ~isempty(ABCstruc.FreeParameters{m})
                    p = ABCstruc.DownParameters{m} + rand(1,length(ABCstruc.FreeParameters{m})).*(ABCstruc.UpParameters{m}-ABCstruc.DownParameters{m});
                end
                
                % Simulate candidate dataset
                Data = [];
                Meas = [];
                for ij = 1:length(ABCstruc.Fit)
                    clear temp1 temp2
                    [temp1 temp2] = LightSim(ABCstruc.control{ij}, p, ij,ABCstruc.Fit{ij}(1));
                    Data = [Data;transpose(temp2(:,1:length(ABCstruc.Fit{ij})))];
                    Meas = [Meas;transpose(ABCstruc.Fit{ij})];
                    
                end
                
                Deps = sqrt(sum((Data(:)-Meas(:)).^2));
                %pause
                
                % Accept/Reject
                if Deps < epsilon(i)
                    
                    
                    % assign a weight to the particle
                    w = 1; % possibly altered from simulation
                    
                    % save particle
                    Particle{j}.m = m;
                    Particle{j}.p = p;
                    %Particle{j}.ic = ic(ABCstruc.FreeIC{m});
                    Particle{j}.w = w;
                    
                    % break from while loop
                    break
                    
                end
                
            end
            
        end
        
    else
        
        for j = 1:ABCstruc.Particles
            
            if mod(j,10) == 0
                j
            end
            
            notaccepted = 1;
            while notaccepted
                
                clear m1 m p1 p Data meas
                
                % Sample model from proposal distribution
                m = 1;
                
                % Sample parameters/ICs from proposal parameter distribution
                randp = randsrc(1,1,[params{m}.int;params{m}.dist]);              
                
                % Sample new parameter set according to perterbation kernel
                p = Particle{randp}.p + ABCstruc.KernelParameters{m}.*randn(1,length(ABCstruc.FreeParameters{m}));                
                
                % Check parameter bounds
                while min(p>=ABCstruc.DownParameters{m})==0 || min(ABCstruc.UpParameters{m}>=p)==0
                    p = Particle{randp}.p + ABCstruc.KernelParameters{m}.*randn(1,length(ABCstruc.FreeParameters{m}));
                end                
                
                % Simulate candidate dataset
                Data = [];
                Meas = [];
                for ij = 1:length(ABCstruc.Fit)
                    clear temp1 temp2
                    [temp1 temp2] = LightSim(ABCstruc.control{ij}, p, ij,ABCstruc.Fit{ij}(1));
                    Data = [Data;transpose(temp2(:,1:length(ABCstruc.Fit{ij})))];
                    Meas = [Meas;transpose(ABCstruc.Fit{ij})];
                end
                
                Deps = sqrt(sum((Data(:)-Meas(:)).^2));
                
                % Accept/Reject
                if Deps < epsilon(i)
                    
                    % assign a weight to the particle
                    if resample == 1
                        w = 1;
                    else
                        S1 = 0;
                        for k = 1:ABCstruc.Models
                            S1 = S1 + model(k)*ABCstruc.KernelModel(k,m);
                        end
                        S2 = 0;
                        for k = 1:length(params{m}.int)
                            tempK = 1;
                            for h = 1:length(ABCstruc.KernelParameters{m})
                                tempK = tempK*ABCstruc.Kernel(p(h),Particle{params{m}.int(k)}.p(h),ABCstruc.KernelParameters{m}(h));
                            end
%                             if ~isempty(ABCstruc.IC{m})
%                                 if ~isempty(ABCstruc.FreeIC{m})
%                                     for h = 1:length(ABCstruc.KernelIC{m})
%                                         tempK = tempK*ABCstruc.Kernel(ic(ABCstruc.FreeIC{m}(h)),Particle{params{m}.int(k)}.ic(h),ABCstruc.KernelIC{m}(h));
%                                     end
%                                 end
%                             end
                            S2 = S2 + Particle{params{m}.int(k)}.w/model(m)*tempK;
                        end
                        S = S1*S2;
                        % obtain b by M
                        M = 1; b = zeros(M,1);
                        if M > 1
                            for ai = 1:M
                                measb = Data.statevalues(2:length(ABCstruc.TimePoints)+1,ABCstruc.MeasStates) + ...
                                    [ones(length(ABCstruc.TimePoints),1)*ABCstruc.SimNoise].*randn(length(ABCstruc.TimePoints),length(ABCstruc.MeasStates));
                                Depsb = ABCstruc.Dist(measb(:),Measurements(:));
                                if Depsb < epsilon(i)
                                    b(ai) = 1;
                                end
                            end
                        else
                            b = 1;
                        end
                
                        w = mean(b)/prod([ABCstruc.UpParameters{m}]-[ABCstruc.DownParameters{m}])/S; % possibly altered from simulation
                    end
                    
                    % save particle
                    ParticleNew{j}.m = m;
                    ParticleNew{j}.p = p;
                    %ParticleNew{j}.ic = ic(ABCstruc.FreeIC{m});
                    ParticleNew{j}.w = w;
                    
                                        
                    % break from while loop
                    break
                    
                end
                
            end
            
        end
        
        Particle = ParticleNew;
        
    end
    
    % Build Marginal Posterior Distributions
    clear model params
    
    % Marginal Posterior for the model
    model = zeros(1,ABCstruc.Models);
    for k = 1:ABCstruc.Particles
        model(Particle{k}.m) = model(Particle{k}.m) + Particle{k}.w;
    end
    % Normalize the weights
    for k = 1:ABCstruc.Particles
        Particle{k}.w = Particle{k}.w/sum(model);
    end
    model = model/sum(model);
    % Compute the effective sample size
    ESS(i) = 0; maxw = 0; ESS2(i) = 0;
    for k = 1:ABCstruc.Particles
        ESS(i) = ESS(i) + (ABCstruc.Particles*Particle{k}.w-1)^2;
        ESS2(i) = ESS2(i) + Particle{k}.w^2;
        maxw = max(maxw,Particle{k}.w);
    end
    maxw = maxw
    ESS(i) = ABCstruc.Particles/(1 + (1/ABCstruc.Particles)*ESS(i))
    ESS2(i) = 1/ESS2(i)
    if ESS(i) < ABCstruc.Resample
        resample = 1;
    else
        resample = 0;
    end
    if resample == 0
        % Adaptively compute the next epsilon
        for k = 1:ABCstruc.Particles
            m = Particle{k}.m;
            %ic(ABCstruc.FreeIC{m}) = Particle{k}.ic;
            % Simulate candidate dataset
            %Data = SBPDsimulate(ABCstruc.Name,[0 ABCstruc.TimePoints],ic,...
            %   [ABCstruc.SetParameters{m},ABCstruc.FreeParameters{m}],[ABCstruc.ValParameters{m},Particle{k}.p]);
            % Simulate candidate dataset
                Data = [];
                Meas = [];
                for ij = 1:length(ABCstruc.Fit)
                    clear temp1 temp2
                    [temp1 temp2] = LightSim(ABCstruc.control{ij}, Particle{k}.p, ij,ABCstruc.Fit{ij}(1));
                    Data = [Data;transpose(temp2(:,1:length(ABCstruc.Fit{ij})))];
                    Meas = [Meas;transpose(ABCstruc.Fit{ij})];
                end
                PA(k) = sqrt(sum((Data(:)-Meas(:)).^2));
        end
        alive = find(PA<epsilon(i));
        sortalive = sort(PA(alive));
        targetalive = ceil(ABCstruc.EpsAccept*length(sortalive));
        if ~isempty(targetalive) && targetalive ~=0
        epsilon_new = sortalive(targetalive)
        else
            %break
        end
    end
    
    % Marginal Posterior for the parameters
    for k = 1:ABCstruc.Models
        params{k}.dist = [];
        params{k}.int = [];
    end
    for k = 1:ABCstruc.Particles
        params{Particle{k}.m}.dist = [params{Particle{k}.m}.dist, Particle{k}.w];
        params{Particle{k}.m}.int = [params{Particle{k}.m}.int, k];
    end
    for k = 1:ABCstruc.Models
        params{k}.dist = params{k}.dist/sum(params{k}.dist);
        numpartsinmodel(k) = length(params{k}.dist);
    end
    numpartsinmodel = numpartsinmodel
    
    [t1 t2] = min(PA)
    p = Particle{t2}.p
    % plot
    close all
    b = ['b','k','c','g','k','c','g','k','c','g','r','k','m--','y--','k:-','r:-'];
    for ij = 1:length(ABCstruc.Fit)
        for ddd = 1:2
        Data = [];
        Meas = [];
        figure(ddd)
        hold on
        clear temp1 temp2
        [Data Meas] = LightSim(ABCstruc.control{ij}, p, ij,ABCstruc.Fit{ij}(1));
        %subplot(2,1,1)
        hold on
        length(Data)
        if ddd == 1
        plot(0:5:(length(Data)-1)*5,Data(4,:))
        plot(0:30:(length(ABCstruc.Fit{ij})-1)*30,ABCstruc.Fit{ij}(1,:),'-r*')
        %subplot(2,1,2)
        hold on
        plot(0:30:(length(ABCstruc.Fit{ij})-1)*30,Meas(1,1:length(ABCstruc.Fit{ij})),'--bo')
        plot(0:30:(length(ABCstruc.Fit{ij})-1)*30,ABCstruc.Fit{ij}(1,:),'--r*')
        else
            plot(0:5:(length(Data)-1)*5,Data(4,:)*0.0192)
        plot(0:30:(length(ABCstruc.Fit{ij})-1)*30,ABCstruc.Fit{ij}(1,:)*0.0192,'-r*')
        %subplot(2,1,2)
        hold on
        plot(0:30:(length(ABCstruc.Fit{ij})-1)*30,Meas(1,1:length(ABCstruc.Fit{ij}))*0.0192,'--bo')
        plot(0:30:(length(ABCstruc.Fit{ij})-1)*30,ABCstruc.Fit{ij}(1,:)*0.0192,'--r*')
        end
        %pause
        end
    end
    %close all
    %end
    pause(1)
close all

    if epsilon(i) < ABCstruc.Breakpoint && ESS(i) > ABCstruc.Resample
        break
    end
end

    % Visualize Parameters
    for modd = 1:ABCstruc.Models
        if isempty(params{modd}.int) == 0
            figure
            for k = 1:length(ABCstruc.FreeParameters{modd})
                subplot(1,length(ABCstruc.FreeParameters{modd}),k)
                for h = 1:1000
                    ptest(h,1) = Particle{randsrc(1,1,[params{modd}.int;params{modd}.dist])}.p(k);
                    Distribute(h,:) = Particle{randsrc(1,1,[params{modd}.int;params{modd}.dist])}.p;
                end
                hist(ptest,30)
%                xlabel(sprintf('%s (%d, %d)',ABCstruc.FreeParameters{modd}{k},ABCstruc.DownParameters{modd}(k),ABCstruc.UpParameters{modd}(k)));
                ylabel('samples')
                Distribute(:,k) = ptest;
            end
        end
        
    end
    
    % Visualize tolerance and effective sample size
    figure
    plot(ESS)
    figure
    plot(epsilon)

output.Particle = Particle;
output.Dist = Distribute;
