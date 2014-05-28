function [data logdata normlogdata normdata Fscatterdata Sscatterdata histcount x_axis gate_fr] = process_plate(plate,rows,cols,gate)
% process the cell array 'plate' in order to extract the data needed. 'rows' is the number of rows measured on the plate 
%(not necessarily consecutive) and 'cols' is the number of measured wells
% per row (consecutive)

%% Ellipse data for gating
A = [[5.4127   -1.9511];[-1.9511    6.2387]]; %ellipse given by x'Ax
off = [4.2026 4.1240]; %ellipse center (x1 = Fscatter, x2 = Sscater)
level = 1.5;
%%


data = cell(rows,cols);
logdata = cell(rows,cols);
normlogdata = cell(rows,cols);
normdata = cell(rows,cols);
histcount = cell(rows,cols);
x_axis = cell(rows,cols);
Fscatterdata = cell(rows,cols);
Sscatterdata = cell(rows,cols);
for i = 1:rows
    for j = 1:cols
    if(~isempty(plate{i,j}))      
    data{i,j} = plate{i,j}.DATA(:,3);         % select which channel to process (channel 3 --> FITC)
    Fscatterdata{i,j} = plate{i,j}.DATA(:,1);
    Sscatterdata{i,j} = plate{i,j}.DATA(:,2);
    cellsizedata =  plate{i,j}.DATA(:,1);
    zerodata = find(data{i,j}==0);
    zerodata = [zerodata;find(cellsizedata==0)];
    if(gate==1)
          N1 = numel(Fscatterdata{i,j});
          gatedata = find(diag([log10(Fscatterdata{i,j})-off(1) log10(Sscatterdata{i,j})-off(2)]*A*[log10(Fscatterdata{i,j})-off(1) log10(Sscatterdata{i,j})-off(2)]')>level);
          N2 = N1 - numel(gatedata);
          gate_fr(i,j) = N2/N1;
          zerodata = [zerodata;gatedata];
    elseif (gate==2)
        N1 = numel(Fscatterdata{i,j});
        Y1s = quantile(log10(Sscatterdata{i,j}),0.25);
        Y2s = quantile(log10(Sscatterdata{i,j}),0.5);
        Y1f = quantile(log10(Fscatterdata{i,j}),0.25);
        Y2f = quantile(log10(Fscatterdata{i,j}),0.5);
        gatedata = find(log10(Fscatterdata{i,j})>Y2f);find(log10(Sscatterdata{i,j})>Y2s);find(log10(Fscatterdata{i,j})<Y1f);find(log10(Sscatterdata{i,j})<Y1s);
        N2 = N1 - numel(gatedata);
        gate_fr(i,j) = N2/N1;
        zerodata = [zerodata;gatedata];
    else
          gate_fr(i,j) = 1;
    end
    Fscatterdata{i,j}(zerodata) = [];
    Sscatterdata{i,j}(zerodata) = [];
    data{i,j}(zerodata) = [];
    cellsizedata(zerodata) = [];
%     gatedata = [find(log10(Fscatterdata{i,j})<4.4);find(log10(Fscatterdata{i,j})>4.8);find(log10(Sscatterdata{i,j})<4.4);find(log10(Sscatterdata{i,j})>4.8)];
%     Fscatterdata{i,j}(gatedata) = [];
%     Sscatterdata{i,j}(gatedata) = [];
%     data{i,j}(gatedata) = [];
%     cellsizedata(gatedata) = [];
    logdata{i,j} = log10(data{i,j});          % unnormalized log-data
    normdata{i,j} = data{i,j}./cellsizedata;  % normalized data
    normlogdata{i,j} = log10(normdata{i,j});  % normalized log-data
    histcount{i,j} = histc(normlogdata{i,j},min(normlogdata{i,j}):0.02:max(normlogdata{i,j}));
    x_axis{i,j} = min(normlogdata{i,j}):0.02:max(normlogdata{i,j});
    end
    end
end

%FSC/SSC cutoff: 5.3