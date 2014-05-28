function meandata = RunNowClosed(step,samples)
fss = 0.018;
FileFolder = '.\FCS\';
files = dir(strcat(FileFolder,'*.fcs'));
n = length(files)/samples;
if step <= n
    tot_cols = n
    tot_rows = samples;
    [plate_ordered rows cols] = fcs_plate_read2(FileFolder,tot_cols);
    [data logdata normlogdata normdata] = process_plate(plate_ordered,tot_rows,tot_cols,1);
    for i = 1:rows
        j = step;
        meandata(i,1) = median(normdata{i,j})/fss;
    end
else
    meandata = [];
end
