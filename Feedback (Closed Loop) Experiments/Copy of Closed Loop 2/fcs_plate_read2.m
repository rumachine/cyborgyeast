function [plate_ordered rows cols] = fcs_plate_read2(folder,cols)
% FCS_PLATE_READ2   Reads a directory containing .FCS files
% FCS_PLATE_READ2(folder) loads all wells in a plate in a directory into 
% a struct array containing each well's HEADER, PARAMS and DATA.
% 'cols' specifies how many wells per row have been measured. Rows do not
% need to be consecutive, but the wells in each row do.

files = dir(strcat(folder,'*.fcs'));
n = length(files)
letters = cell({'A','B','C','D','E','F','G','H'});
if cols > 12
badletters = cell({'B','D','F','H'});
end


for i=1:n
    file_name = strcat(folder, char(files(i).name));                   
    [DATA,PARAMS,HEADER]=fcsread(file_name);
    plate(i) = struct('DATA',{DATA.*(DATA>0)},'PARAMS',{PARAMS},'HEADER',{HEADER})    ;
end

for i = 1:n
    w_p_str = fcs_read_header(plate(i).HEADER,'TUBE NAME');
    row{i} = w_p_str(1);
    if cols > 12
        if i > 1
            w_p_str_old = fcs_read_header(plate(i-1).HEADER,'TUBE NAME');
            rowold{i} = w_p_str_old(1);
        end
        if ismember(row{i},badletters)
            row{i} = rowold{i};
            newnum = num2str(12+str2num(w_p_str(2:end)));
            plate(i).HEADER(find(ismember(plate(i).HEADER,'TUBE NAME')==1),2);
            plate(i).HEADER(find(ismember(plate(i).HEADER,'TUBE NAME')==1),2) = cell({strcat(row{i},newnum)});
        end
    end
end
rows = unique(row);
k = length(rows);

for i = 1:k
    for j=1:cols
        w_p_str = fcs_read_header(plate((i-1)*cols+j).HEADER,'TUBE NAME')     ;% String containing the position of the well.
        row = i                ;% Computes the row of the well being read
        col = str2num(w_p_str(2:end))                              ;% Computes the column of the well being read
        plate_ordered{i,col}= plate((i-1)*cols+j);
    end
end
[rows cols] = size(plate_ordered);

