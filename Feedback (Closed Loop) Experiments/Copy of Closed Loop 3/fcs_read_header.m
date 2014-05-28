function value = fcs_read_header(HEADER,key)
% READ_HEADER Reads a given entry in a .FCS header array
%   READ_HEADER(HEADER) extracts a given parameter stored in a header by
%   name. HEADER is a string and the value returned is a string as well.

value = char(HEADER(find(ismember(HEADER,key)==1),2));