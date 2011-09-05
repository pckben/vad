function [ S ] = read_wrd( fn )
%READ_WRD Summary of this function goes here
%   Detailed explanation goes here


% Read annotation file
fid = fopen(fn,'rt');
S = textscan(fid,'%d %d %s');
S{1} = S{1} +1; % MATLAB indices start from 1
S{2} = S{2} +1;
fclose(fid);

SS = [];
ii = 1;

for i=1:size(S{3},1)
    SS.start(ii,1) = S{1}(i); 
    SS.end  (ii,1) = S{2}(i);
    SS.phone(ii,1) = S{3}(i);
    ii = ii + 1;
end

S = SS;

end

