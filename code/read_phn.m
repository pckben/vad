%PHNREAD Read TIMIT phoneme transcript file, only returns location of the
%        specified phonemes in P
%
% @params
%   fn  path to .PHN file
%   P   cell array contains the list of interested phonemes
%
% @return
%   S   structure 
%       .start [nx1]
%       .end   [nx1]
%       .phone {nx1}
%
% @author
%   Pham Chau Khoa
%   Created: 12 March 2011

function S = read_phn( fn, P )

if nargin < 2,
    P = { 'iy';'ih';'eh';'ey';'ae';'aa';'aw';'ay';'ah';'ao';...
          'oy';'ow';'uh';'uw';'ux';'er';'ax';'ix';'axr';'ax-h';... % vowels
          'b';'d';'g';'p';'t';'k';'dx';'q';...      % stops
          'jh';'ch';...                             % affricatives
          'bcl';'dcl';'gcl';'pcl';'tck';'kcl';'tcl';...   % stop closures
          's';'sh';'z';'zh';'f';'th';'v';'dh';...   % fricatives
          'm';'n';'ng';'em';'en';'eng';'nx';...     % nasals
          'l';'r';'w';'y';'hh';'hv';'el';...        % semivowels and glides
          };          
end

if (~isa(P,'cell')), error('Input argument "P" must be a cell.'); end

% Read annotation file
fid = fopen(fn,'rt');
S = textscan(fid,'%d %d %s');
S{1} = S{1} +1; % MATLAB indices start from 1
S{2} = S{2} +1;
fclose(fid);

SS = [];
ii = 1;

for i=1:size(S{3},1)
    % if this segment's phoneme is in P, add the row
    %
    if sum(strcmp(S{3}(i),P)) > 0
        SS.start(ii,1) = S{1}(i); 
        SS.end  (ii,1) = S{2}(i);
        SS.phone(ii,1) = S{3}(i);
        ii = ii + 1;
    else
%         fprintf('%s ',S{3}{i});
    end
end

S = SS;

end

