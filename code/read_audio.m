
function x = read_audio(filepath,type)
    
    if nargin<2, type='sphere'; end
    fid = fopen(filepath,'r');
    
    if strcmpi(type,'aurora')==1
        % AURORA:    
        x = fread(fid,inf,'int16',0,'ieee-be');
    elseif strcmpi(type,'sphere')==1
        % SPHERE:
        x = fread(fid,inf,'int16');
        x = x(513:end); % remove SPHERE header
    elseif strcmpi(type,'wav')==1
        % WAV:
        x = wavread(filepath);
    elseif strcmpi(type,'pcm')==1
        % PCM:
        x = fread(fid,inf,'int16');
    end
    
    fclose(fid);
end