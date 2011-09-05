% VAD HANGOVER for performance calculation.
% min_s : minimum speech duration 
% min_ns: minimum non-speech duration
% durations are in term of number of samples
function newlabel = vadhangover(label,min_ns)
    
i = find(label==0,1,'first');
while ~isempty(i) && i<length(label)
    j = find(label(i+1:end)==1,1,'first');
    if j<min_ns
        label(i:i+j-1)=1;
    end
    i = i+j-1+find(label(i+j:end)==0,1,'first');
end

newlabel = label;

end