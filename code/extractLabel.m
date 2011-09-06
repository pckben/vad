function L = extractLabel(len,labels)

L = false(len,1);

for i=1:length(labels.start)
    L(labels.start(i):labels.end(i)) = true;
end

end

