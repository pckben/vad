% calculate the Bhattacharyya coefficient to determine the closeness of the
% two distribution.
% The implementation follows Wikipedia article on the same topic.
%
function d = bhattacharyya(x,n)

N = 100; % number of partitions
    
m = min([x;n]);
M = max([x;n]);
step = (M-m)/N;
BC = 0; % bhattacharyya coefficient

for i=m:step:M
    xi = sum(x>=i & x<i+step);
    ni = sum(n>=i & n<i+step);
    
    BC = BC+sqrt(xi*ni);
end

d = -log(BC);

end