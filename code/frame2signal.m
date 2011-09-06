function x = frame2signal(r,Nw,Nsh)

Lr = length(r);
Lx = Nsh*(Lr-1)+Nw;
x = zeros(Lx,1);

for i=1:Lr
    sidx = (i-1)*Nsh+1;
    eidx = (i-1)*Nsh+Nw;
    x(sidx:eidx) = x(sidx:eidx)+ones(Nw,1)*r(i);
end

x = x>=round(max(x)/2);

end