function E = energy(x,Nw,Nsh)

T = floor((length(x)-Nw)/Nsh);
E = zeros(T,1);

for t = 1:T
    xt   = x((1:Nw) + (t-1)*Nsh);
    E(t) = sum(xt.^2);
%     E(t) = sum(abs(xt));
end

end

