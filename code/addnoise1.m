function [y,n_gain] = addnoise1(x,n,SNR,ex,x_mask,Nw,Nsh)

Px = sum(ex(x_mask)) / sum(x_mask);

% making n same length as x
L = length(x);
while length(n)<L, n = [n;n]; end
n = n(1:L);

% calculate noise power
N = floor((L-Nw)/Nsh)+1;
en = zeros(N,1);
for i=1:N
    ni = n(Nsh*(i-1)+(1:Nw));
    en(i) = sum(ni.^2);
end

sorted_en = sort(en);
bottom20  = sorted_en(round(0.2*N));
top20     = sorted_en(round(0.8*N));
en_thresh = bottom20 + 0.1*(top20-bottom20);
n_mask   = en > en_thresh;

Pn = sum(en(n_mask)) / sum(n_mask);

Pn1 = Px/(10^(SNR/10));
n_gain = sqrt(Pn1/Pn);
y = x + n_gain*n;

end