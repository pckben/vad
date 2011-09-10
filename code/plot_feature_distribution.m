function plot_feature_distribution(Fx,Fn)

    N = 1000; % number of bins
    
    [cx,xx] = hist(Fx,N);
    [cn,nn] = hist(Fn,N);
    
    % normalize cx, cn such that their integral area =1
    norm_cx = cx/trapz(cx);
    norm_cn = cn/trapz(cn);
    
    assert(trapz(norm_cx)-1 <= eps);
    assert(trapz(norm_cn)-1 <= eps);
    
    plot(xx,norm_cx,'b',nn,norm_cn,'r');
    
    legend('signal+noise','noise');
    
end