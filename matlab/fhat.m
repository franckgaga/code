function xhatnext = fhat(xhat, u, par, Ts)
    xd = xhat(1:2);
    xs = xhat(3);
    uhat = u + xs;
    xdnext = f(xd, uhat, par, Ts);
    xsnext = xs;
    xhatnext = [xdnext; xsnext];
end