function xnext = f(x, u, par, Ts)
    k1 = pendulum(par, x, u);
    k2 = pendulum(par, x + 0.5*Ts*k1, u);
    k3 = pendulum(par, x + 0.5*Ts*k2, u);
    k4 = pendulum(par, x + Ts*k3, u);
    xnext = x + Ts/6*(k1 + 2*k2 + 2*k3 + k4);
end

function x = pendulum(par, x, u)
    g = par(1); L=par(2); K=par(3); m=par(4);
    theta = x(1); omega = x(2);
    tau = u(1);
    x = [omega; -g/L*sin(theta) - K/m*omega + tau/m/L^2];
end
