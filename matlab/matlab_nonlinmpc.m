Ts = 0.1;
par = [9.8; 0.4; 1.2; 0.3];
par_plant = par;
par_plant(3) = par(3)*1.25;

x0 = [0; 0];
u = 0.5;
Y = test_pendulum(par, Ts, u, x0);
T = Ts*(0:length(Y)-1);
figure()
plot(T,Y)

xhat0 = [0; 0; 0];
ukf = unscentedKalmanFilter(@fhat, @hhat, xhat0);
ukf.Alpha = 0.01;
ukf.StateCovariance = [
    (1/2)^2  0.0      0.0;
    0.0      (1/2)^2  0.0;
    0.0      0.0      (1)^2;
];
ukf.ProcessNoise = [
    (0.1)^2  0.0      0.0;
    0.0      (0.5)^2  0.0;
    0.0      0.0      (0.1)^2;
];
ukf.MeasurementNoise = (0.5)^2;

x0 = [0; 0];
xhat0 = [0; 0; 0];
ynoise = 0.5;
[Y, Yhat, ~, X, Xhat] = test_ukf(ukf, ...
    par_plant, par, Ts, u, x0, xhat0, ynoise);

figure();
subplot(3,2, [1 3 5])
plot(T, Y, T, Yhat)
subplot(3,2,2)
plot(T, X(1,:), T, Xhat(1,:))
subplot(3,2,4)
plot(T, X(2,:), T, Xhat(2,:))
subplot(3,2,6)
plot(T, Xhat(3,:))

myf = @() test_ukf(ukf, ...
    par_plant, par, Ts, u, x0, xhat0, ynoise);

%btime_ukf = repeated_timeit(myf, 100);
%display(btime_ukf)




mympc = nlmpc(3,1,1);

mympc.PredictionHorizon = 20;
mympc.ControlHorizon = 2;
mympc.Weights.ManipulatedVariablesRate = sqrt(2.5);
mympc.Weights.OutputVariables = sqrt(0.5);
mympc.Model.StateFcn = @fhat;
mympc.Model.IsContinuousTime = false;
mympc.Model.OutputFcn = @hhat;
mympc.Model.NumberOfParameters = 2;
mympc.ManipulatedVariables.Min = -1.5;
mympc.ManipulatedVariables.Max = +1.5;
%opt = optimoptions( ...
%    'fmincon', 'Algorithm', ...
%    'interior-point', ...
%    'Display', 'none', ...
%    'SpecifyObjectiveGradient', true, ...
%    'SpecifyConstraintGradient', true ...
%);
%mympc.Optimization.SolverOptions = opt;

validateFcns(mympc, [0; 0; 0], 0.5, [], {par, Ts});

opt = nlmpcmoveopt;
opt.Parameters = {par, Ts};


[R, Y, ~, U, ~, ~] = test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, 0);
figure();
subplot(121)
plot(T, Y, T, R, '--')
subplot(122)
stairs(T,U)
hold on
plot(T, -1.5*ones(1, length(T)),T,+1.5*ones(1, length(T)))
hold off

myf = @() test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, 0);

btime_mpc_track = repeated_timeit(myf, 10);
display(btime_mpc_track)


x0 = [pi; 0];
xhat0 = [pi; 0; 0];
[R, Y, Yhat, U, X, Xhat] = test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, 10);
figure();
subplot(121)
plot(T, Y, T, R, '--')
subplot(122)
stairs(T,U)
hold on
plot(T, -1.5*ones(1, length(T)),T,+1.5*ones(1, length(T)))
hold off

myf = @() test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, 10);

btime_mpc_track = repeated_timeit(myf, 10);
display(btime_mpc_track)

function x = pendulum(par, x, u)
    g = par(1); L=par(2); K=par(3); m=par(4);
    theta = x(1); omega = x(2);
    tau = u(1);
    x = [omega; -g/L*sin(theta) - K/m*omega + tau/m/L^2];
end

function xnext = f(x, u, par, Ts)
    k1 = pendulum(par, x, u);
    k2 = pendulum(par, x + 0.5*Ts*k1, u);
    k3 = pendulum(par, x + 0.5*Ts*k2, u);
    k4 = pendulum(par, x + Ts*k3, u);
    xnext = x + Ts/6*(k1 + 2*k2 + 2*k3 + k4);
end

function y = h(x)
    y = 180/pi*x(1);
end

function xhatnext = fhat(xhat, u, par, Ts)
    xd = xhat(1:2);
    xs = xhat(3);
    uhat = u + xs;
    xdnext = f(xd, uhat, par, Ts);
    xsnext = xs;
    xhatnext = [xdnext; xsnext];
end

function yhat = hhat(xhat, ~, ~, ~)
    xhatd = xhat(1:2);
    yd = h(xhatd);
    yhat = yd;
end

function Y = test_pendulum(par, Ts, u, x0)
    N = 35;
    Y = zeros(1, N);
    x = x0;
    for i=1:N
        y = h(x);
        Y(:, i) = y;
        x = f(x, u, par, Ts);
    end
end

function [Y, Yhat, U, X, Xhat] = test_ukf(ukf, ...
    par_plant, par, Ts, u, x0, xhat0, ynoise)
    N = 35;
    Y = zeros(1, N);
    Yhat = zeros(1, N);
    U = zeros(1, N);
    X = zeros(2, N);
    Xhat = zeros(3, N);
    x = x0;
    xhat = xhat0;
    for i=1:N
        y = h(x) + ynoise*randn(1, 1);
        yhat = hhat(xhat);
        Y(:, i) = y;
        Yhat(:, i) = yhat;
        U(:, i) = u;
        X(:, i) = x;
        Xhat(:, i) = xhat;
        correct(ukf, y, u, par, Ts);
        xhat = predict(ukf, u, par, Ts);
        x = f(x, u, par_plant, Ts);
    end
end

function [R, Y, Yhat, U, X, Xhat] = test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, ystep)
    N = 35;
    R = zeros(1,N);
    Y = zeros(1, N);
    Yhat = zeros(1, N);
    U = zeros(1, N);
    X = zeros(2, N);
    Xhat = zeros(3, N);
    x = x0;
    xhat = xhat0;
    r = 180;
    for i=1:N
        y = h(x) + ystep;
        yhat = hhat(xhat);
        u = nlmpcmove(mympc, xhat, u, r, [], opt);
        R(:, i) = r;
        Y(:, i) = y;
        Yhat(:, i) = yhat;
        U(:, i) = u;
        X(:, i) = x;
        Xhat(:, i) = xhat;
        correct(ukf, y, u, par, Ts);
        xhat = predict(ukf, u, par, Ts);
        x = f(x, u, par_plant, Ts);
    end
end

