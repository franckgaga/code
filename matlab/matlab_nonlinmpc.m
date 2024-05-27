% Exectued on:
% MATLAB Version: 24.1.0.2603908 (R2024a) Update 3
% MATLAB License Number: 968398
% Operating System: Linux 6.8.0-76060800daily20240311-generic 

Ts = 0.1;
par = [9.8; 0.4; 1.2; 0.3];
par_plant = par;
par_plant(3) = par(3)*1.25;

xhat0 = [0; 0; 0];
ukf = unscentedKalmanFilter(@fhat, @hhat, xhat0);
ukf.Alpha = 0.01;
Phat0 = [
    (1/2)^2  0.0      0.0;
    0.0      (1/2)^2  0.0;
    0.0      0.0      (1)^2;
];
ukf.StateCovariance = Phat0;
ukf.ProcessNoise = [
    (0.1)^2  0.0      0.0;
    0.0      (0.5)^2  0.0;
    0.0      0.0      (0.1)^2;
];
ukf.MeasurementNoise = (0.5)^2;

x0 = [0; 0];
xhat0 = [0; 0; 0];
u = 0.5;
ynoise = 0.5;
[Y, Yhat, ~, X, Xhat] = test_ukf(ukf, ...
    par_plant, par, Ts, u, x0, xhat0, Phat0, ynoise);
T = Ts*(0:length(Y)-1);

figure();
subplot(3,2, [1 3 5])
plot(T, Y, T, Yhat)
subplot(3,2,2)
plot(T, X(1,:), T, Xhat(1,:))
subplot(3,2,4)
plot(T, X(2,:), T, Xhat(2,:))
subplot(3,2,6)
plot(T, Xhat(3,:))

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

x0 = [0; 0];
xhat0 = [0; 0; 0];
u = 0.0;
[R, Y, ~, U, ~, ~] = test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, Phat0, 0);
T = Ts*(0:length(Y)-1);
figure();
subplot(121)
plot(T, Y, T, R, '--')
subplot(122)
stairs(T,U)
hold on
plot(T, -1.5*ones(1, length(T)),T,+1.5*ones(1, length(T)))
hold off

mympc.Optimization.SolverOptions.Algorithm = "interior-point";
myf = @() test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, Phat0, 0);
btime_nmpc_track_solver_IP = repeated_timeit(myf, 10) %#ok<NOPTS> 

mympc.Optimization.SolverOptions.Algorithm = "sqp";
myf = @() test_mpc(mympc, opt, ukf, ...
   par_plant, par, Ts, u, x0, xhat0, Phat0, 0);
btime_nmpc_track_solver_SQ = repeated_timeit(myf, 10) %#ok<NOPTS> 

x0 = [pi; 0];
xhat0 = [pi; 0; 0];
u = 0.0;
[R, Y, Yhat, U, X, Xhat] = test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, Phat0, 10);
figure();
subplot(121)
plot(T, Y, T, R, '--')
subplot(122)
stairs(T,U)
hold on
plot(T, -1.5*ones(1, length(T)),T,+1.5*ones(1, length(T)))
hold off

mympc.Optimization.SolverOptions.Algorithm = "interior-point";
myf = @() test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, Phat0, 0);
btime_nmpc_regul_solver_IP = repeated_timeit(myf, 10) %#ok<NOPTS> 

mympc.Optimization.SolverOptions.Algorithm = "sqp";
myf = @() test_mpc(mympc, opt, ukf, ...
    par_plant, par, Ts, u, x0, xhat0, Phat0, 0);
btime_nmpc_regul_solver_SQ = repeated_timeit(myf, 10) %#ok<NOPTS> 

function [Y, Yhat, U, X, Xhat] = test_ukf(ukf, ...
    par_plant, par, Ts, u, x0, xhat0, Phat0, ynoise)
    N = 35;
    Y = zeros(1, N);
    Yhat = zeros(1, N);
    U = zeros(1, N);
    X = zeros(2, N);
    Xhat = zeros(3, N);
    x = x0;
    xhat = xhat0;
    ukf.StateCovariance = Phat0;
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
    par_plant, par, Ts, u, x0, xhat0, Phat0, ystep)
    N = 35;
    R = zeros(1,N);
    Y = zeros(1, N);
    Yhat = zeros(1, N);
    U = zeros(1, N);
    X = zeros(2, N);
    Xhat = zeros(3, N);
    x = x0;
    xhat = xhat0;
    ukf.StateCovariance = Phat0;
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

