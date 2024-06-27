% Exectued on:
% MATLAB Version: 24.1.0.2603908 (R2024a) Update 4
% Operating System: Linux 6.9.5 

Ts = 0.1;
par = [9.8; 0.4; 1.2; 0.3];
par_plant = par;
par_plant(3) = par(3)*1.25;

syms omega theta tau;
x_sym = [theta; omega];
u_sym = tau;
f_sym = f(x_sym, u_sym, par, Ts);
h_sym = h(x_sym);
A_fun = matlabFunction(jacobian(f_sym, x_sym));
B_fun = matlabFunction(jacobian(f_sym, u_sym));
C_fun = matlabFunction(jacobian(h_sym, x_sym));

[A, B, C] = getSSmatrices(A_fun, B_fun, C_fun, [0; 0], 0);
[Ahat, Bhat, Chat] = getAugSSmatrices(A, B, C);

mpcmodel = ss(Ahat,Bhat,Chat,[],Ts);
nominal = struct('X', [0; 0], 'U', 0, 'Y', 0, 'DX', [0; 0]);

Phat0 = [
    (1/2)^2  0.0      0.0;
    0.0      (1/2)^2  0.0;
    0.0      0.0      (1)^2;
];
Qhat = [
    (0.1)^2  0.0      0.0;
    0.0      (0.5)^2  0.0;
    0.0      0.0      (0.1)^2;
];
Rhat = (0.5)^2;

mpcverbosity off;
mympc = mpc(mpcmodel,Ts,20,2);
setEstimator(mympc,"custom");
setoutdist(mympc, "model", []);
mympc.MV.Min = -1.5;
mympc.MV.Max = +1.5;
mympc.Model.Nominal = nominal;
mympc.Weights.ManipulatedVariablesRate = sqrt(2.5);
mympc.Weights.OutputVariables = sqrt(0.5);
myestim = mpcstate(mympc);

mympc.Optimizer.Algorithm = "active-set";

myestim.Plant = [];
myestim.Disturbance = [];
myestim.LastMove = [];
x0 = [0; 0];
xhat0 = [0; 0; 0];
[U_data, Y_data, R_data] = test_mpc(mympc, myestim, mpcmodel, nominal, ...
    xhat0, Phat0, Qhat, Rhat, @f, @h, par_plant, x0, Ts, ...
    A_fun, B_fun, C_fun, par, 0);

T_data= Ts*(0:length(Y_data)-1);
figure();
subplot(121)
plot(T_data, Y_data, T_data, R_data, '--')
subplot(122)
stairs(T_data,U_data)
hold on
plot(T_data, -1.5*ones(1, length(T_data)),T_data,+1.5*ones(1, length(T_data)))
hold off

myf = @() test_mpc(mympc, myestim, mpcmodel, nominal, ...
    xhat0, Phat0, Qhat, Rhat, @f, @h, par_plant, x0, Ts, ...
    A_fun, B_fun, C_fun, par, 0);
btime_slmpc_track_solver_AS = repeated_timeit(myf, 1) %#ok<NOPTS> 

myestim.Plant = [];
myestim.Disturbance = [];
myestim.LastMove = [];
x0 = [pi; 0];
xhat0 = [pi; 0; 0];
[U_data, Y_data, R_data] = test_mpc(mympc, myestim, mpcmodel, nominal, ...
    xhat0, Phat0, Qhat, Rhat, @f, @h, par_plant, x0, Ts, ...
    A_fun, B_fun, C_fun,  par, 10);

T_data= Ts*(0:length(Y_data)-1);
figure();
subplot(121)
plot(T_data, Y_data, T_data, R_data, '--')
subplot(122)
stairs(T_data,U_data)
hold on
plot(T_data, -1.5*ones(1, length(T_data)),T_data,+1.5*ones(1, length(T_data)))
hold off

myf = @() test_mpc(mympc, myestim, mpcmodel, nominal, ...
    xhat0, Phat0, Qhat, Rhat, @f, @h, par_plant, x0, Ts, ...
    A_fun, B_fun, C_fun, par, 10);
btime_slmpc_regul_solver_AS = repeated_timeit(myf, 1) %#ok<NOPTS> 

function [U_data, Y_data, R_data] = test_mpc(mympc, myestim, mpcmodel, nominal, ...
    xhat0, Phat, Qhat, Rhat, f, h, par, x0, Ts, ...
    A_fun, B_fun, C_fun, par_model, ystep)
    N = 35;
    r = 180;
    Y_data = zeros(1,N);
    U_data = zeros(1,N);
    R_data = zeros(1,N);
    u = 0;
    x = x0;
    xhat = xhat0;
    xhatop = zeros(size(xhat));
    dxhatop = zeros(size(xhat));
    xhat_d = xhat(1:2);
    [A, B, C] = getSSmatrices(A_fun, B_fun, C_fun, xhat_d, u);
    [Ahat, Bhat, Chat] = getAugSSmatrices(A, B, C);
    xhatop(1:2) = xhat_d;
    uop = u;
    yop = h(xhat_d);
    dxhatop(1:2) = f(xhat_d, u, par_model, Ts) - xhat_d;
    mpcmodel.A = Ahat;
    mpcmodel.B = Bhat;
    mpcmodel.C = Chat;
    nominal.X = xhatop;
    nominal.U = uop;
    nominal.Y = yop;
    nominal.DX = dxhatop;
    for i=1:N
        y = h(x) + ystep;
        myestim.Plant = xhat;
        u = mpcmoveAdaptive(mympc, myestim, mpcmodel, nominal, [], r);
        xhat_d = xhat(1:2);
        [A, B, C] = getSSmatrices(A_fun, B_fun, C_fun, xhat_d, u);
        [Ahat, Bhat, Chat] = getAugSSmatrices(A, B, C);
        xhatop(1:2) = xhat_d;
        uop = u;
        yop = h(xhat_d);
        dxhatop(1:2) = f(xhat_d, u, par_model, Ts) - xhat_d;
        mpcmodel.A = Ahat;
        mpcmodel.B = Bhat;
        mpcmodel.C = Chat;
        nominal.X = xhatop;
        nominal.U = uop;
        nominal.Y = yop;
        nominal.DX = dxhatop;
        U_data(:,i) = u;
        Y_data(:,i) = y;
        R_data(:,i) = r;
        x = f(x, u, par, Ts);
        [xhat, Phat] = updateKF(Ahat, Bhat, Chat, xhat, xhatop, uop, yop, dxhatop, ...
                                u, y, Phat, Qhat, Rhat);
    end
end

function [A, B, C] = getSSmatrices(A_fun, B_fun, C_fun, xop, uop)
    theta_op = xop(1);
    omega_op = xop(2);
    tau_op = uop;
    A = A_fun(omega_op, tau_op, theta_op);
    B = B_fun(omega_op, tau_op, theta_op);
    C = C_fun();
end

function [Ahat, Bhat, Chat] = getAugSSmatrices(A, B, C)
    Ahat = [A B; zeros(1, size(A, 2)), 1];
    Bhat = [B; zeros(size(B, 2), 1)];
    Chat = [C zeros(size(C, 1), 1)];
end

function [xhatNext, PhatNext] = updateKF(Ahat, Bhat, Chat, xhat, ...
                                         xhatop, uop, yop, dxhatop, ...
                                         u, y, Phat, Qhat, Rhat)
    M = (Phat*Chat')/(Chat*Phat*Chat'+Rhat);
    Khat  = Ahat*M;
    xhat0 = xhat-xhatop;
    yhat  = Chat*xhat0 + yop;
    xhatNext = Ahat*xhat0 + Bhat*(u-uop) + Khat*(y-yhat) + dxhatop + xhatop;
    PhatNext = Ahat*(Phat-M*Chat*Phat)*Ahat' + Qhat;
end