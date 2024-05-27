# ==========================================
# ========== STATE ESTIMATOR ===============
# ==========================================
using ModelPredictiveControl
function pendulum(par, x, u)
    g, L, K, m = par     # [m/s²], [m], [kg/s], [kg]
    θ, ω = x[1], x[2]    # [rad], [rad/s]
    τ = u[1]             # [Nm]
    return [ω, -g/L*sin(θ) - K/m*ω + τ/m/L^2]
end
par = (9.8, 0.4, 1.2, 0.3); Ts = 0.1; nu = 1; nx = 2; ny = 1
f(x, u, _ ) = pendulum(par, x, u)
h(x, _ )    = [180/π*x[1]]  # [°]
model = NonLinModel(f, h, Ts, nu, nx, ny)
vu = ["\$τ\$ (Nm)"]
vx = ["\$θ\$ (rad)", "\$ω\$ (rad/s)"]
vy = ["\$θ\$ (°)"]
model = setname!(model; u=vu, x=vx, y=vy)

## =========================================
α=0.01; σQ = [0.1, 0.5]; σR=[0.5]; nint_u=[1]; σQint_u=[0.1]
estim = UnscentedKalmanFilter(model; α, σQ, σR, nint_u, σQint_u)

## =========================================
const par_plant = (par[1], par[2], 1.25*par[3], par[4])
f_plant(x, u, _) = pendulum(par_plant, x, u)
plant = NonLinModel(f_plant, h, Ts, nu, nx, ny)
N = 35; u=[0.5]; res = sim!(estim, N, u; plant, y_noise=[0.5])
using Plots; plot(res, plotu=false, plotxwithx̂=true)

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures
theme(:default)
default(fontfamily="Computer Modern");
p = plot(res, plotu=false, plotxwithx̂=true, size=(475, 300))
yticks!(p[2], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
yticks!(p[3], [-0.3, 0.0, 0.3, 0.6, 0.9, 1.2])
yticks!(p[4], [-0.075, -0.05, -0.025, 0.0, 0.025, 0.05, 0.075])
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_NonLinMPC1.pdf")

## ==========================================
## ========== CONTROLLER ====================
## ==========================================

## =========================================
Hp, Hc, Mwt, Nwt = 20, 2, [0.5], [2.5]
mpc = NonLinMPC(estim; Hp, Hc, Mwt, Nwt, Cwt=Inf)
mpc = setconstraint!(mpc, umin=[-1.5], umax=[+1.5]) 

## =========================================
x_0 = [0, 0]; x̂_0 = [0, 0, 0]
using JuMP; unset_time_limit_sec(mpc.optim) # hide
res_ry = sim!(mpc, N, [180]; plant, x_0, x̂_0)
plot(res_ry)

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures 
theme(:default)
default(fontfamily="Computer Modern")
p = plot(res_ry, size=(425, 200), bottom_margin=10px)
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_NonLinMPC2.pdf")

## =========================================
## ========= Benchmark =====================
## =========================================
using BenchmarkTools
x_0 = [0, 0]; x̂_0 = [0, 0, 0]
bm = @benchmark sim!($mpc, $N, [180]; plant=$plant, x_0=$x_0, x̂_0=$x̂_0) samples= seconds=5*60
display(bm)

## =========================================
x_0 = [π, 0]; x̂_0 = [π, 0, 0]
res_yd = sim!(mpc, N, [180.0]; plant, x_0, x̂_0, y_step=[10])
plot(res_yd)

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures
theme(:default)
default(fontfamily="Computer Modern")
p = plot(res_yd, size=(425, 200), bottom_margin=10px)
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_NonLinMPC3.pdf")

## =========================================
## ========= Benchmark =====================
## =========================================
using BenchmarkTools
x_0 = [π, 0]; x̂_0 = [π, 0, 0]
bm = @benchmark sim!($mpc, $N, [180.0]; plant=$plant, x_0=$x_0, x̂_0=$x̂_0, y_step=[10]) seconds=30
display(bm)

#=
## =========================================
## ========= EMPC ==========================
## =========================================

h2(x, _ ) = [180/π*x[1], x[2]]
nu, nx, ny = 1, 2, 2
model2 = NonLinModel(f      , h2, Ts, nu, nx, ny)
plant2 = NonLinModel(f_plant, h2, Ts, nu, nx, ny)
estim2 = UnscentedKalmanFilter(model2; σQ, σR, nint_u, σQint_u, i_ym=[1])

function JE(UE, ŶE, _ )
    τ, ω = UE[1:end-1], ŶE[2:2:end-1]
    return Ts*sum(τ.*ω)
end
empc = NonLinMPC(estim2, Hp=20, Hc=2, Mwt=[0.5, 0], Nwt=[2.5], Cwt=Inf, Ewt=3.5e3, JE=JE)
empc = setconstraint!(empc, umin=[-1.5], umax=[+1.5])

using JuMP; unset_time_limit_sec(empc.optim) # hide
res2_ry = sim!(empc, N, [180, 0], plant=plant2, x0=[0, 0], x̂0=[0, 0, 0])
plot(res2_ry)

function calcW(res)
    τ, ω = res.U_data[1, 1:end-1], res.X_data[2, 1:end-1]
    return Ts*sum(τ.*ω)
end
Dict(:W_nmpc => calcW(res_ry), :W_empc => calcW(res2_ry))

res2_yd = sim!(empc, N, [180; 0]; plant=plant2, x0=[π, 0], x̂0=[π, 0, 0], y_step=[10, 0])
plot(res2_yd)

Dict(:W_nmpc => calcW(res_yd), :W_empc => calcW(res2_yd))
=#
