# ==========================================
# ========== STATE ESTIMATOR ===============
# ==========================================
using ModelPredictiveControl
function pendulum!(ẋ, x, u, par)
    g, L, K, m = par     # [m/s²], [m], [kg/s], [kg]
    θ, ω = x[1], x[2]    # [rad], [rad/s]
    τ = u[1]             # [Nm]
    ẋ[1] = ω
    ẋ[2] = -g/L*sin(θ) - K/m*ω + τ/m/L^2
    return nothing
end
const par = (9.8, 0.4, 1.2, 0.3)
f!(ẋ, x, u, _ ) = pendulum!(ẋ, x, u, par)
h!(y, x, _ ) = (y[1] = 180/π*x[1])   # [°]
nu = 1; nx = 2; ny = 1; Ts = 0.1
model = NonLinModel(f!, h!, Ts, nu, nx, ny)
vu = ["\$τ\$ (Nm)"]
vx = ["\$θ\$ (rad)", "\$ω\$ (rad/s)"]
vy = ["\$θ\$ (°)"]
model = setname!(model; u=vu, x=vx, y=vy)

## =========================================
σQ = [0.1, 0.5]; σR=[0.5]; nint_u=[1]; σQint_u=[0.1]
estim = UnscentedKalmanFilter(model; σQ, σR, nint_u, σQint_u)

## =========================================
const par_plant = (par[1], par[2], 1.25*par[3], par[4])
f_plant!(ẋ, x, u, _) = pendulum!(ẋ, x, u, par_plant)
plant = NonLinModel(f_plant!, h!, Ts, nu, nx, ny)
N = 35; u=[0.5]; res = sim!(estim, N, u; plant, y_noise=[0.5])
using Plots; plot(res, plotu=false, plotxwithx̂=true)

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures
theme(:default)
default(fontfamily="Computer Modern");
p = plot(res, plotu=false, plotxwithx̂=true, size=(425, 275))
yticks!(p[2], [0.0, 0.25, 0.5, 0.75])
yticks!(p[3], [-0.5, 0, 0.5, 1.0, 1.5])
yticks!(p[4], [-0.10, -0.05, 0, 0.05, 0.1])
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_NonLinMPC1.pdf")

## ==========================================
## ========== NONLINEAR MPC =================
## ==========================================

## =========================================
Hp, Hc, Mwt, Nwt, Cwt = 20, 2, [0.5], [2.5], Inf
nmpc = NonLinMPC(estim; Hp, Hc, Mwt, Nwt, Cwt)
umin, umax = [-1.5], [+1.5]
nmpc = setconstraint!(nmpc; umin, umax)

## =========================================
using JuMP; unset_time_limit_sec(nmpc.optim)

## =========================================
x_0 = [0, 0]; x̂_0 = [0, 0, 0]; ry = [180]
res_ry = sim!(nmpc, N, ry; plant, x_0, x̂_0)
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
using JuMP, Ipopt, KNITRO

optim = JuMP.Model(Ipopt.Optimizer, add_bridges=false)
nmpc_ipopt = NonLinMPC(estim; Hp, Hc, Mwt, Nwt, Cwt, optim)
nmpc_ipopt = setconstraint!(nmpc_ipopt; umin, umax)
JuMP.unset_time_limit_sec(nmpc_ipopt.optim)

optim = JuMP.Model(KNITRO.Optimizer, add_bridges=false)
set_attribute(optim, "algorithm", "4") # 4th algorithm is SQP
nmpc_knitro = NonLinMPC(estim; Hp, Hc, Mwt, Nwt, Cwt, optim)
nmpc_knitro = setconstraint!(nmpc_knitro; umin, umax)
JuMP.unset_time_limit_sec(nmpc_knitro.optim)

bm = @benchmark(
        sim!($nmpc_ipopt, $N, $ry; plant=$plant, x_0=$x_0, x̂_0=$x̂_0),
        samples=50, 
        seconds=10*60
    )
@show btime_NMPC_track_solver_IP = median(bm)

bm = @benchmark(
        sim!($nmpc_knitro, $N, $ry; plant=$plant, x_0=$x_0, x̂_0=$x̂_0),
        samples=50,
        seconds=10*60
    )
@show btime_NMPC_track_solver_SQ = median(bm)


## =========================================
x_0 = [π, 0]; x̂_0 = [π, 0, 0]; y_step = [10]
res_yd = sim!(nmpc, N, [180.0]; plant, x_0, x̂_0, y_step)
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
bm = @benchmark(
        sim!($nmpc_ipopt, $N, $[180.0]; plant=$plant, x_0=$x_0, x̂_0=$x̂_0, y_step=$y_step),
        samples=50,
        seconds=10*60
    )
@show btime_NMPC_regul_solver_IP = median(bm)

bm = @benchmark(
        sim!($nmpc_knitro, $N, $[180.0]; plant=$plant, x_0=$x_0, x̂_0=$x̂_0, y_step=$y_step),
        samples=50,
        seconds=10*60
    )
@show btime_NMPC_regul_solver_SQ = median(bm)

# ==========================================
# ========== ECONOMIC MPC ==================
# ==========================================
h2!(y, x, _ ) = (y[1] = 180/π*x[1]; y[2]=x[2])
nu, nx, ny = 1, 2, 2
model2 = NonLinModel(f!      , h2!, Ts, nu, nx, ny)
plant2 = NonLinModel(f_plant!, h2!, Ts, nu, nx, ny)
model2 = setname!(model2, u=vu, x=vx, y=[vy; vx[2]])
plant2 = setname!(plant2, u=vu, x=vx, y=[vy; vx[2]])
estim2 = UnscentedKalmanFilter(model2; σQ, σR, nint_u, σQint_u, i_ym=[1])


## =========================================
function JE(UE, ŶE, _ )
    τ, ω = UE[1:end-1], ŶE[2:2:end-1]
    return Ts*sum(τ.*ω)
end
Mwt2, Ewt = [Mwt; 0.0], 3.5e3
empc = NonLinMPC(estim2; Hp, Hc, Nwt, Mwt=Mwt2, Cwt, JE, Ewt)
empc = setconstraint!(empc; umin, umax)

## =========================================
using JuMP; unset_time_limit_sec(empc.optim)

## =========================================
x_0 = [0, 0]; x̂_0 = [0, 0, 0]; ry = [180; 0]
res2_ry = sim!(empc, N, ry; plant=plant2, x_0, x̂_0)
plot(res2_ry, ploty=[1])

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures 
theme(:default)
default(fontfamily="Computer Modern")
p = plot(res2_ry, ploty=[1], size=(425, 200), bottom_margin=10px)
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_EconomMPC1.pdf")

## =========================================
function calcW(res)
    τ, ω = res.U_data[1, 1:end-1], res.X_data[2, 1:end-1]
    return Ts*sum(τ.*ω)
end
display(Dict(:W_nmpc => calcW(res_ry), :W_empc => calcW(res2_ry)))

## =========================================
## ========= Benchmark =====================
## =========================================
using BenchmarkTools
using JuMP, Ipopt, KNITRO

optim = JuMP.Model(Ipopt.Optimizer, add_bridges=false)
empc_ipopt = NonLinMPC(estim2; Hp, Hc, Nwt, Mwt=Mwt2, Cwt, JE, Ewt, optim)
empc_ipopt = setconstraint!(empc_ipopt; umin, umax)
JuMP.unset_time_limit_sec(empc_ipopt.optim)

optim = JuMP.Model(KNITRO.Optimizer, add_bridges=false)
set_attribute(optim, "algorithm", "4") # 4th algorithm is SQP
empc_knitro = NonLinMPC(estim2; Hp, Hc, Nwt, Mwt=Mwt2, Cwt, JE, Ewt, optim)
empc_knitro = setconstraint!(empc_knitro; umin, umax)
JuMP.unset_time_limit_sec(empc_knitro.optim)

bm = @benchmark(
        sim!($empc_ipopt, $N, $ry; plant=$plant2, x_0=$x_0, x̂_0=$x̂_0),
        samples=50, 
        seconds=10*60
    )
@show btime_EMPC_track_solver_IP = median(bm)

bm = @benchmark(
        sim!($empc_knitro, $N, $ry; plant=$plant2, x_0=$x_0, x̂_0=$x̂_0),
        samples=50,
        seconds=10*60
    )
@show btime_EMPC_track_solver_SQ = median(bm)

## =========================================
x_0 = [π, 0]; x̂_0 = [π, 0, 0]; y_step = [10; 0]
res2_yd = sim!(empc, N, ry; plant=plant2, x_0, x̂_0, y_step)
plot(res2_yd, ploty=[1])

## =========================================
display(Dict(:W_nmpc => calcW(res_yd), :W_empc => calcW(res2_yd)))

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures
theme(:default)
default(fontfamily="Computer Modern")
p = plot(res2_yd, ploty=[1], size=(425, 200), bottom_margin=10px)
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_EconomMPC2.pdf")

## =========================================
## ========= Benchmark =====================
## =========================================
bm = @benchmark(
        sim!($empc_ipopt, $N, $ry; plant=$plant2, x_0=$x_0, x̂_0=$x̂_0, y_step=$y_step),
        samples=50,
        seconds=10*60
    )
@show btime_EMPC_regul_solver_IP = median(bm)

bm = @benchmark(
        sim!($empc_knitro, $N, $ry; plant=$plant2, x_0=$x_0, x̂_0=$x̂_0, y_step=$y_step),
        samples=50,
        seconds=10*60
    )
@show btime_EMPC_regul_solver_SQ = median(bm)

## ==========================================
## ====== SUCCESSIVE LINEARIZATION MPC ======
## ==========================================
using Pkg; Pkg.add(["JuMP","DAQP"]) # install JuMP and DAQP
using JuMP, DAQP
optim = JuMP.Model(DAQP.Optimizer, add_bridges=false)

## ==========================================
linmodel = linearize(model, x=[0, 0], u=[0])
kf = KalmanFilter(linmodel; σQ, σR, nint_u, σQint_u)
mpc3 = LinMPC(kf; Hp, Hc, Mwt, Nwt, Cwt, optim)
mpc3 = setconstraint!(mpc3; umin, umax)

## ==========================================
function sim2!(mpc, nlmodel, N, ry, plant, x_0, x̂_0, y_step)
    U, Y, Ry = zeros(1, N), zeros(1, N), zeros(1, N)
    u, x̂ = [0], x̂_0
    initstate!(mpc, u, plant())
    setstate!(plant, x_0); setstate!(mpc, x̂_0)
    linmodel = linearize(nlmodel; u, x=x̂[1:2])
    setmodel!(mpc, linmodel)
    for i = 1:N
        y = plant() + y_step
        u = mpc(ry)
        linearize!(linmodel, nlmodel; u, x=x̂[1:2])
        setmodel!(mpc, linmodel) 
        U[:,i], Y[:,i], Ry[:,i] = u, y, ry
        x̂ = updatestate!(mpc, u, y)
        updatestate!(plant, u)
    end
    U_data, Y_data, Ry_data = U, Y, Ry
    return SimResult(mpc, U_data, Y_data; Ry_data)
end

## ==========================================
x_0 = [0, 0]; x̂_0 = [0, 0, 0]; ry = [180]; y_step=[0]
res3_ry = sim2!(mpc3, model, N, ry, plant, x_0, x̂_0, y_step)
plot(res3_ry)

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures 
theme(:default)
default(fontfamily="Computer Modern")
p = plot(res3_ry, size=(425, 200), bottom_margin=10px)
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_SuccLinMPC1.pdf")

## =========================================
## ========= Benchmark =====================
## =========================================
using BenchmarkTools

x_0 = [0, 0]; x̂_0 = [0, 0, 0]; ry = [180]; y_step=[0]
bm = @benchmark(
        sim2!($mpc3, $model, $N, $ry, $plant, $x_0, $x̂_0, $y_step),
        samples=500, 
        seconds=10*60
    )
@show btime_SLMPC_track_solver_AS = median(bm)

## =========================================
x_0 = [π, 0]; x̂_0 = [π, 0, 0]; ry = [180]; y_step=[10]
res3_yd = sim2!(mpc3, model, N, ry, plant, x_0, x̂_0, y_step)
plot(res3_yd)

## =========================================
## ========= Plot PDF ======================
## =========================================
using PlotThemes, Plots.PlotMeasures 
theme(:default)
default(fontfamily="Computer Modern")
p = plot(res3_yd, size=(425, 200), bottom_margin=10px)
display(p)
savefig(p, "$(@__DIR__())/../../fig/plot_SuccLinMPC2.pdf")

## =========================================
## ========= Benchmark =====================
## =========================================
x_0 = [π, 0]; x̂_0 = [π, 0, 0]; ry = [180]; y_step=[10]
bm = @benchmark(
        sim2!($mpc3, $model, $N, $ry, $plant, $x_0, $x̂_0, $y_step),
        samples=500, 
        seconds=10*60
    )
@show btime_SLMPC_track_solver_AS = median(bm)