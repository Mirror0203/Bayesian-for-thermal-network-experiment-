# inference_thermonly.jl
# Thermal-only posterior with learned node-wise σ_T (8).
# Parameter layout (length 43):
#   P[1:27]  -> thermal params (19 R + 8 C), in the SAME order as your models.jl expects
#   P[28:35] -> initial temperatures C_temp0 (8)
#   P[36:43] -> σ_T_vec (8), positive
#
# Other args:
#   f_temp_data       :: 8×N temperature window (°C)
#   time_window_size  :: Real, length of window in model "steps"
#   T_init_fixed      :: length-8 initial temps for the window (from first col of the window)
#   P_air_fixed       :: length-13 fixed airflow/occupancy vector (NOT inferred here)
#
# Returns: log posterior (Float64)

using LinearAlgebra
using Distributions

function evaluateLogPdf_thermonly(P::AbstractVector,
                                  f_temp_data,
                                  time_window_size::Real,
                                  T_init_fixed::AbstractVector,
                                  P_air_fixed::AbstractVector)

    # --- Basic checks & slicing ---
    length(P) == 43 || return -Inf
    length(T_init_fixed) == 8 || return -Inf
    length(P_air_fixed) == 13 || return -Inf

    P_therm   = @view P[1:27]      # 19 R + 8 C  (order must match solutionForward_thermal)
    C_temp0   = @view P[28:35]     # 8 initial temperatures (°C) – will be steered by prior to data[ :,1 ]
    σ_t_vec   = @view P[36:43]     # 8 noise scales (°C), must be > 0

    if any(!isfinite, P_therm) || any(!isfinite, C_temp0) ||
       any(x -> !isfinite(x) || x <= 0.0, σ_t_vec)
        return -Inf
    end

    # --- Forward model over the window ---
    tspan = (0.0, float(time_window_size))
    # solutionForward_thermal expects [R,C] then airflow block (fixed here)
    P_therm_full = vcat(P_therm, P_air_fixed)  # length 27 + 13 = 40
    sol   = solutionForward_thermal(P_therm_full, C_temp0, tspan)
    Umod  = reduce(hcat, sol.u)                          # 8×N
    Udat  = f_temp_data isa AbstractMatrix ? f_temp_data : reduce(hcat, f_temp_data)
    size(Umod) == size(Udat) || return -Inf

    # --- Likelihood with node-specific σ’s (constant over time per node) ---
    N_cols      = size(Umod, 2)
    σ_full      = repeat(σ_t_vec, N_cols)                # length 8N, matches vec order
    loglik      = sum(logpdf.(Normal.(vec(Umod), σ_full), vec(Udat)))

    # --- Priors (mirror your joint code’s spirit; adjust if needed) ---
    logprior = 0.0

    # 19 resistances (R): truncated Normal(2.0, 2.0) on (0,5)
    # 8 capacities (C):   truncated Normal(2000, 2000) on (0,6000)
    for i in 1:19
        logprior += logpdf(truncated(Normal(0.5, 2), 0.05, 8.0), P_therm[i]) # was truncated(Normal(0.5, 2), 0.05, 5.0)
    end
    for i in 1:8
        logprior += logpdf(truncated(Normal(3000.0, 2000.0), 500, 15_000.0), P_therm[19 + i]) # was truncated(Normal(3000.0, 2000.0), 1_000, 10_000.0)
    end

    # Initial temperatures (8): Normal(μ0 = first data col, σ = σ_i)
    for i in 1:8
        μ0 = T_init_fixed[i]
        σi = σ_t_vec[i]
        logprior += logpdf(Normal(μ0, σi), C_temp0[i])
    end

    # σ_T priors (positive): truncated Normal(0.5, 0.3) on (0, ∞)
    prior_σ_t = truncated(Normal(0.1, 0.06), 0.0, 0.5)
    for i in 1:8
        logprior += logpdf(prior_σ_t, σ_t_vec[i])
    end

    return loglik + (isfinite(logprior) ? logprior : -Inf)
end