# main_thermonly.jl — Temperature-only thermal inference (no airflow/occupancy inference)
#
# This script mirrors your air-only pipeline, but for the THERMAL side:
#   • infers 27 thermal parameters (19 R + 8 C)
#   • infers 8 initial temperatures for each window
#   • infers 8 node-wise σ_T
# Airflow/occupancy (13 params) are FIXED here and fed as P_air_fixed.
#
# NOTE:
#   1) By default below we set P_air_fixed to zeros as a placeholder.
#      Replace it with posterior means from your air-only run for best results.
#   2) Ensure models.jl defines solutionForward_thermal(P_therm_full, C_temp0, tspan)
#      where P_therm_full = vcat(P_therm, P_air_fixed).

using Random, LinearAlgebra, Statistics
using RobustAdaptiveMetropolisSampler, CSV, DataFrames, Dates
using Measures, Plots, Colors

include("models.jl")
include("inference_thermonly.jl")

filepath = "D:/OneDrive/桌面/OneDrive - TU Eindhoven/Desktop/Experiment/20251024 no airflow/CO2.csv" # Read the CSV file 
#filepath = "C:/Users/A/OneDrive - TU Eindhoven/Desktop/Experiment/20251024 no airflow/CO2.csv"
df = CSV.read(filepath, DataFrame)
df_CO2 = unstack(df, :date_time, :sensor_id, :CO2)
df_tair = unstack(df, :date_time, :sensor_id, :t_air)

# Map sensor IDs to room names
name_map = Dict(
    :ID3983 => :H2,
    :ID3982 => :E,
    :ID3978 => :H1,
    :ID3987 => :F,
    :ID3985 => :A,
    :ID3979 => :B,
    :ID3980 => :D,
    :ID3981 => :C,
    :ID3986 => :Atm
)

# Rename columns in both wide tables
rename!(df_CO2,  name_map)
rename!(df_tair, name_map)

t_start = Time("11:10:00")
t_end   = Time("11:20:00")

# === Filter to that period ===
df_calib_CO2  = filter(row -> t_start ≤ row.date_time ≤ t_end, df_CO2)
df_calib_tair = filter(row -> t_start ≤ row.date_time ≤ t_end, df_tair)

# === Helper function to compute mean per room (ignoring :date_time) ===
function room_means(df)
    means = Dict{Symbol, Float64}()
    for col in names(df)
        coldata = df[!, col]
        # Only take columns that contain numbers (not Time)
        if eltype(coldata) <: Union{Missing, Number}
            means[Symbol(col)] = mean(skipmissing(coldata))
        end
    end
    return means
end

# === Compute mean per room ===
mean_CO2_by_room  = room_means(df_calib_CO2)
mean_tair_by_room = room_means(df_calib_tair)

# === Compute grand means ===
grand_mean_CO2  = mean(values(mean_CO2_by_room))
grand_mean_tair = mean(values(mean_tair_by_room))

# === Compute offsets ===
co2_offsets_room = Dict(room => mean_CO2_by_room[room] - grand_mean_CO2 for room in keys(mean_CO2_by_room))
t_offsets_room   = Dict(room => mean_tair_by_room[room] - grand_mean_tair for room in keys(mean_tair_by_room))

# Make deep copies so original data remain intact
df_CO2_cal  = deepcopy(df_CO2)
df_tair_cal = deepcopy(df_tair)

for room in keys(co2_offsets_room)
    df_CO2_cal[!, room] .= df_CO2[!, room] .- co2_offsets_room[room]
end

for room in keys(t_offsets_room)
    df_tair_cal[!, room] .= df_tair[!, room] .- t_offsets_room[room]
end

# Choose model room ordering
const MODEL_ROOMS = [:A, :B, :C, :D, :E, :F, :H1, :H2]

# Helper to produce 8×N matrix in model order
function df_temp_to_matrix_8(df_temp)::Matrix{Float64}
    cols_needed = MODEL_ROOMS
    missing_cols = setdiff(cols_needed, propertynames(df_temp))
    if !isempty(missing_cols)
        error("Missing columns in df_tair: $(missing_cols). Did mapping/renaming run?")
    end
    M = reduce(hcat, (Float64.(df_temp[!, r]) for r in cols_needed))
    return permutedims(M)  # 8×N
end

# Drop :Atm and map to model rooms for temperature
select!(df_tair_cal, [:date_time; :A; :B; :C; :D; :E; :F; :H1; :H2])

# Optional: define a time window (copy your choices if needed)
# (You can comment this out to run full series.)
# start_t = Time(11,42,0)
# end_t   = Time(12,0,0)
# mask = (df_tair.date_time .>= start_t) .& (df_tair.date_time .< end_t)
# df_tair = df_tair[mask, :]

f_temp_real  = df_temp_to_matrix_8(df_tair_cal)
f_temp_noise = f_temp_real       # if you separate noise and real, adjust here
N_real = size(f_temp_real, 2)

# -------------------- Fixed airflow/occupancy block (13 params) --------------------
# IMPORTANT: Replace this placeholder with your air-only posterior means for best results.
# Order must match the airflow block your thermal model expects: [8×num_ppl, 5×AF_Atm_*]
P_air_fixed = vcat(zeros(8), zeros(5))  # length 13

# -------------------- Thermal-only inference setup --------------------
# Parameter vector P length 43 = [27 thermal ; 8 T0 ; 8 σT]
# Provide a sensible initial guess (adjust as you like).
R_init = fill(0.53, 19)           # near prior mean
C_init = fill(3_230.0, 8)        # near prior mean
T0_init = vec(f_temp_real[:, 1]) # start at first observed temps
σT_init = fill(0.1, 8)           # near prior mean

P0 = vcat(R_init, C_init, T0_init, σT_init)

# Proposal scales (diagonal RAM initial widths) – tune if needed
sigma_init = vcat(fill(0.10, 19),   # R
                  fill(200.0, 8),   # C
                  fill(0.15, 8),    # T0
                  fill(0.04, 8))    # σT

# Plot helper (optional)
include("visualize.jl")

# -------------------- Sliding-window MCMC --------------------
mkpath("Figures/thermonly_windows")

time_window_size = 200

for i in 100:(N_real - time_window_size)
    window_start = i
    window_end   = window_start + time_window_size
    f_temp_window = f_temp_noise[:, window_start:window_end]

    # Fixed initial state prior mean for the window: first column of the window
    T_init_fixed = f_temp_window[:, 1]

    # Posterior: P = [27 thermal ; 8 T0 ; 8 σT] (length 43)
    function logPosterior_thermonly(P)
        evaluateLogPdf_thermonly(P, f_temp_window, time_window_size, T_init_fixed, P_air_fixed)
    end

    Random.seed!(0)
    n_samples   = 200_000
    output      = RAM_sample(logPosterior_thermonly, P0, sigma_init, n_samples)
    whole_chain = output.chain
    sampling_chain = whole_chain[end-50_000+1:end, :]
    post_mean   = vec(mean(sampling_chain, dims=1))  # length 43

    # Slices
    P_therm_mean = post_mean[1:27]
    T0_mean      = post_mean[28:35]
    σT_mean      = post_mean[36:43]

    # Update init for next window with the whole vector
    P0 = post_mean
    println("Posterior mean (thermal-only):",
            "\n  R,C = ", P_therm_mean,
            "\n  T0  = ", T0_mean,
            "\n  σT  = ", σT_mean)

    # Plot with node-wise ribbons
    plt = plot_thermal_inference_fullseries(
        f_temp_noise, window_start, window_end,
        P_therm_mean, T0_mean, σT_mean
    )
    savefig(plt, "Figures/thermonly_windows/thermonly_window_$(i).png")

    # CSV log: (index, 27 thermal, 8 T0, 8 σT)
    open("Figures/thermonly_windows/thermonly_log.csv", "a") do f
        write(f, join([i; P_therm_mean; T0_mean; σT_mean], ",") * "\n")
    end
end