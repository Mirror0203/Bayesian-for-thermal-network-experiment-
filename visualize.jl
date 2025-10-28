
function plot_data(df1)
    rooms = [:A, :B, :C, :D, :E, :F,:H1, :H2]
    plt = plot(
        xlabel = "Time",
        ylabel = "CO2 (ppm)",#"T (C)", 
        title  =  "CO2 levels — all rooms",#"Tmeperature levels — all rooms",
        legend = :outerright,
        size = (900, 400),
        bottom_margin = 5mm,
        left_margin = 5mm
    )
    for r in rooms
        plot!(plt, df1.date_time, df1[!, r], label = String(r))
    end

    return(plt)
end

# --- in main_thermalonly.jl ---
function plot_thermal_inference_fullseries(
    f_temp_full::AbstractMatrix,
    window_start::Int,
    window_end::Int,
    P_therm_mean::AbstractVector,  # length 27
    T0_mean::AbstractVector,       # length 8 (kept for reference)
    σT_vec::AbstractVector;        # length 8
    labels = ["Room A" "Room B" "Room C" "Room D" "Room E" "Room F" "Room H1" "Room H2"],
    use_observed_u0::Bool = true
)
    N = size(f_temp_full, 2)
    ts_all    = collect(0.0:(N-1))
    tspan_win = (window_start - 1.0, min(window_end - 1.0, N - 1.0))

    # Initial condition for the first segment
    u0_plot = use_observed_u0 ? vec(f_temp_full[:, window_start]) : T0_mean

    # Forward with posterior mean on the chosen window
    P_therm_full = vcat(P_therm_mean, P_air_fixed)  # keep airflow block fixed/zeroed outside
    sol1 = solutionForward_thermal(P_therm_full, u0_plot, tspan_win)

    # Continue the trajectory to the end of the full record
    if window_end < N
        u0_tail = sol1.u[end]
        tspan_tail = (window_end - 1.0, N - 1.0)
        sol2 = solutionForward_thermal(P_therm_full, u0_tail, tspan_tail)

        # Concatenate time and states (drop duplicate junction point)
        ts_mean = vcat(sol1.t, sol2.t[2:end])
        Y_mean  = vcat(permutedims(reduce(hcat, sol1.u)),
                       permutedims(reduce(hcat, sol2.u))[2:end, :])
    else
        ts_mean = sol1.t
        Y_mean  = permutedims(reduce(hcat, sol1.u))
    end

    fT   = permutedims(f_temp_full)
    ymin = minimum(fT) - 1
    ymax = maximum(fT) + 1

    plt = plot(size=(1500, 600),
        xlabel="Time (steps)", ylabel="Temperature (°C)",
        legend=:outerright, left_margin=12mm, right_margin=18mm,
        bottom_margin=14mm, top_margin=6mm, ylim=(ymin, ymax))

    # Mark the inference window
    vline!(plt, [window_start-1, window_end-1]; linestyle=:dash, linewidth=3, color=:black, label="")

    # Constant ±2σ ribbons for the full posterior-mean curve
    @assert length(σT_vec) == 8
    ribbon_full = hcat([fill(2*σT_vec[j], length(ts_mean)) for j in 1:8]...)

    plot!(plt, ts_mean, Y_mean;
          ribbon=ribbon_full, linestyle=:dash, linewidth=2, color=:black,
          fillcolor=:gray80, fillalpha=0.35, label="Model ±2σ")

    # downsample observed points for readability
    skip_all = max(1, round(Int, length(ts_all) / 400))
    idx_all  = 1:skip_all:length(ts_all)

    scatter!(plt, ts_all[idx_all], fT[idx_all, :];
        ms=2, markerstrokewidth=0, alpha=0.6,
        colors = [
            colorant"#E69F00",  # Room A
            colorant"#56B4E9",  # Room B
            colorant"#009E73",  # Room C
            colorant"#F0E442",  # Room D
            colorant"#0072B2",  # Room E
            colorant"#D55E00",  # Room F
            colorant"#CC79A7",  # Room H1
            colorant"#999999"   # Room H2
        ],
        label = labels)

    return plt
end
