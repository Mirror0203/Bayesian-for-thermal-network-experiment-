using CSV, DataFrames, Dates, Plots, Colors
using Statistics

#filepath = "D:/OneDrive/桌面/OneDrive - TU Eindhoven/Desktop/Experiment/20251024 no airflow/CO2.csv" # Read the CSV file 
filepath = "C:/Users/A/OneDrive - TU Eindhoven/Desktop/Experiment/20251024 no airflow/CO2.csv"
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



function plot_CO2_absolute_time(df_df)
    plot_interval = 10 # minutes
    rooms = [:Atm, :H1, :H2, :A, :B, :C, :D, :E, :F]

    # 9 distinguishable colors (Okabe–Ito palette + 1 extra)
    colors = [
        colorant"#000000",   # black
        colorant"#E69F00",  # orange
        colorant"#56B4E9",  # sky blue
        colorant"#009E73",  # green
        colorant"#F0E442",  # yellow
        colorant"#0072B2",  # dark blue
        colorant"#D55E00",  # reddish orange
        colorant"#CC79A7",  # pink/purple
        colorant"#999999"  # grey
        
    ]

    t0 = df_df.date_time[1]
    rel_min = Float64.((df_df.date_time .- t0) ./ Minute(plot_interval))

    plt = plot(
        xlabel = "*$plot_interval minutes from start",
        ylabel = "CO₂ concentration (ppm)",
        size = (800, 450),
        legend = :outerright,
        xlim = (0, maximum(rel_min)),
        xticks = 0:ceil(Int, maximum(rel_min))
    )

    for (i, r) in enumerate(rooms)
        plot!(plt, rel_min, df_df[!, r], label = String(r), color = colors[i], lw=2)
    end

    return plt
end


plot_CO2_absolute_time(df_tair_cal)
plot_CO2_absolute_time(df_tair_cal)