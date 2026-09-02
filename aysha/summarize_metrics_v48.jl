using CSV, DataFrames, Statistics

df = CSV.read("summaries/1D_v48/analysis_results_fitted_energy_accounting_1D_v48.csv", DataFrame)

println("--- Summary by Phase ---")
for gdf in groupby(df, :phase)
    p = gdf.phase[1]
    r = round(mean(gdf.rmse_K), digits=2)
    s = round(mean(abs.(gdf.steady_error_K)), digits=2)
    t = round(mean(abs.(gdf.t90_error_s)), digits=2)
    println("Phase $p: Mean RMSE = $r K, Mean |Steady Err| = $s K, Mean |t90 Err| = $t s")
end

println("\n--- Summary by Sensor (Heating) ---")
df_heat = filter(r -> r.phase == "heating", df)
for gdf in groupby(df_heat, :sensor)
    sn = gdf.sensor[1]
    r = round(mean(gdf.rmse_K), digits=2)
    s = round(mean(abs.(gdf.steady_error_K)), digits=2)
    println("Sensor $sn: Mean RMSE = $r K, Mean |Steady Err| = $s K")
end

println("\n--- Summary by Sensor (Cooling) ---")
df_cool = filter(r -> r.phase == "cooling", df)
for gdf in groupby(df_cool, :sensor)
    sn = gdf.sensor[1]
    r = round(mean(gdf.rmse_K), digits=2)
    s = round(mean(abs.(gdf.steady_error_K)), digits=2)
    println("Sensor $sn: Mean RMSE = $r K, Mean |Steady Err| = $s K")
end
