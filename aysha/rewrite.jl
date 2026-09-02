
code = read("1D_v42.jl", String)

old_logic = "    scale = absorbed_power_scale_v42(irradiance, p)\r\n    Qcore_solar = ETA_ABS_FIXED_v42 * scale * Q_aperture\r\n    Qperim_solar = scale * Q_spill_total * perimeter_spillover_capture_v42(p)"
new_logic = "    M = absorbed_power_scale_v42(irradiance, p)\r\n    Q_delivered = M * Q_aperture\r\n    chi = core_source_fraction_v42(p)\r\n    Qcore_solar = chi * Q_delivered\r\n    Qperim_solar = (1.0 - chi) * Q_delivered"
code = replace(code, old_logic => new_logic)

old_logic2 = "    scale = absorbed_power_scale_v42(irradiance, p)\n    Qcore_solar = ETA_ABS_FIXED_v42 * scale * Q_aperture\n    Qperim_solar = scale * Q_spill_total * perimeter_spillover_capture_v42(p)"
code = replace(code, old_logic2 => new_logic)

code = replace(code, "perimeter_spillover_capture_v42(p)" => "core_source_fraction_v42(p)")

lines = split(code, "\n")
in_lb = false
for i in 1:length(lines)
    global in_lb
    line = lines[i]
    if occursin("lb_full_v42 = [", line)
        in_lb = true
    elseif in_lb && occursin("]", line)
        in_lb = false
    elseif in_lb
        if occursin("150.0,", line)
            lines[i] = replace(line, "150.0," => "100.0,")
        end
        if occursin("80.0,", line)
            lines[i] = replace(line, "80.0," => "50.0,")
        end
    end
end

write("1D_v42.jl", join(lines, "\n"))

