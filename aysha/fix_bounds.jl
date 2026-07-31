lines = readlines("1D_v35.jl")
for i in 1:length(lines)
    if occursin("12000.0, # p[12] upper bound", lines[i])
        lines[i] = "        2000.0, # p[12] upper bound"
    end
    if occursin("50.0,", lines[i]) && occursin("pnew_v35[13]", lines[i+1])
        lines[i] = "        2000.0,"
    end
end

in_lb = false
in_ub = false
for i in 1:length(lines)
    if occursin("lb_full_v35 = [", lines[i])
        global in_lb = true
    end
    if occursin("ub_full_v35 = [", lines[i])
        global in_ub = true
    end
    if in_lb && occursin("0.0,", lines[i]) && occursin("0.10,", lines[i+2])
        lines[i] = "        300.0,"
        global in_lb = false
    end
    if in_ub && occursin("300.0,", lines[i]) && occursin("20.0,", lines[i+2])
        # already 300
        global in_ub = false
    end
end

write("1D_v35.jl", join(lines, "\n") * "\n")
