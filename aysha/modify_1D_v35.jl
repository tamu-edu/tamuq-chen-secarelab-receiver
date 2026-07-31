str = read("1D_v35.jl", String)
str = replace(str, "_v34" => "_v35")
str = replace(str, "v34" => "v35")

# 1. Conductivity
old_cond = r"perimeter_axial_conductivity_v35\(T, p\) =[\r\n\s]+max\(p\[12\], 0\.0\) \* \(property_temperature\(T\) / 900\.0\)\^3"
new_cond = "perimeter_axial_conductivity_v35(T, p) = max(p[12], 0.0)"
str = replace(str, old_cond => new_cond)

# 2. Upper bound for p[12]
str = replace(str, r"50\.0,([\r\n\s]+pnew_v35\[13\],)" => s"2000.0,\1")

# 3. lb_full_v35 array
old_lb = "    lb_full_v35 = [\n        0.01,\n        0.0,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        0.5,\n        0.5,\n        0.5,\n        2.0,\n        150.0,\n        0.0,\n        pnew_v35[13],\n        0.0,\n        0.0,"
new_lb = "    lb_full_v35 = [\n        0.01,\n        0.0,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        0.5,\n        0.5,\n        0.5,\n        2.0,\n        150.0,\n        0.0,\n        pnew_v35[13],\n        0.0,\n        300.0, # p[15] beta_perim fixed"
str = replace(str, old_lb => new_lb)

# 3. ub_full_v35 array
old_ub = "    ub_full_v35 = [\n        25.00,\n        0.60,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        1.5,\n        1.5,\n        1.5,\n        80.0,\n        230.0,\n        2000.0,\n        pnew_v35[13],\n        1.00,\n        300.0,"
new_ub = "    ub_full_v35 = [\n        25.00,\n        0.60,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        1.5,\n        1.5,\n        1.5,\n        80.0,\n        230.0,\n        2000.0,\n        pnew_v35[13],\n        1.00,\n        300.0,"
# wait, wait! The p[12] upper bound replacement will change the string to 2000.0, let's just make sure.
write("1D_v35.jl", str)
