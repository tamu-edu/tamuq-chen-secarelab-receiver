str = read("1D_v35.jl", String)

# Replace the array elements directly
# lb_full_v35
old_lb = "    lb_full_v35 = [\n        0.01,\n        0.0,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        0.5,\n        0.5,\n        0.5,\n        2.0,\n        150.0,\n        0.0,\n        pnew_v35[13],\n        0.0,\n        0.0,"
new_lb = "    lb_full_v35 = [\n        0.01,\n        0.0,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        0.5,\n        0.5,\n        0.5,\n        2.0,\n        150.0,\n        0.0,\n        pnew_v35[13],\n        0.0,\n        300.0,"
str = replace(str, old_lb => new_lb)

# ub_full_v35
old_ub = "    ub_full_v35 = [\n        25.00,\n        0.60,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        1.5,\n        1.5,\n        1.5,\n        80.0,\n        230.0,\n        2000.0,\n        pnew_v35[13],\n        1.00,\n        300.0,"
new_ub = "    ub_full_v35 = [\n        25.00,\n        0.60,\n        pnew_v35[3],\n        pnew_v35[4],\n        pnew_v35[5],\n        pnew_v35[6],\n        1.5,\n        1.5,\n        1.5,\n        80.0,\n        230.0,\n        2000.0,\n        pnew_v35[13],\n        1.00,\n        300.0,"
# wait, wait! The p[12] upper bound is STILL 50.0 in the file! I checked the output!
# Let's just do it manually with regexes that actually work.
