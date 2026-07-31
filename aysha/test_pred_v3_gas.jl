include("1D_v3.jl"); model, result, exp = solve_case_v3(pnew, "E81"); println("T3 Model v3: ", model[end, 4], " | T3 Exp: ", exp[end, 4]);
