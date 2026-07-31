include("1D_v31.jl"); model, result, exp = solve_case_v31(pnew_v31, "E81"); println("T8 Model v31: ", model[end, 1], " | T8 Exp: ", exp[end, 1]);
