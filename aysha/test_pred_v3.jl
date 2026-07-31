include("1D_v3.jl"); model, result, exp = solve_case_v3(pnew, "E81"); println("T8 Model v3: ", model[end, 1], " | T8 Exp: ", exp[end, 1]);
