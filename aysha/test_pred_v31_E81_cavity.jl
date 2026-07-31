include("1D_v31.jl"); model, result, exp = solve_case_v31(pnew_v31, "E81"); println("Tcavity Model v31 E81: ", model[end, 7], " | T2 Exp E81: ", exp[end, 7]);
