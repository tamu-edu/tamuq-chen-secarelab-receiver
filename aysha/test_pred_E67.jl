include("1D_v31.jl"); model, result, exp = solve_case_v31(pnew_v31, "E67"); println("T8 Model v31 E67: ", model[end, 1], " | T8 Exp E67: ", exp[end, 1]);
