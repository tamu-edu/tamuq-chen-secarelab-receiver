include("1D_v3.jl"); model, result, exp = solve_case_v3(pnew, "E67"); println("T8 Model v3 E67: ", model[end, 1], " | T8 Exp E67: ", exp[end, 1]);
