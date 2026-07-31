content = read("1D_v36.jl", String)

content = replace(content, "function active_flow_fraction_v36(reynolds, p)\n    phi0 = p[4]\n    m = max(p[5], 0.0)\n    ratio = max(reynolds, eps(Float64)) / REYNOLDS_REFERENCE_v36\n    return clamp(phi0 * ratio^m, 0.1, 1.0)\nend" => "function active_flow_fraction_v36(Tcore_avg, p)\n    return clamp(p[4] * (Tcore_avg / 295.0)^(-p[5]), 0.1, 1.0)\nend")

content = replace(content, "function receiver_active_flow_fraction_v36(flow_lpm, Tin, p)\n    flow_lpm <= 1e-12 && return 1.0\n    mdot = m_dot_standard_v36(flow_lpm)\n    mu = mu_f_f(Tin)\n    reynolds = mdot * Dh / (A_flow * mu)\n    return active_flow_fraction_v36(reynolds, p)\nend" => "function receiver_active_flow_fraction_v36(flow_lpm, Tcore_avg, p)\n    flow_lpm <= 1e-12 && return 1.0\n    return active_flow_fraction_v36(Tcore_avg, p)\nend")

content = replace(content, "phi_act = receiver_active_flow_fraction_v36(flow, Tin, p)" => "Tcore_avg = sum(Tcore) / length(Tcore)\n    phi_act = receiver_active_flow_fraction_v36(flow, Tcore_avg, p)")

content = replace(content, "        phi_act = receiver_active_flow_fraction_v36(flow, Tin, p)" => "        Tcore_avg = sum(Tcore_now) / length(Tcore_now)\n        phi_act = receiver_active_flow_fraction_v36(flow, Tcore_avg, p)")

content = replace(content, "        0.05,\n        0.0,\n        pnew_v36[6]," => "        0.1,\n        0.0,\n        pnew_v36[6],")

write("1D_v36.jl", content)

