content = read(\"1D_v38.jl\", String)

# Replace gas_profile_v38!
old_gas_str = \"function gas_profile_v38!(Tg, Qgas, hcell, Qgas_rear, hrear,\\n                          Tcore, Ttube, time, p, operating, z, dx, z_rear, dx_rear)\"
new_gas_str = \"function gas_profile_v38!(Tg, Qgas, hcell, Qgas_rear, hrear,\\n                          Tcore, Tperim, Ttube, time, p, operating, z, dx, z_rear, dx_rear, Tgas_perim, Qgas_perim)\"
content = replace(content, old_gas_str => new_gas_str)

# We also need to replace the body of gas_profile_v38!.
# Let's use a regex to replace the entire body.
# Wait, it's safer to just provide the full new function.
