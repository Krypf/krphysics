function AccelSchwarzschild(v0::Vector{Float64}, x0::Vector{Float64}, args::Tuple)
    # theta = pi / 2
    rs = args[1]
    A = zero(x0)
    # t = x0[1];
    r = x0[2];
    # ph =  x0[3];
    tdiff = v0[1]; rdiff = v0[2]; phdiff =  v0[3];
    A[1] = rs * rdiff * tdiff / ((rs - r) * r)
    A[2] = -(2 * (rs-r)^2 * r^3 * (phdiff)^2 -rs*(rs-r)^2 * tdiff^2 + rs * r^2 * rdiff^2 ) / (2 * (rs-r) * r^3 )# theta = pi / 2 修正
    A[3] = -2 * rdiff * phdiff / r
    return A
end
#%%
function AccelKerr(v0::Vector{Float64}, x0::Vector{Float64}, args::Tuple)
    # theta = pi / 2
    a, M = args[1], args[2]
    A = zero(x0)
    t = x0[1]; r = x0[2]; ph =  x0[3];
    tdiff = v0[1]; rdiff = v0[2]; phdiff =  v0[3];
    A[1] = (2 * M * (a^3 * phdiff - a^2 * tdiff + 3 * a * r^2 * phdiff - r^2 * tdiff) * rdiff) / ((- 2 * M * r + a^2 + r^2) * r^2)

    A[2] = (M * (2 * a * phdiff - tdiff) * (- 2 * M * r + a^2 + r^2)^2 * r * tdiff + (- 2 * M * r + a^2 + r^2)^2 * (- a^2 * (2 * M * r - a^2 - r^2) - (a^2 + r^2)^2 + (- a^2 * (- M + r) + 2 * (a^2 + r^2) * r) * r) * (phdiff)^2 + (2 * M * r - a^2 - (M - r) * r - r^2) * r^4 * (rdiff)^2) / ((- 2 * M * r + a^2 + r^2) * r^5)

    A[3] = (2 * (M * a^2 * phdiff - M * a * tdiff + 2 * M * r^2 * phdiff - r^3 * phdiff) * rdiff) / ((- 2 * M * r + a^2 + r^2) * r^2)
    return A
end