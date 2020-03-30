#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
export precession_nutation, precession_nutation_matrix

function precession_nutation(::IAU2000, ep::Epoch, δψ, δϵ)
    # IAU 2000 precession-rate adjustments.
    δψ_pr, δϵ_pr = precession(iau2000, ep)
    # Mean obliquity, consistent with IAU 2000 precession-nutation.
    ϵ = obliquity(iau1980, ep) + δϵ_pr
    # Frame bias and precession matrices and their product.
    rb, rp, rbp = bias_precession_matrix(iau2000, ep)
    # Nutation matrix
    rn = nutation_matrix(ϵ, δψ, δϵ)
    return ϵ, rb, rp, rbp, rn, rn * rbp
end

function precession_nutation(iau::IAU2000Model, ep::Epoch)
    n = nutation(iau, ep)
    pn = precession_nutation(iau2000, ep, n...)
    return n..., pn...
end

function precession_nutation_matrix(iau::IAU2000Model, ep::Epoch)
    pn = precession_nutation(iau, ep)
    return pn[end]
end

function nutation_matrix(iau::IAU2000Model, ep::Epoch)
    pn = precession_nutation(iau, ep)
    return pn[end-1]
end

function precession_nutation_matrix(::IAU2006A, ep::Epoch)
    γ, ϕ, ψ, ϵ = fukushima_williams(iau2006, ep)
    δψ, δϵ = nutation(iau2006a, ep)
    return fukushima_williams_matrix(γ, ϕ, ψ + δψ, ϵ + δϵ)
end

