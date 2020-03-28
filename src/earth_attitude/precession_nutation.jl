export precession_nutation

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

