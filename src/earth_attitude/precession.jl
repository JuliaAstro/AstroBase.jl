using ..Time: Epoch, TT, centuries, julian_period
using ..Util: sec2rad

export
    bias,
    bias_precession_matrix,
    fukushima_williams,
    fukushima_williams_matrix,
    precession

function bias(::IAU2000)
    δψ_b = sec2rad(-0.041775)
    δϵ_b = sec2rad(-0.0068192)
    δra₀ = sec2rad(-0.0146)
    return δψ_b, δϵ_b, δra₀
end

function precession(::IAU2000, ep::Epoch; scale=TT)
    t = julian_period(ep; scale=scale, unit=centuries, raw=true)
    precor = sec2rad(-0.29965)
    oblcor = sec2rad(-0.02524)
    return precor * t, oblcor * t
end

function fukushima_williams(::IAU2006, ep::Epoch)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)

    γ = @evalpoly(
        t,
        -0.052928,
        10.556378,
        0.4932044,
        -0.00031238,
        -0.000002788,
        0.0000000260,
    ) |> sec2rad
    ϕ = @evalpoly(
        t,
        84381.412819,
        -46.811016,
        0.0511268,
        0.00053289,
        -0.000000440,
        -0.0000000176,
    ) |> sec2rad
    ψ = @evalpoly(
        t,
        -0.041775,
        5038.481484,
        1.5584175,
        -0.00018522,
        -0.000026452,
        -0.0000000148,
    ) |> sec2rad
    ϵ = obliquity(iau2006, ep)
    return γ, ϕ, ψ, ϵ
end

function fukushima_williams_matrix(γ, ϕ, ψ, ϵ)
    return compose_rotation(
        angleaxis_to_dcm(γ, [0.0, 0.0, 1.0]),
        angle_to_dcm(ϕ, -ψ, -ϵ, :XZX),
    )
end

function bias_precession_matrix(::IAU2000, ep::Epoch)
    t = julian_period(ep; scale=TT, unit=centuries, raw=true)

    ϵ₀ = sec2rad(84381.448)

    δψ_b, δϵ_b, δra₀ = bias(iau2000)

    ψ_a77 = sec2rad(@evalpoly(t, 0.0, 5038.7784, -1.07259, -0.001147))
    ω_a77  = ϵ₀ + sec2rad(@evalpoly(t, 0.0, 0.0, 0.05127, -0.007726))
    χ_a   = sec2rad(@evalpoly(t, 0.0, 10.5526, -2.38064, -0.001125))

    δψ_pr, δϵ_pr = precession(iau2000, ep)
    ψ_a = ψ_a77 + δψ_pr
    ω_a  = ω_a77 + δϵ_pr

    rb = angle_to_dcm(δra₀, δψ_b * sin(ϵ₀), -δϵ_b, :ZYX)
    rp = compose_rotation(
        angleaxis_to_dcm(ϵ₀, [1.0, 0.0, 0.0]),
        angle_to_dcm(-ψ_a, -ω_a, χ_a, :ZXZ),
    )

    return rb, rp, compose_rotation(rb, rp)
end

