using AstroBase
using Base.Iterators: product
using Base.Threads: nthreads, @threads

μ = 3.986004418e14
n = 2

a_rng = range(7e6, 9e8, length=n)
e_rng = range(0.001, 0.999, length=n)
i_rng = range(0.001, 0.999π, length=n)
Ω_rng = range(0.001, 1.999π, length=n)
ω_rng = range(0.001, 1.999π, length=n)
ν_rng = range(-π, π, length=n)
iter = collect(product(a_rng, e_rng, i_rng, Ω_rng, ω_rng, ν_rng))
println("Testing $(length(iter)) combinations with $(nthreads()) threads.")
#= @threads for i in eachindex(iter) =#
#=     ele₀ = iter[i] =#
#=     a, e, i, Ω, ω, ν = ele₀ =#
#=     ele₁ = keplerian(cartesian(a, e, i, Ω, ω, ν, μ)..., μ) =#
#=     if !(isapprox(ele₁[1], ele₀[1], rtol=1.0) && =#
#=          isapprox(ele₁[2], ele₀[2], atol=sqrt(eps())) && =#
#=          isapprox(ele₁[3], ele₀[3], atol=sqrt(eps())) && =#
#=          isapprox(ele₁[4], ele₀[4], atol=sqrt(eps())) && =#
#=          isapprox(ele₁[5], ele₀[5], atol=sqrt(eps())) && =#
#=          isapprox(abs(ele₁[6]), abs(ele₀[6]), atol=sqrt(eps()))) =#
#=         println("a: ", ele₁[1], " ", ele₀[1]) =#
#=         println("e: ", ele₁[2], " ", ele₀[2]) =#
#=         println("i: ", ele₁[3], " ", ele₀[3]) =#
#=         println("Ω: ", ele₁[4], " ", ele₀[4]) =#
#=         println("ω: ", ele₁[5], " ", ele₀[5]) =#
#=         println("ν: ", ele₁[6], " ", ele₀[6]) =#
#=         break =#
#=     end =#
#= end =#

a_rng = range(-2e6, -9e8, length=n)
e_rng = range(1.001, 5, length=n)
i_rng = range(0.001, 0.999π, length=n)
Ω_rng = range(0.001, 1.999π, length=n)
ω_rng = range(0.001, 1.999π, length=n)
ν_rng = range(-π, π, length=n)
iter = collect(product(a_rng, e_rng, i_rng, Ω_rng, ω_rng, ν_rng))
println("Testing $(length(iter)) combinations with $(nthreads()) threads.")
for i in eachindex(iter)
    ele₀ = iter[i]
    a, e, i, Ω, ω, ν = ele₀
    ele₁ = keplerian(cartesian(a, e, i, Ω, ω, ν, μ)..., μ)
    #= ele₁ = AstroBase.TwoBody.keplerian_orekit(AstroBase.TwoBody.cartesian_orekit(a, e, i, Ω, ω, ν, μ)..., μ) =#
    if !(isapprox(ele₁[1], ele₀[1], rtol=1.0) &&
         isapprox(ele₁[2], ele₀[2], atol=sqrt(eps())) &&
         isapprox(ele₁[3], ele₀[3], atol=sqrt(eps())) &&
         isapprox(ele₁[4], ele₀[4], atol=sqrt(eps())) &&
         isapprox(ele₁[5], ele₀[5], atol=sqrt(eps())) &&
         isapprox(abs(ele₁[6]), abs(ele₀[6]), atol=sqrt(eps())))
        Core.println("a₁: ", ele₁[1]," a₀: ", ele₀[1])
        Core.println("e₁: ", ele₁[2]," e₀: ", ele₀[2])
        Core.println("i₁: ", ele₁[3]," i₀: ", ele₀[3])
        Core.println("Ω₁: ", ele₁[4]," Ω₀: ", ele₀[4])
        Core.println("ω₁: ", ele₁[5]," ω₀: ", ele₀[5])
        Core.println("ν₁: ", ele₁[6]," ν₀: ", ele₀[6])
        break
    end
end
