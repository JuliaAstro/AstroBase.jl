using Base.Iterators: product
using Distributed

@everywhere using AstroBase

μ = 3.986004418e14

procs = nprocs()
n = 4
abort = false

if !isempty(ARGS)
    try
        global n
        n = something(tryparse.(Int, ARGS)...)
    catch
        # pass
    end

    abort = "abort" in ARGS
else
    if procs > 1
        n = (procs - 1) * 2
    end
end

@everywhere function keplerian_roundtrip(ele₀, μ; abort=false)
    ele₁ = keplerian(cartesian(ele₀..., μ)..., μ)
    failed = !(isapprox(ele₁[1], ele₀[1], rtol=1.0) &&
               isapprox(ele₁[2], ele₀[2], atol=sqrt(eps())) &&
               isapprox(ele₁[3], ele₀[3], atol=sqrt(eps())) &&
               isapprox(ele₁[4], ele₀[4], atol=sqrt(eps())) &&
               isapprox(ele₁[5], ele₀[5], atol=sqrt(eps())) &&
               isapprox(abs(ele₁[6]), abs(ele₀[6]), atol=sqrt(eps())))
    if failed && abort
        error("Results do not match!\n",
              "a₁: ", rpad(ele₁[1], 24), "a₀:", ele₀[1], "\n",
              "e₁: ", rpad(ele₁[2], 24), "e₀:", ele₀[2], "\n",
              "i₁: ", rpad(ele₁[3], 24), "i₀:", ele₀[3], "\n",
              "Ω₁: ", rpad(ele₁[4], 24), "Ω₀:", ele₀[4], "\n",
              "ω₁: ", rpad(ele₁[5], 24), "ω₀:", ele₀[5], "\n",
              "ν₁: ", rpad(ele₁[6], 24), "ν₀:", ele₀[6])
    end
    failed, ele₁, ele₀
end

function run_tests(rngs...)
    iter = product(rngs...)
    total = length(iter)
    println("Testing $total combinations with $procs processes.")

    results = pmap(x -> keplerian_roundtrip(x, μ; abort=abort), iter)
    failed = reduce(+, map(first, results))
    passed = total - failed
    percent_failed = round(Int, failed / total * 100.0)
    percent_passed = round(Int, passed / total * 100.0)
    println("Passed $passed ($percent_passed%), failed $failed ($percent_failed%).")
end

println("Elliptic Orbits")
println("===============")

a_rng = range(7e6, 9e8, length=n)
e_rng = range(0.001, 0.999, length=n)
i_rng = range(0.001, 0.999π, length=n)
Ω_rng = range(0.001, 1.999π, length=n)
ω_rng = range(0.001, 1.999π, length=n)
ν_rng = range(-π, π, length=n)

t = @elapsed run_tests(a_rng, e_rng, i_rng, Ω_rng, ω_rng, ν_rng)
println("In $t seconds.\n")

println("Hyperbolic Orbits")
println("=================")

a_rng = range(-2e6, -9e8, length=n)
e_rng = range(1.001, 5, length=n)
i_rng = range(0.001, 0.999π, length=n)
Ω_rng = range(0.001, 1.999π, length=n)
ω_rng = range(0.001, 1.999π, length=n)
ν_rng = range(-0.999π, 0.999π, length=n)

t = @elapsed run_tests(a_rng, e_rng, i_rng, Ω_rng, ω_rng, ν_rng)
println("In $t seconds.\n")

