macro test_elements(name, μ, ele, pos, vel)
    quote
        @testset $name begin
            μ = $μ
            a, e, i, Ω, ω, ν = $ele
            pos = $pos
            vel = $vel

            pos₁, vel₁ = cartesian(a, e, i, Ω, ω, ν, μ)
            a₁, e₁, i₁, Ω₁, ω₁, ν₁ = keplerian(pos, vel, μ)

            @testset "Cartesian" begin
                @testset for i = 1:3
                    @test pos₁[i] ≈ pos[i]
                    @test vel₁[i] ≈ vel[i]
                end
            end

            @testset "Keplerian" begin
                @test a ≈ a₁
                @test e ≈ e₁ atol=sqrt(eps())
                @test i ≈ i₁
                @test Ω ≈ Ω₁
                @test ω ≈ ω₁
                @test ν ≈ ν₁
            end
        end
    end
end

macro test_elements(name, μ, ele)
    quote
        @testset $name begin
            μ = $μ
            a, e, i, Ω, ω, ν = $ele

            pos₁, vel₁ = cartesian(a, e, i, Ω, ω, ν, μ)
            a₁, e₁, i₁, Ω₁, ω₁, ν₁ = keplerian(pos₁, vel₁, μ)

            @test a ≈ a₁
            @test e ≈ e₁ atol=sqrt(eps())
            @test i ≈ i₁
            @test Ω ≈ Ω₁
            @test ω ≈ ω₁
            @test ν ≈ ν₁
        end
    end
end

@testset "Two Body" begin
    @testset "Orbital Elements" begin
        @testset "Perifocal" begin
            p = 1.13880762905224e7
            ecc = 0.7311
            ν = 0.44369564302687126
            μ = 3.9860047e14

            pos = [6194863.12535486, 2944437.90016286, 0.0]
            vel = [-2539.71254827, 9668.69568539, 0.0]
            pos₁, vel₁ = perifocal(p, ecc, ν, μ)
            @testset for i = 1:3
                @test pos₁[i] ≈ pos[i]
                @test vel₁[i] ≈ vel[i]
            end
        end
        # Reference values from Orekit
        @test_elements("Elliptic",
                       3.9860047e14,
                       [24464560.0, 0.7311, 0.122138, 1.00681, 3.10686, 0.44369564302687126],
                       [-0.107622532467967e+07, -0.676589636432773e+07, -0.332308783350379e+06],
                       [0.935685775154103e+04, -0.331234775037644e+04, -0.118801577532701e+04])
        # Reference values from poliastro
        @test_elements("Circular",
                       3.986004418e14,
                       [6778136.6, 0.0, deg2rad(15), deg2rad(20), 0.0, deg2rad(30)],
                       [4396398.60746266, 5083838.45333733,  877155.42119322],
                       [-5797.06004014,  4716.60916063,  1718.86034246])
        @test_elements("Circular Orekit",
                       3.9860047e14,
                       # Orekit uses an argument of pericenter of 3.10686 here
                       [24464560.0, 0.0, 0.122138, 1.00681, 0.0, 0.048363])
        @test_elements("Hyperbolic Orekit",
                       3.9860047e14,
                       [-24464560.0, 1.7311, 0.122138, 1.00681,  3.10686, 0.12741601769795755])
        @test_elements("Equatorial",
                       3.9860047e14,
                       [24464560.0, 0.7311, 0.0, 0.0, 3.10686, 0.44369564302687126])
        @test_elements("Circular-Equatorial",
                       3.9860047e14,
                       [24464560.0, 0.0, 0.0, 0.0, 0.0, 0.44369564302687126])
    end
    @testset "Anomalies" begin
        anomalies = -π : 0.1π : π
        anomalies2π = 0.0 : 0.1π : 2π
        @testset "Elliptic" begin
            eccentricities = 0.1:0.1:0.9
            @testset "E <-> ν" begin
                @testset for E = anomalies, ecc = eccentricities
                    ν = transform(elliptic, eccentric_anomaly, true_anomaly, E, ecc)
                    E₁ = transform(elliptic, true_anomaly, eccentric_anomaly, ν, ecc)
                    @test E₁ ≈ E
                end
                @testset for E = anomalies2π, ecc = eccentricities
                    ν = transform(elliptic, eccentric_anomaly, true_anomaly, E, ecc)
                    E₁ = transform(elliptic, true_anomaly, eccentric_anomaly, ν, ecc)
                    @test E₁ ≈ E
                end
            end
            @testset "E <-> M" begin
                @testset for E = anomalies, ecc = eccentricities
                    M = transform(elliptic, eccentric_anomaly, mean_anomaly, E, ecc)
                    E₁ = transform(elliptic, mean_anomaly, eccentric_anomaly, M, ecc)
                    @test E₁ ≈ E
                end
                @testset for E = anomalies2π, ecc = eccentricities
                    M = transform(elliptic, eccentric_anomaly, mean_anomaly, E, ecc)
                    E₁ = transform(elliptic, mean_anomaly, eccentric_anomaly, M, ecc)
                    @test E₁ ≈ E
                end
            end
            @testset "M <-> ν" begin
                @testset for M = anomalies, ecc = eccentricities
                    ν = transform(elliptic, mean_anomaly, true_anomaly, M, ecc)
                    M₁ = transform(elliptic, true_anomaly, mean_anomaly, ν, ecc)
                    @test M₁ ≈ M
                end
                @testset for M = anomalies2π, ecc = eccentricities
                    ν = transform(elliptic, mean_anomaly, true_anomaly, M, ecc)
                    M₁ = transform(elliptic, true_anomaly, mean_anomaly, ν, ecc)
                    @test M₁ ≈ M
                end
            end
        end
        @testset "Hyperbolic" begin
            eccentricities = 1.1:0.1:10.0
            @testset "E <-> ν" begin
                @testset for E = anomalies, ecc = eccentricities
                    ν = transform(hyperbolic, eccentric_anomaly, true_anomaly, E, ecc)
                    E₁ = transform(hyperbolic, true_anomaly, eccentric_anomaly, ν, ecc)
                    @test E₁ ≈ E
                end
                @testset for E = anomalies2π, ecc = eccentricities
                    ν = transform(hyperbolic, eccentric_anomaly, true_anomaly, E, ecc)
                    E₁ = transform(hyperbolic, true_anomaly, eccentric_anomaly, ν, ecc)
                    @test E₁ ≈ E
                end
            end
            @testset "E <-> M" begin
                @testset for E = anomalies, ecc = eccentricities
                    M = transform(hyperbolic, eccentric_anomaly, mean_anomaly, E, ecc)
                    E₁ = transform(hyperbolic, mean_anomaly, eccentric_anomaly, M, ecc)
                    @test E₁ ≈ E
                end
                @testset for E = anomalies2π, ecc = eccentricities
                    M = transform(hyperbolic, eccentric_anomaly, mean_anomaly, E, ecc)
                    E₁ = transform(hyperbolic, mean_anomaly, eccentric_anomaly, M, ecc)
                    @test E₁ ≈ E
                end
            end
            @testset "M <-> ν" begin
                @testset for M = anomalies, ecc = eccentricities
                    ν = transform(hyperbolic, mean_anomaly, true_anomaly, M, ecc)
                    M₁ = transform(hyperbolic, true_anomaly, mean_anomaly, ν, ecc)
                    @test M₁ ≈ M
                end
                @testset for M = anomalies2π, ecc = eccentricities
                    ν = transform(hyperbolic, mean_anomaly, true_anomaly, M, ecc)
                    M₁ = transform(hyperbolic, true_anomaly, mean_anomaly, ν, ecc)
                    @test M₁ ≈ M
                end
            end
        end
        @testset "Parabolic" begin
            @testset "E <-> ν" begin
                @testset for E = anomalies
                    ν = transform(parabolic, eccentric_anomaly, true_anomaly, E, 1.0)
                    E₁ = transform(parabolic, true_anomaly, eccentric_anomaly, ν, 1.0)
                    @test E₁ ≈ E
                end
                @testset for E = anomalies2π
                    ν = transform(parabolic, eccentric_anomaly, true_anomaly, E, 1.0)
                    E₁ = transform(parabolic, true_anomaly, eccentric_anomaly, ν, 1.0)
                    @test E₁ ≈ E
                end
            end
            @testset "E <-> M" begin
                @testset for E = anomalies
                    M = transform(parabolic, eccentric_anomaly, mean_anomaly, E, 1.0)
                    E₁ = transform(parabolic, mean_anomaly, eccentric_anomaly, M, 1.0)
                    @test E₁ ≈ E
                end
                @testset for E = anomalies2π
                    M = transform(parabolic, eccentric_anomaly, mean_anomaly, E, 1.0)
                    E₁ = transform(parabolic, mean_anomaly, eccentric_anomaly, M, 1.0)
                    @test E₁ ≈ E
                end
            end
            @testset "M <-> ν" begin
                @testset for M = anomalies
                    ν = transform(parabolic, mean_anomaly, true_anomaly, M, 1.0)
                    M₁ = transform(parabolic, true_anomaly, mean_anomaly, ν, 1.0)
                    @test M₁ ≈ M
                end
                @testset for M = anomalies2π
                    ν = transform(parabolic, mean_anomaly, true_anomaly, M, 1.0)
                    M₁ = transform(parabolic, true_anomaly, mean_anomaly, ν, 1.0)
                    @test M₁ ≈ M
                end
            end
        end
    end
    @testset "Kepler" begin
        # Source: Vallado, Fundamentals of Astrodynamics and Applications, 4th edition, p. 94-95
        μ = 3.986004418e5
        r₀exp = [1131.340, -2282.343, 6672.423]
        v₀exp = [-5.64305, 4.30333, 2.42879]
        Δt = 40*60
        r₁exp = [-4219.7527, 4363.0292, -3958.7666]
        v₁exp = [3.689866, -1.916735, -6.112511]
        r₁, v₁ = kepler(μ, r₀exp, v₀exp, Δt)
        @test r₁ ≈ r₁exp rtol=1e-4
        @test v₁ ≈ v₁exp rtol=1e-5
        r₀, v₀ = kepler(μ, r₁, v₁, -Δt)
        @test r₀ ≈ r₀exp
        @test v₀ ≈ v₀exp
        r₁, v₁ = kepler(μ, r₀exp, v₀exp, eps())
        @test r₁ == r₀exp
        @test v₁ == v₀exp
        @test_throws ErrorException kepler(μ, r₀exp, v₀exp, Δt, 1)
    end
end

