@testset "Two Body" begin
    @testset "Keplerian <-> Cartesian" begin
        #= μ = 3.9860047e14 =#
        #=  =#
        #= a = 24464560.0 =#
        #= e = 0.7311 =#
        #= i = 0.122138 =#
        #= Ω = 1.00681 =#
        #= ω = 3.10686 =#
        #= v = 0.048363 =#
        #=  =#
        #= pos = [-0.107622532467967e+07, =#
        #=      -0.676589636432773e+07, =#
        #=      -0.332308783350379e+06, =#
        #=     ] =#
        #=  =#
        #= vel = [0.935685775154103e+04, =#
        #=      -0.331234775037644e+04, =#
        #=      -0.118801577532701e+04, =#
        #=     ] =#
        #=  =#
        #= pos₁, vel₁ = cartesian(a, e, i, Ω, ω, v, μ) =#
        #= a₁, e₁, i₁, Ω₁, ω₁, v₁ = keplerian(pos, vel, μ) =#

        #= @testset for i = 1:3 =#
        #=     @test pos₁[i] .≈ pos[i] =#
        #=     @test vel₁[i] .≈ vel[i] =#
        #= end =#
        #=  =#
        #= @test a ≈ a₁ =#
        #= @test e ≈ e₁ =#
        #= @test i ≈ i₁ =#
        #= @test Ω ≈ Ω₁ =#
        #= @test ω ≈ ω₁ =#
        #= @test v ≈ v₁ =#
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
end
