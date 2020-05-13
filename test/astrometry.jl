@testset "Astrometry" begin
    @testset "ab" begin
        pnat = [-0.76321968546737951,-0.60869453983060384,-0.21676408580639883]
        v = [2.1044018893653786e-5,-8.9108923304429319e-5,-3.8633714797716569e-5]
        s = 0.99980921395708788
        bm1 = 0.99999999506209258
        exp = ERFA.ab(pnat, v, s, bm1)
        act = aberration(pnat, v, s, bm1)
        @test act[1] ≈ exp[1] atol=1e-12
        @test act[2] ≈ exp[2] atol=1e-12
        @test act[3] ≈ exp[3] atol=1e-12
    end
end
