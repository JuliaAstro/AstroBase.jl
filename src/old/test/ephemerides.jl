struct DefunctEphemeris <: Ephemeris end

@testset "Ephemerides" begin
    load_ephemeris!(DefunctEphemeris())
    ep = TTEpoch(2000, 1, 1)
    @test_throws ErrorException state(ep, Earth)
    @test_throws ErrorException position(ep, Earth)
    @test_throws ErrorException velocity(ep, Earth)
end
