@testset "Bodies" begin
    @testset for body in [PLANETS; SATELLITES]
        α = Dict(
            "Mercury" => 4.904544502878169,
            "Venus" => 4.760560067739733,
            "Earth" => 1.5314935667739428e-7,
            "Mars" => 5.544586951293573,
            "Jupiter" => 4.678480799964803,
            "Saturn" => 0.7084116986931903,
            "Uranus" => 4.4909241515991285,
            "Neptune" => 5.224359154533481,
            "Moon" => 4.6575372394964125,
        )
        δα = Dict(
            "Mercury" => -1.8140416085321446e-13,
            "Venus" => 0.0,
            "Earth" => -3.545123997161905e-12,
            "Mars" => -5.867982154428676e-13,
            "Jupiter" => -1.325795519124575e-13,
            "Saturn" => -1.9910212776572318e-13,
            "Uranus" => 0.0,
            "Neptune" => 3.5324603829189327e-12,
            "Moon" => -3.9302170086288704e-11,
        )
        δ = Dict(
            "Mercury" => 1.071881743978274,
            "Venus" => 1.1721631256393916,
            "Earth" => 1.5707964598747588,
            "Mars" => 0.9230435694063646,
            "Jupiter" => 1.1256642372977634,
            "Saturn" => 1.4579956981941933,
            "Uranus" => 6.018331593189447,
            "Neptune" => 0.749625184019545,
            "Moon" => 1.1457089606668125,
        )
        δδ = Dict(
            "Mercury" => -2.710001183477899e-14,
            "Venus" => 0.0,
            "Earth" => -3.0805523657085508e-12,
            "Mars" => -3.3681443280368176e-13,
            "Jupiter" => 2.3698823822682944e-14,
            "Saturn" => -2.212245864063591e-14,
            "Uranus" => 0.0,
            "Neptune" => 9.665269067238681e-14,
            "Moon" => 1.3081139255787021e-9,
        )
        w = Dict(
            "Mercury" => 5.698135863118986,
            "Venus" => 2.8089448431930744,
            "Earth" => 0.16849737156984945,
            "Mars" => 0.020664853154317875,
            "Jupiter" => 3.6596412821667506,
            "Saturn" => 6.169792117398191,
            "Uranus" => 1.6474170907902204,
            "Neptune" => 6.022111783780965,
            "Moon" => 0.6040059356860075,
        )
        δw = Dict(
            "Mercury" => 1.2398442877841507e-6,
            "Venus" => -2.9924494208700665e-7,
            "Earth" => 7.292115373194001e-5,
            "Mars" => 7.088218066303858e-5,
            "Jupiter" => 0.00017585323445765458,
            "Saturn" => 0.0001637849901848791,
            "Uranus" => -0.00010123719558981861,
            "Neptune" => 0.0001083382503473229,
            "Moon" => 2.6619514809869067e-6,
        )
        ep = TDBEpoch(2000, 1, 1)
        b = string(body)
        @eval begin
            @test right_ascension($body, $ep) ≈ $α[$b]
            @test right_ascension_rate($body, $ep) ≈ $δα[$b]/s
            @test declination($body, $ep) ≈ $δ[$b]
            @test declination_rate($body, $ep) ≈ $δδ[$b]/s
            @test rotation_angle($body, $ep) ≈ $w[$b]
            @test rotation_rate($body, $ep) ≈ $δw[$b]/s
        end
    end
end
