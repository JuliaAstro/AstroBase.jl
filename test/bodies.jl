# Reference values from SPICE
@testset "Bodies" begin
    @testset for body in [PLANETS; SATELLITES]
        ref = Dict(
            Mercury => (mod2pi.([0.1920568791315425, 0.49892931917252165, 1.4040611034965877]),
                [-1.8140416085263552e-13, 2.7100011833133014e-14, 1.240120048823325e-6]),
            Venus => (mod2pi.([0.04817108735504291, 0.398633201155505, -2.833859444572255]),
                [0.0, 0.0, -2.992449420870066e-7]),
            Earth => (mod2pi.([1.5688687239990866, 0.001674999621307277, -0.187235143671636]),
                [-3.5451239819165454e-12, 3.0805523657133963e-12, 7.292115373193997e-5]),
            Mars => (mod2pi.([0.8318788836923737, 0.6479359092404205, 3.0426176782275434]),
                [-5.867982099364807e-13, 3.3681443313368714e-13, 7.088218066303858e-5]),
            Jupiter => (mod2pi.([-0.03390662297737303, 0.44510949929258364, -1.5234175859907497]),
                [6.457809574650152e-14, -1.632996546363025e-14, 0.0001758532344576546]),
            Saturn => (mod2pi.([2.279099758383522, 0.11281265827898987, -1.7982494449054445]),
                # SPICE values for the Euler derivatives for Saturn's body-fixed frame are too
                # imprecise due to being extracted from the rotation matrix
                [-1.9910212776572318e-13, 2.2122458640635911e-14, 0.00016378499018487907]),
            Uranus => (mod2pi.([-0.2214648287855613, 1.835650040785036, -1.797580054095306]),
                [0.0, 0.0, -0.0001012371955898186]),
            Neptune => (mod2pi.([0.5138887918755141, 0.821113615356254, 0.41331795993755804]),
                [3.5096143643114534e-12, -3.0780073978414767e-13, 0.00010833825036298873]),
            Moon => (mod2pi.([-0.029263330369691223, 0.4330241293442422, 2.818632893330232]),
                [5.972450880751731e-11, -1.6820380693583266e-9, 2.6616471668504273e-6]),
        )

        ep = TDBEpoch(UTCEpoch(2017, 3, 25, 17, 04, 23.789))
        ref_angles = ref[body][1]
        ref_derivatives = ref[body][2]
        angles = mod2pi.(collect(euler_angles(body, ep)))
        derivatives = collect(euler_derivatives(body, ep))
        @testset "$body angle $i" for i = 1:3
            @test angles[i] ≈ ref_angles[i]
        end
        @testset "$body derivative $i" for i = 1:3
            @test derivatives[i] ≈ ref_derivatives[i]
        end
    end
end
