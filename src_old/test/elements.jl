using DelimitedFiles

@testset "Elements" begin
    # Reference values from Orekit
    μ = 3.986004418e5
    rexp = [1131.340, -2282.343, 6672.423]
    vexp = [-5.64305, 4.30333, 2.42879]
    sma_exp = 7200.470581180567
    ecc_exp = 0.008100116890743586
    inc_exp = 1.7208944567902595
    node_exp = 5.579892976386111
    peri_exp = 1.2370820968712155
    ano_exp = 7.194559370904103E-5
    sma, ecc, inc, node, peri, ano = keplerian(rexp, vexp, μ)
    @test sma ≈ sma_exp
    @test ecc ≈ ecc_exp
    @test inc ≈ inc_exp
    @test node ≈ node_exp
    @test peri ≈ peri_exp
    @test ano ≈ ano_exp
    r, v = cartesian(sma, ecc, inc, node, peri, ano, μ)
    @test r ≈ rexp
    @test v ≈ vexp

    @test !isprograde(inc)
    @test isretrograde(inc)
    @test !ispolar(inc)

    ref_rv = readdlm(joinpath("data", "rv.csv"), ',')
    ref_el = readdlm(joinpath("data", "elements.csv"), ',')

    n, _ = size(ref_rv)
    for i = 1:n
        r_exp = ref_rv[i, 1:3]
        v_exp = ref_rv[i, 4:6]
        el_exp = ref_el[i, :]
        r_act, v_act = cartesian(el_exp, μ)
        @test r_act ≈ r_exp
        @test v_act ≈ v_exp
        el_act = collect(keplerian(r_exp, v_exp, μ))
        @test el_act ≈ el_exp
    end
end
