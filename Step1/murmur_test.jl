include("murmur.jl")

using Test

@testset "compute_p" begin
    @test compute_p(50, 100) == 0.75
    @test compute_p(0, 100) == 0
    @test compute_p(1, 1e10) != 0
    @test_throws ArgumentError compute_p(100, 50)
    @test_throws ArgumentError compute_p(50, 0)
end

@testset "sum_body" begin
    @test sum_body(2, 0.5, 4) == 0.375
    @test abs(sum_body(1, 0.51, 9)-0.0299096375126409) < 1e-17
    @test abs(sum_body(1, 0.51, 9)-0.0299096375126409)/sum_body(1, 0.51, 9) < 1e-15
    @test abs(sum_body(2, 0.51, 9)-0.001655951531561063859516483) < 1e-18
    @test abs(sum_body(2, 0.51, 9)-0.001655951531561063859516483)/sum_body(2, 0.51, 9) < 1e-15
    @test abs(sum_body(3, 0.51, 9)-0.000222745391052210891629437964054484) < 1e-19
    @test abs(sum_body(3, 0.51, 9)-0.000222745391052210891629437964054484)/sum_body(3, 0.51, 9) < 1e-15
    @test abs(sum_body(4, 0.51, 9)-0.0000802217525874537526203420827542224126) < 1e-19
    @test abs(sum_body(4, 0.51, 9)-0.0000802217525874537526203420827542224126)/sum_body(4, 0.51, 9) < 1e-15
    @test abs(sum_body(5, 0.51, 9)-0.0000802217525874537526203420827542224126) < 1e-19
    @test abs(sum_body(5, 0.51, 9)-0.0000802217525874537526203420827542224126)/sum_body(5, 0.51, 9) < 1e-15
end

@testset "murmur_totality" begin
    #@test  == 1.6677181699666569
    @test abs(murmur_totality(3, 9, 10, 0.1)-0.0319487779404282722563866057295629288252) < 1e-15
    @test abs(murmur_totality(3, 9, 10, 0.1)-0.0319487779404282722563866057295629288252)/murmur_totality(3, 9, 10, 0.1) < 1e-13
    @test_throws ArgumentError murmur_totality(200, 10, 100, 1.1)
    @test_throws ArgumentError murmur_totality(200, 10, 100, -0.1)
end
