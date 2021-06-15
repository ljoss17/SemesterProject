include("sieve.jl")

using Test

@testset "binomial_k" begin
    @test abs(binomial_k(5, 0.5, 2)-0.3125) < 1e-16
end

@testset "sum_binomial" begin
    @test abs(sum_binomial(0, 4, 10, 0.2)-0.9672065024)/sum_binomial(0, 4, 10, 0.2) < 1e-15
    @test abs(sum_binomial(4, 10, 10, 0.2)-0.1208738816)/sum_binomial(4, 10, 10, 0.2) < 1e-15
    @test abs(sum_binomial(0, 6, 20, 0.2)-0.91330748643)/sum_binomial(0, 6, 20, 0.2) < 1e-10
    @test abs(sum_binomial(6, 20, 20, 0.2)-0.19579221454)/sum_binomial(6, 20, 20, 0.2) < 1e-10
    @test abs(sum_binomial(0, 15, 30, 0.2)-0.99994761271)/sum_binomial(0, 15, 30, 0.2) < 1e-10
    @test abs(sum_binomial(15, 30, 30, 0.2)-0.00023122561)/sum_binomial(15, 30, 30, 0.2) < 1e-8
    @test abs(sum_binomial(0, 30, 100, 0.2)-0.99394066452)/sum_binomial(0, 30, 100, 0.2) < 1e-8
    @test abs(sum_binomial(30, 100, 100, 0.2)-0.01124897872)/sum_binomial(30, 100, 100, 0.2) < 1e-8
end
