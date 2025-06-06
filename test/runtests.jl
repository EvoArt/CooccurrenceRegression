using CooccurrenceRegression, Random
using Test
X = 1 .- rand(LKJ(50,1))
Y = rand(Bool,3,50)
@testset "CooccurrenceRegression.jl" begin
    Random.seed!(123)
    chn1 =cooccurrence_regression(X,Y)
    Random.seed!(123)
    chn2 =cooccurrence_regression([X],Y)
    @test chn1[:α] == chn2[:α]
end
