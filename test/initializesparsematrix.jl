using SparseArrays
using scran
                                               
@testset "Initialization tests, 32-bit indices" begin
    function spawner(type::Type, m::Integer, n::Integer, density::AbstractFloat) 
        raw = abs.(sprand(type, m, n, density))
        return SparseMatrixCSC(m, n, Vector{Int32}(raw.colptr), Vector{Int32}(raw.rowval), raw.nzval)
    end

    @testset "Int32 values, Int32 indices" begin
        x = spawner(Int32, 8, 12, 0.2)
        mat = scran.initializesparsematrix(x)
        @test mat isa scran.ScranMatrix
        @test size(mat) == (8, 12)
        @test size(mat, 1) == 8
        @test size(mat, 2) == 12
        @test isapprox(scran.extractcolumn(mat, 1), Vector{Float64}(x[:,1]))
        @test isapprox(scran.extractrow(mat, 2), Vector{Float64}(x[2,:]))
    end

    @testset "Int64 values, Int32 indices" begin
        x = spawner(Int64, 50, 11, 0.2)
        mat = scran.initializesparsematrix(x; forceint = false) # otherwise, it's forced to a 32-bit int.
        @test mat isa scran.ScranMatrix
        @test size(mat) == (50, 11)
        @test size(mat, 1) == 50
        @test size(mat, 2) == 11
        @test isapprox(scran.extractcolumn(mat, 10), Vector{Float64}(x[:,10]))
        @test isapprox(scran.extractrow(mat, 24), Vector{Float64}(x[24,:]))
    end

    @testset "Float32 values, Int32 indices" begin
        x = spawner(Float32, 20, 20, 0.2)
        mat = scran.initializesparsematrix(x) 
        @test mat isa scran.ScranMatrix
        @test size(mat) == (20, 20)
        @test size(mat, 1) == 20
        @test size(mat, 2) == 20
        @test isapprox(scran.extractcolumn(mat, 9), trunc.(Vector{Float64}(x[:,9])))
        @test isapprox(scran.extractrow(mat, 14), trunc.(Vector{Float64}(x[14,:])))
    end

    @testset "Float64 values, Int32 indices" begin
        x = spawner(Float64, 25, 30, 0.2)
        mat = scran.initializesparsematrix(x) 
        @test mat isa scran.ScranMatrix
        @test size(mat) == (25, 30)
        @test size(mat, 1) == 25
        @test size(mat, 2) == 30
        @test isapprox(scran.extractcolumn(mat, 22), trunc.(Vector{Float64}(x[:,22])))
        @test isapprox(scran.extractrow(mat, 5), trunc.(Vector{Float64}(x[5,:])))
    end
end

@testset "Initialization tests, 64-bit indices" begin
    @testset "Int32 values, Int64 indices" begin
        x = abs.(sprand(Int32, 8, 12, 0.2))
        mat = scran.initializesparsematrix(x)
        @test mat isa scran.ScranMatrix
        @test size(mat) == (8, 12)
        @test size(mat, 1) == 8
        @test size(mat, 2) == 12
        @test isapprox(scran.extractcolumn(mat, 1), Vector{Float64}(x[:,1]))
        @test isapprox(scran.extractrow(mat, 2), Vector{Float64}(x[2,:]))
    end

    @testset "Int64 values, Int64 indices" begin
        x = abs.(sprand(Int64, 50, 11, 0.2)) 
        mat = scran.initializesparsematrix(x; forceint = false) # otherwise, it's forced to a 32-bit int.
        @test mat isa scran.ScranMatrix
        @test size(mat) == (50, 11)
        @test size(mat, 1) == 50
        @test size(mat, 2) == 11
        @test isapprox(scran.extractcolumn(mat, 10), Vector{Float64}(x[:,10]))
        @test isapprox(scran.extractrow(mat, 24), Vector{Float64}(x[24,:]))
    end

    @testset "Float32 values, Int64 indices" begin
        x = abs.(sprand(Float32, 20, 20, 0.2)) 
        mat = scran.initializesparsematrix(x) 
        @test mat isa scran.ScranMatrix
        @test size(mat) == (20, 20)
        @test size(mat, 1) == 20
        @test size(mat, 2) == 20
        @test isapprox(scran.extractcolumn(mat, 9), trunc.(Vector{Float64}(x[:,9])))
        @test isapprox(scran.extractrow(mat, 14), trunc.(Vector{Float64}(x[14,:])))
    end

    @testset "Float64 values, Int64 indices" begin
        x = abs.(sprand(Float64, 25, 30, 0.2)) 
        mat = scran.initializesparsematrix(x) 
        @test mat isa scran.ScranMatrix
        @test size(mat) == (25, 30)
        @test size(mat, 1) == 25
        @test size(mat, 2) == 30
        @test isapprox(scran.extractcolumn(mat, 22), trunc.(Vector{Float64}(x[:,22])))
        @test isapprox(scran.extractrow(mat, 5), trunc.(Vector{Float64}(x[5,:])))
    end
end

@testset "Initialization tests, more options" begin
    @testset "Floats with no forced integers" begin
        x = abs.(sprand(Float64, 25, 30, 0.2)) 
        mat = scran.initializesparsematrix(x; forceint = false) 
        @test isapprox(scran.extractcolumn(mat, 22), Vector{Float64}(x[:,22]))
        @test isapprox(scran.extractrow(mat, 5), Vector{Float64}(x[5,:]))
    end

    @testset "Avoiding a copy" begin
        x = abs.(sprand(Int32, 25, 30, 0.2)) 
        mat = scran.initializesparsematrix(x; nocopy = true) 
        @test isapprox(scran.extractcolumn(mat, 12), trunc.(Vector{Float64}(x[:,12])))
        @test isapprox(scran.extractrow(mat, 15), trunc.(Vector{Float64}(x[15,:])))
    end
end
