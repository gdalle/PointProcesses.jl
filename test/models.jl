@testset verbose = true "Models" begin

    @testset "Multivariate Poisson - built-in" begin
        pp = MultivariatePoissonProcess(rand(10))
        h = rand(pp, 0.0, 1000.0)
        pp_est = fit_mle(MultivariatePoissonProcess, h)
        error = mean(abs, pp_est.λ - pp.λ)
        @test error < 0.1
    end

    @testset "Multivariate Poisson - generic" begin
        pp = PoissonProcess(λ=1., mark_dist=Categorical(randprobvec(10)))
        h = rand(pp, 0.0, 1000.0)
        pp_est = fit_mle(PoissonProcess{Categorical}, h)
        λ_error = mean(abs, pp_est.λ - pp.λ)
        @test λ_error < 0.1
        p_error = mean(abs, pp_est.mark_dist.p - pp.mark_dist.p)
        @test p_error < 0.1
    end

    @testset "Multivariate Poisson - naive" begin
        struct MyPoissonProcess <: TemporalPointProcess{Int,Float64}
            λ::Vector{Float64}
        end

        function PointProcesses.mark_distribution(pp::MyPoissonProcess, t=nothing, h=nothing)
            return Categorical(pp.λ / sum(pp.λ))
        end

        function PointProcesses.intensity(pp::MyPoissonProcess, m, t=nothing, h=nothing)
            return pp.λ[m]
        end

        function PointProcesses.ground_intensity(pp::MyPoissonProcess, t=nothing, h=nothing)
            return sum(pp.λ)
        end

        function PointProcesses.ground_intensity_bound(pp::MyPoissonProcess, t=nothing, h=nothing)
            return ground_intensity(pp, t, h), Inf
        end

        function PointProcesses.integrated_ground_intensity(pp::MyPoissonProcess, h, a, b)
            par = [pp]
            f = (t, par) -> ground_intensity(par[1], h, t)
            prob = QuadratureProblem(f, a, b, par)
            sol = solve(prob, HCubatureJL())
            return sol.u
        end

        pp = MyPoissonProcess(rand(10))
        h = rand(pp, 0.0, 1000.0)
        pp_est = fit_mle(MultivariatePoissonProcess, h)
        error = mean(abs, pp_est.λ - pp.λ)
        @test error < 0.1
        logL = logdensityof(pp, h)
        logL_ref = logdensityof(MultivariatePoissonProcess(pp.λ), h)
        @test logL ≈ logL_ref
    end

end
