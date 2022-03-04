@testset verbose = true "History" begin
    history = History([0.2, 0.8, 1.1], ["a", "b", "c"], 0.0, 2.0);
    @test duration(history) == 2.
    @test nb_events(history) == 3
    @test nb_events(history, 1.0, 2.0) == 1
    @test has_events(history)
    @test !has_events(history, 1.5, 2.0)
    push!(history, 1.7, "d")
    @test has_events(history, 1.5, 2.0)
end
