using Test
using SBM_Bioreactor

# Fast default suite with plotting and slow scientific checks disabled.
ENV["SBM_RUN_SLOW_TESTS"] = "0"
ENV["SBM_RUN_PLOTS_TESTS"] = "0"

include(joinpath(@__DIR__, "..", "test", "runtests.jl"))
