using PackageCompiler
using Libdl

const ROOT = normpath(joinpath(@__DIR__, ".."))
const WORKLOAD = joinpath(ROOT, "scripts", "precompile_sysimage_workload.jl")
const OUTPUT = joinpath(ROOT, "artifacts", "SBM_Bioreactor_sysimage.$(Libdl.dlext)")

mkpath(dirname(OUTPUT))

create_sysimage(
    ["SBM_Bioreactor"];
    project = ROOT,
    sysimage_path = OUTPUT,
    precompile_execution_file = WORKLOAD,
    incremental = true,
)

println(OUTPUT)
