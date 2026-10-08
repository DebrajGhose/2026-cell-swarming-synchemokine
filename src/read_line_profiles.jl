include("model.jl")
include("params.jl")
include("xlsx_read.jl")

const LPDIR  = datapath("line_profiles")
const CENTRE = 250

const PX_UM = 14000.0 / 1188.0
mm_from_centre(p) = (p - CENTRE) * PX_UM / 1000

SOURCES = [("Par",  "LineProfiles_Par_circuit.xlsx",  "Off-target (parental)"),
           ("CD19", "LineProfiles_CD19_circuit.xlsx", "Target (CD19)")]

PROF = Dict{Tuple{String,Int,Int}, Vector{Float64}}()
POS  = Float64[]
for (key, file, _) in SOURCES
    want = ["Position", "Average"]
    sheets = read_xlsx(joinpath(LPDIR, file); want = want)
    println("$(file): $(length(sheets)) sheets")
    for iii in sheets
        rep = parse(Int, match(r"Replicate(\d+)", iii.name).captures[1])
        hr  = parse(Int, match(r"_(\d+)h$", iii.name).captures[1])
        a   = iii.cols["Average"]; p = iii.cols["Position"]
        isempty(POS) && (global POS = p)
        PROF[(key, rep, hr)] = a
        asym = sqrt(mean((a .- reverse(a)) .^ 2))
        println("  $(rpad(iii.name, 24)) n=$(length(a))  ",
                "centre=$(round(a[CENTRE+1]; digits=2))  ",
                "max=$(round(maximum(a); digits=2)) @$(Int(p[argmax(a)]))  ",
                "background=$(round((a[1]+a[end])/2; digits=2))  ",
                "asymmetry RMS=$(round(asym; digits=6))")
    end
end
REPS  = sort(unique(ii for (_, ii, _) in keys(PROF)))
HOURS = sort(unique(ii for (_, _, ii) in keys(PROF)))
meanprof(key, hr) = mean(PROF[(key, ii, hr)] for ii in REPS)

println("\n=== replicate mean per condition and timepoint ===")
for (key, _, lab) in SOURCES, iii in HOURS
    m = meanprof(key, iii)
    println("  $(rpad(lab, 22)) $(rpad(iii, 3)) h  peak = $(round(maximum(m); digits=2))  ",
            "background = $(round((m[1] + m[end]) / 2; digits=2))")
end

open(datapath("infiltration_measured.csv"), "w") do io
    println(io, "condition,replicate,hours,position_px,intensity")
    for ((key, rep, hr), a) in sort(collect(PROF); by = first), iii in eachindex(a)
        println(io, key, ",", rep, ",", hr, ",", Int(POS[iii]), ",",
                round(a[iii]; digits=4))
    end
end
open(datapath("infiltration_measured_mean.csv"), "w") do io
    println(io, "condition,hours,position_px,mean_intensity,mm_from_centre")
    for (key, _, _) in SOURCES, iii in HOURS
        m = meanprof(key, iii)
        for iiii in eachindex(m)
            println(io, key, ",", iii, ",", Int(POS[iiii]), ",",
                    round(m[iiii]; digits=4), ",",
                    round(mm_from_centre(POS[iiii]); digits=4))
        end
    end
end
println("SAVED infiltration_measured.csv, infiltration_measured_mean.csv")
