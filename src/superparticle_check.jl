include("model.jl")
include("params.jl")

FP = Dict{String,Float64}()
open(resultpath("fit_parameters.csv")) do io
    readline(io)
    for ii in eachline(io)
        at = findlast(',', ii); at === nothing && continue
        v = tryparse(Float64, strip(ii[at+1:end]))
        v === nothing || (FP[strip(ii[1:at-1])] = v)
    end
end
const Q_PRIM      = FP["Q_prim (nM um2/min) FIXED"]
const Q_CELL_10K  = FP["Q_cell (nM um2/min) ELISA-CONSTRAINED"]
const K_CELL      = FP["k_cell (1/min)"]
const KD_AUTO     = FP["Kd_auto (nM)"]
const BIAS_PRIM   = FP["bias_prim (tumor chemokine)"]
const BIAS_AUTO   = FP["bias_auto (orthokine)"]
const ROT_DIFF    = FP["rot_diffusion (rad²/min)"]
const D_FIXED     = FP["D (um2/min) FIXED"]
const K_FIXED     = FP["k tumor (1/min) FIXED"]
const KD_TUMOR    = FP["Kd tumor (nM) FIXED"]

const WELL_R  = well_radius_um("24-well")
const TUMOR_R = 1500.0
const RUN_T   = 72 * 60.0
const GRID_UM = 50.0
const GRID_R  = TUMOR_R + 1100.0
const STRAT   = :receptor
const PX_UM   = 14000.0 / 1188.0
const CELL_R  = CELL_DIAMETER_UM / 2
const CENTRE  = 250
const NPOS    = 2 * CENTRE + 1
const PX_AREA_MM2 = (PX_UM / 1000)^2
px_to_mm(p) = p * PX_UM / 1000

const N_REAL = 200_000
scale_for(n) = N_REAL / n
q_cell_for(n) = Q_CELL_10K * 10_000 / n
const ARMS = [(n = 10_000,  seeds = [11, 23, 37, 49, 61, 73, 89, 101],
               lab = "10,000 particles × 20", col = BP_GOLD,   lw = 5.0),
              (n = 200_000, seeds = [11, 23, 37],
               lab = "200,000 particles × 1", col = BP_INDIGO, lw = 1.8)]
const CONDS = [(key = "Par",  lab = "Off-target (parental)", feedback = false),
               (key = "CD19", lab = "Target (CD19)",         feedback = true)]
const HOURS = [24, 48, 72]

x_arr, dx, nx, mask = build_well_grid(GRID_R, GRID_UM)
center = reshape([0.0, 0.0], 1, 2)

println("forward check of the super-particle weighting, fitted parameters held fixed")
println("  Q_prim  = $(round(Q_PRIM; sigdigits=6)) nM·µm²/min")
for ii in ARMS
    println("  N = $(lpad(ii.n, 7))  SCALE = $(rpad(round(scale_for(ii.n); digits=2), 6)) ",
            "Q_cell = $(round(q_cell_for(ii.n); sigdigits=6)) nM·µm²/min  ",
            "→ total = $(round(q_cell_for(ii.n) * ii.n; sigdigits=6)) per fully-recruited well")
end
println("  mesh $(nx)×$(nx) at $(GRID_UM) µm, well r = $(Int(WELL_R)) µm, $(Int(RUN_T/60)) h")

const CELL_SAMPLES = let n = 16, o = Tuple{Float64,Float64}[]
    for ii in 1:n, iii in 1:n
        ox = -CELL_R + (ii - 0.5) * 2CELL_R / n
        oy = -CELL_R + (iii - 0.5) * 2CELL_R / n
        ox^2 + oy^2 <= CELL_R^2 && push!(o, (ox, oy))
    end
    o
end
const CELL_W = 1.0 / length(CELL_SAMPLES)

function render!(img, pos)
    fill!(img, 0.0)
    c = CENTRE + 1
    @inbounds for ii in 1:size(pos, 1)
        px = pos[ii, 1]; py = pos[ii, 2]
        for (ox, oy) in CELL_SAMPLES
            ix = round(Int, (px + ox) / PX_UM) + c
            iy = round(Int, (py + oy) / PX_UM) + c
            (1 <= ix <= NPOS && 1 <= iy <= NPOS) || continue
            img[ix, iy] += CELL_W
        end
    end
    img
end

@inline function bilerp(img, fx, fy)
    i0 = floor(Int, fx); j0 = floor(Int, fy)
    (1 <= i0 && i0 + 1 <= NPOS && 1 <= j0 && j0 + 1 <= NPOS) || return 0.0
    tx = fx - i0; ty = fy - j0
    @inbounds ((1 - tx) * (1 - ty) * img[i0, j0] + tx * (1 - ty) * img[i0+1, j0] +
               (1 - tx) * ty * img[i0, j0+1] + tx * ty * img[i0+1, j0+1])
end

fold(p) = [(p[CENTRE + 1 + ii] + p[CENTRE + 1 - ii]) / 2 for ii in 0:CENTRE]

function line_profile(img)
    c = CENTRE + 1
    acc = zeros(NPOS)
    for ii in 0:179
        θ = ii * pi / 180
        ux = cos(θ); uy = sin(θ)
        @inbounds for iii in 0:(NPOS - 1)
            t = iii - CENTRE
            acc[iii + 1] += bilerp(img, c + t * ux, c + t * uy)
        end
    end
    acc ./= 180
    fold(acc)
end

function run_one(ncells, Qc, seed)
    run_multisource_recruitment(;
        well_radius=WELL_R, island_centers=center, island_radius=TUMOR_R,
        island_Q=[Qc], n_cells=ncells, D=D_FIXED, k=K_FIXED, Kd=KD_TUMOR,
        ck_sensitivity=10.0, cell_speed=CELL_SPEED, rot_diffusion=ROT_DIFF,
        sensing_radius=SENSING_RADIUS, bias_weight=BIAS_PRIM, bias_auto=BIAS_AUTO,
        k_auto=K_CELL, Kd_auto=KD_AUTO, dt=DT, total_time=RUN_T,
        n_snapshots=length(HOURS), mask=mask, x_arr=x_arr, dx=dx, nx=nx,
        strategy=STRAT, n_seed_per_island=0, island_primary_Q=[Q_PRIM], seed=seed)
end

snap_at(snaps, hr) = argmin(abs(ii.t - hr * 60.0) for ii in snaps)

PROF = Dict{Tuple{Int,String,Int},Vector{Float64}}()
FRAC = Dict{Tuple{Int,String,Int},Float64}()
PEAK = Dict{Tuple{Int,String,Int},Float64}()

for ii in ARMS
    jobs = [(iii, iiii) for iii in CONDS for iiii in ii.seeds]
    prof = Vector{Dict{Int,Vector{Float64}}}(undef, length(jobs))
    frac = Vector{Dict{Int,Float64}}(undef, length(jobs))
    peak = Vector{Dict{Int,Float64}}(undef, length(jobs))
    t0 = time()
    Threads.@threads for iii in eachindex(jobs)
        c, sd = jobs[iii]
        img = zeros(NPOS, NPOS)
        snaps = run_one(ii.n, c.feedback ? q_cell_for(ii.n) : 0.0, sd)
        p = Dict{Int,Vector{Float64}}(); f = Dict{Int,Float64}(); k = Dict{Int,Float64}()
        for iiii in HOURS
            s = snaps[snap_at(snaps, iiii)]

            p[iiii] = line_profile(render!(img, s.pos)) .* scale_for(ii.n) ./ PX_AREA_MM2
            f[iiii] = s.n_rec[1] / ii.n
            k[iiii] = maximum(s.C_auto)
        end
        prof[iii] = p; frac[iii] = f; peak[iii] = k
    end
    for iii in CONDS
        idx = [iiii for iiii in eachindex(jobs) if jobs[iiii][1].key == iii.key]
        for iiii in HOURS
            PROF[(ii.n, iii.key, iiii)] = sum(prof[iiiii][iiii] for iiiii in idx) ./ length(idx)
            FRAC[(ii.n, iii.key, iiii)] = mean(frac[iiiii][iiii] for iiiii in idx)
            PEAK[(ii.n, iii.key, iiii)] = mean(peak[iiiii][iiii] for iiiii in idx)
        end
    end
    println("  N = $(ii.n): $(length(jobs)) runs in ",
            round((time() - t0) / 60; digits=1), " min")
end

const SKIP = 10
println("\nrecruited fraction at 72 h, and agreement over the fitted radii:")
for ii in CONDS
    for iii in HOURS
        a = PROF[(10_000,  ii.key, iii)][(SKIP + 1):end]
        b = PROF[(200_000, ii.key, iii)][(SKIP + 1):end]
        denom = maximum(b)
        rms = sqrt(mean((a .- b) .^ 2)) / denom * 100
        mx  = maximum(abs.(a .- b)) / denom * 100
        println("  $(rpad(ii.lab, 24)) $(iii) h  ",
                "frac $(round(FRAC[(10_000, ii.key, iii)]; digits=4)) vs ",
                "$(round(FRAC[(200_000, ii.key, iii)]; digits=4))   ",
                "peak orthokine $(round(PEAK[(10_000, ii.key, iii)]; sigdigits=3)) vs ",
                "$(round(PEAK[(200_000, ii.key, iii)]; sigdigits=3)) nM   ",
                "profile RMS $(round(rms; digits=2))%  max $(round(mx; digits=2))% of peak")
    end
end

open(resultpath("superparticle_check.csv"), "w") do io
    println(io, "n_particles,scale,condition,hours,mm_from_centre,density_real_cells_per_mm2")
    for ii in ARMS, iii in CONDS, iiii in HOURS
        v = PROF[(ii.n, iii.key, iiii)]
        for iiiii in 0:CENTRE
            println(io, ii.n, ",", scale_for(ii.n), ",", iii.key, ",", iiii, ",",
                    px_to_mm(iiiii), ",", v[iiiii + 1])
        end
    end
end
println("SAVED superparticle_check.csv")

ytop = ceil(maximum(maximum(PROF[ii]) for ii in keys(PROF)) / 100) * 100
panels = []
for ii in CONDS, iii in HOURS
    p = plot(; title="$(ii.lab)  ·  $(iii) h", titlefontsize=10, legend=false,
        xlabel=(ii.key == CONDS[end].key ? "Distance from tumor centre (mm)" : ""),
        ylabel=(iii == HOURS[1] ? "T cell density (cells / mm²)" : ""),
        xlim=(0, px_to_mm(CENTRE)), ylim=(0, ytop), size=(400, 300),
        tickfontsize=8, guidefontsize=9, grid=false)
    vspan!(p, [0, px_to_mm(SKIP)]; color=RGB(0.85, 0.85, 0.85), alpha=0.9, label="")
    for iiii in ARMS
        plot!(p, px_to_mm.(0:CENTRE), PROF[(iiii.n, ii.key, iii)];
              color=iiii.col, lw=iiii.lw, label="")
    end
    vline!(p, [TUMOR_R / 1000]; color=:black, lw=1.2, ls=:dot, label="")
    a = PROF[(10_000, ii.key, iii)][(SKIP + 1):end]
    b = PROF[(200_000, ii.key, iii)][(SKIP + 1):end]
    annotate!(p, px_to_mm(CENTRE) * 0.97, ytop * 0.93,
        text("RMS " * string(round(sqrt(mean((a .- b).^2)) / maximum(b) * 100; digits=2)) *
             "% of peak", 8, :right, :black))
    annotate!(p, px_to_mm(CENTRE) * 0.97, ytop * 0.84,
        text("recruited " * string(round(FRAC[(10_000, ii.key, iii)] * 100; digits=1)) * "% vs " *
             string(round(FRAC[(200_000, ii.key, iii)] * 100; digits=1)) * "%", 8, :right, :black))
    push!(panels, p)
end

leg = legend_strip([(ARMS[1].col, ARMS[1].col, ARMS[1].lab * " (as fitted)"),
                    (ARMS[2].col, ARMS[2].col, ARMS[2].lab * " (one per real cell)"),
                    (:white, :black, "tumor edge, 1.5 mm")];
                   width=1220, height=42, msize=8, fontsize=9)
fig = plot(leg, panels...; layout=Plots.@layout([a{0.07h}; Plots.grid(2, 3)]),
    size=(1220, 720),
    plot_title="Super-particle check · same fitted parameters, same total orthokine · " *
               "curves are REAL cells/mm²",
    plot_titlefontsize=12, bottom_margin=4Plots.mm, left_margin=8Plots.mm,
    right_margin=4Plots.mm, top_margin=2Plots.mm)
savefig_both(fig, "superparticle_check")
println("SAVED superparticle_check.png")
