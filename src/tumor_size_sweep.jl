include("model.jl")
include("params.jl")

FP = Dict{String,Float64}()
open(resultpath("fit_parameters.csv")) do io
    readline(io)
    for ii in eachline(io)
        at = findlast(',', ii)
        at === nothing && continue
        v = tryparse(Float64, strip(ii[at+1:end]))
        v === nothing || (FP[strip(ii[1:at-1])] = v)
    end
end
const Q_PRIM    = FP["Q_prim (nM um2/min) FIXED"]
const Q_CELL    = FP["Q_cell (nM um2/min) ELISA-CONSTRAINED"]
const K_CELL    = FP["k_cell (1/min)"]
const KD_AUTO   = FP["Kd_auto (nM)"]
const BIAS_PRIM = FP["bias_prim (tumor chemokine)"]
const BIAS_AUTO = FP["bias_auto (orthokine)"]
const ROT_DIFF  = FP["rot_diffusion (rad²/min)"]
const D_FIXED   = FP["D (um2/min) FIXED"]
const K_FIXED   = FP["k tumor (1/min) FIXED"]
const KD_TUMOR  = FP["Kd tumor (nM) FIXED"]
const LAMBDA    = sqrt(D_FIXED / K_FIXED)
const LAMBDA_C  = sqrt(D_FIXED / K_CELL)

const WELL_R  = well_radius_um("24-well")
const RUN_T   = 72 * 60.0
const STRAT   = :receptor

const N_CELLS = parse(Int, get(ENV, "FIT_NCELLS", "10000"))
const GRID_UM = 50.0
const PX_UM   = 14000.0 / 1188.0
const CELL_R  = CELL_DIAMETER_UM / 2

const RADII = [750.0, 1500.0, 2250.0, 3500.0]
const SEEDS = [11, 23, 37, 49, 61, 73, 89, 101,
               113, 127, 139, 151, 163, 177, 191, 203]

const SCEN = [(key = "tumor",  Qc = 0.0,     lab = "Tumor only",
               col = BP_GRAY,   short = "tumor only"),
              (key = "feedback", Qc = Q_CELL, lab = "Tumor + T cell positive feedback",
               col = BP_INDIGO, short = "+ T cell feedback")]

const HALF_MM = 5.5
const CENTRE  = round(Int, HALF_MM * 1000 / PX_UM)
const NPOS    = 2 * CENTRE + 1
const PX_AREA_MM2 = (PX_UM / 1000)^2
px_to_mm(p) = p * PX_UM / 1000

println("ELISA-constrained parameters in use:")
println("  Q_prim = $(round(Q_PRIM; sigdigits=6))   Q_cell = $(round(Q_CELL; sigdigits=6))")
println("  bias_prim = $(round(BIAS_PRIM; digits=4))   ",
        "bias_auto = $(round(BIAS_AUTO; digits=4))")
println("  rot_diffusion = $(round(ROT_DIFF; sigdigits=6)) rad²/min")
println("  CXCL12:    k = $(round(K_FIXED; sigdigits=4)) /min, ",
        "λ = $(round(LAMBDA; digits=1)) µm, Kd = $(KD_TUMOR) nM")
println("  orthokine: k_cell = $(round(K_CELL; sigdigits=4)) /min, ",
        "λ_cell = $(round(LAMBDA_C; digits=1)) µm, Kd_auto = $(round(KD_AUTO; sigdigits=4)) nM")
println("tumor radii swept: $(Int.(RADII)) µm; $(length(SEEDS)) seeds; $(N_CELLS) cells")

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

function run_one(R, Qc, seed; nframes)
    x_arr, dx, nx, mask = build_well_grid(R + 1100.0, GRID_UM)
    run_multisource_recruitment(;
        well_radius=WELL_R, island_centers=reshape([0.0, 0.0], 1, 2),
        island_radius=R, island_Q=[Qc], n_cells=N_CELLS,
        D=D_FIXED, k=K_FIXED, Kd=KD_TUMOR, ck_sensitivity=10.0,
        cell_speed=CELL_SPEED, rot_diffusion=ROT_DIFF,
        sensing_radius=SENSING_RADIUS, bias_weight=BIAS_PRIM, bias_auto=BIAS_AUTO,
        k_auto=K_CELL, Kd_auto=KD_AUTO, dt=DT, total_time=RUN_T, n_snapshots=nframes,
        mask=mask, x_arr=x_arr, dx=dx, nx=nx, strategy=STRAT,
        n_seed_per_island=0, island_primary_Q=[Q_PRIM], seed=seed)
end

PROF = Dict{Tuple{Float64,String},Vector{Float64}}()
NREC = Dict{Tuple{Float64,String},Float64}()
for ii in RADII, iii in SCEN

    profs = Vector{Vector{Float64}}(undef, length(SEEDS))
    nrs   = Vector{Float64}(undef, length(SEEDS))
    Threads.@threads for iiii in eachindex(SEEDS)
        img = zeros(NPOS, NPOS)
        f = run_one(ii, iii.Qc, SEEDS[iiii]; nframes=1)[end]
        profs[iiii] = line_profile(render!(img, f.pos))
        nrs[iiii]   = f.n_rec[1]
    end
    acc = sum(profs); nr = nrs
    acc ./= length(SEEDS)
    PROF[(ii, iii.key)] = acc ./ PX_AREA_MM2
    NREC[(ii, iii.key)] = mean(nr)
    println("  R = $(Int(ii)) µm  $(rpad(iii.lab, 36)) ",
            "recruited $(round(NREC[(ii, iii.key)]; digits=1)) / $(N_CELLS)",
            "  peak density $(round(maximum(PROF[(ii, iii.key)][11:end]); digits=1)) cells/mm²")
end

open(resultpath("tumor_size_profiles.csv"), "w") do io
    println(io, "tumor_radius_um,scenario,mm_from_centre,density_cells_per_mm2")
    for ii in RADII, iii in SCEN
        p = PROF[(ii, iii.key)]
        for iiii in 0:CENTRE
            println(io, Int(ii), ",", iii.key, ",", round(px_to_mm(iiii); digits=4), ",",
                    round(p[iiii + 1]; digits=4))
        end
    end
end
println("SAVED tumor_size_profiles.csv")

const SKIP = 10
RMM  = [px_to_mm(ii) for ii in SKIP:CENTRE]
ytop = ceil(maximum(maximum(PROF[(ii, iii.key)][SKIP+1:end])
                    for ii in RADII for iii in SCEN) / 50) * 50

panels = []
for (ii, R) in enumerate(RADII)
    p = plot(; title="Tumor radius $(Int(R)) µm" * (R == 1500 ? "  (as fitted)" : ""),
        titlefontsize=11, legend=false,
        xlabel="Distance from tumor centre (mm)",
        ylabel=(ii == 1 ? "T cell density (cells / mm²)" : ""),
        xlim=(0, HALF_MM), ylim=(0, ytop), size=(430, 380),
        tickfontsize=9, guidefontsize=10, grid=false)
    for iii in SCEN
        plot!(p, RMM, PROF[(R, iii.key)][SKIP+1:end]; color=iii.col, lw=2.2, label="")
    end
    vline!(p, [R / 1000]; color=:black, lw=1.2, ls=:dot, label="")
    push!(panels, p)
end

const PROF_W = 420 * length(RADII)
leg = legend_strip([(SCEN[1].col, SCEN[1].col, SCEN[1].lab),
                    (SCEN[2].col, SCEN[2].col, SCEN[2].lab),
                    (:white, :black, "tumor edge")];
                   width=PROF_W, height=42, msize=8, fontsize=10)
fig = plot(leg, panels...;
    layout=Plots.@layout([a{0.09h}; Plots.grid(1, length(RADII))]),
    size=(PROF_W, 460),
    plot_title="T cell distribution at 72 h vs tumor size · " *
               "ELISA-constrained parameters · ΔFRO",
    plot_titlefontsize=12, bottom_margin=5Plots.mm, left_margin=6Plots.mm,
    top_margin=2Plots.mm)
savefig_both(fig, "tumor_size_profiles")
println("SAVED tumor_size_profiles.png")
