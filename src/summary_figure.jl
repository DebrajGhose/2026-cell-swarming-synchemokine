include("model.jl")
include("params.jl")
import CairoMakie as CM


const N_CELLS = parse(Int, get(ENV, "FIT_NCELLS", "10000"))

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
const K_CELL    = FP["k_cell (1/min)"]
const KD_AUTO   = FP["Kd_auto (nM)"]
const BIAS_PRIM = FP["bias_prim (tumor chemokine)"]
const BIAS_AUTO = FP["bias_auto (orthokine)"]
const ROT_DIFF  = FP["rot_diffusion (rad²/min)"]
const Q_PRIM    = FP["Q_prim (nM um2/min) FIXED"]
const Q_CELL    = FP["Q_cell (nM um2/min) ELISA-CONSTRAINED"]
const D_FIXED   = FP["D (um2/min) FIXED"]
const K_FIXED   = FP["k tumor (1/min) FIXED"]
const KD_TUMOR  = FP["Kd tumor (nM) FIXED"]
const GAIN_A    = FP["intensity gain a"]
const BKG_B     = FP["intensity background b"]
const R2_STORED = FP["R2 all profiles"]

const WELL_R  = well_radius_um("24-well")
const TUMOR_R = 1500.0
const RUN_T   = 72 * 60.0
const STRAT   = :receptor
const GRID_UM = 50.0
const PX_UM   = 14000.0 / 1188.0
const CELL_R  = CELL_DIAMETER_UM / 2
px_to_mm(p) = p * PX_UM / 1000

const CENTRE      = 250
const NPOS        = 2 * CENTRE + 1
const CENTRE_SKIP = 10
const RBIN        = 5
const FIT_HOURS   = [24, 47, 72]
const N_SNAPS     = round(Int, RUN_T / 60.0)
const FINAL_SEEDS = [101, 113, 127, 139, 151, 163, 177, 191]
const CONDS       = ["Par", "CD19"]

const RADII    = [750.0, 1500.0, 2250.0, 3500.0]
const GIF_SEED = 11
const ZOOM_D   = 4300.0
const ZOOM_C   = 2400.0
const TS_HOURS = collect(0:8:72)
const HALF_MM  = 5.5
const SCEN = [(key = "tumor",    Qc = 0.0,    lab = "Tumor only",
               col = BP_GRAY,   short = "tumor only",
               rowlab = "Tumor only\n(no feedback)"),
              (key = "feedback", Qc = Q_CELL, lab = "Tumor + T cell positive feedback",
               col = BP_INDIGO, short = "+ T cell feedback",
               rowlab = "Tumor +\nT cell feedback")]

const COL_NAIVE = PC_TEAL
const COL_REC   = ISLAND_COLORS[1]
const DOT_ALPHA = 0.45
const BIN_UM    = 150.0

println("ELISA-constrained parameters in use ($(N_CELLS) simulated T cells):")
println("  k_cell = $(round(K_CELL; sigdigits=4)) /min   Kd_auto = $(round(KD_AUTO; sigdigits=4)) nM")
println("  bias_prim = $(round(BIAS_PRIM; digits=4))   bias_auto = $(round(BIAS_AUTO; digits=4))")
println("  rot_diffusion = $(round(ROT_DIFF; sigdigits=6)) rad²/min")
println("  Q_prim = $(round(Q_PRIM; sigdigits=6))   Q_cell = $(round(Q_CELL; sigdigits=6))")
println("  gain a = $(round(GAIN_A; sigdigits=5))   background b = $(round(BKG_B; digits=3))")

rawmean = Dict{Tuple{String,Int},Vector{Float64}}()
open(datapath("infiltration_measured_mean.csv")) do io
    readline(io)
    for ii in eachline(io)
        f = split(strip(ii), ',')
        hr = parse(Int, f[2]); hr in FIT_HOURS || continue
        key = (String(f[1]), hr)
        haskey(rawmean, key) || (rawmean[key] = Float64[])
        push!(rawmean[key], parse(Float64, f[4]))
    end
end
const PANELS = [(ii, iii) for ii in CONDS for iii in FIT_HOURS]
fold(p) = [(p[CENTRE + 1 + ii] + p[CENTRE + 1 - ii]) / 2 for ii in 0:CENTRE]
OBSHALF = Dict(ii => fold(rawmean[ii]) for ii in PANELS)
const RBINS = [(ii, min(ii + RBIN - 1, CENTRE)) for ii in (CENTRE_SKIP + 1):RBIN:CENTRE]
binmean(h, (a, b)) = mean(h[(a + 1):(b + 1)])
OBSBIN = Dict(ii => [binmean(OBSHALF[ii], iii) for iii in RBINS] for ii in PANELS)
const RCENT = [sum(ii) / 2 for ii in RBINS]

const CELL_SAMPLES = let n = 16, o = Tuple{Float64,Float64}[]
    for ii in 1:n, iii in 1:n
        ox = -CELL_R + (ii - 0.5) * 2CELL_R / n
        oy = -CELL_R + (iii - 0.5) * 2CELL_R / n
        ox^2 + oy^2 <= CELL_R^2 && push!(o, (ox, oy))
    end
    o
end
const CELL_W = 1.0 / length(CELL_SAMPLES)
const PX_AREA_MM2 = (PX_UM / 1000)^2

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

function line_at(img, deg)
    c = CENTRE + 1
    θ = deg * pi / 180; ux = cos(θ); uy = sin(θ)
    [bilerp(img, c + (ii - CENTRE) * ux, c + (ii - CENTRE) * uy) for ii in 0:(NPOS - 1)]
end

function line_mean(img)
    acc = zeros(NPOS)
    for ii in 0:179
        acc .+= line_at(img, ii)
    end
    acc ./ 180
end

line_profile(img) = fold(line_mean(img))

const CENTER = reshape([0.0, 0.0], 1, 2)

function run_one(R, Qc, seed; nframes, ncells=N_CELLS)
    x_arr, dx, nx, mask = build_well_grid(R + 1100.0, GRID_UM)
    run_multisource_recruitment(;
        well_radius=WELL_R, island_centers=CENTER, island_radius=R,
        island_Q=[Qc], n_cells=ncells, D=D_FIXED, k=K_FIXED, Kd=KD_TUMOR,
        ck_sensitivity=10.0, cell_speed=CELL_SPEED, rot_diffusion=ROT_DIFF,
        sensing_radius=SENSING_RADIUS, bias_weight=BIAS_PRIM, bias_auto=BIAS_AUTO,
        k_auto=K_CELL, Kd_auto=KD_AUTO, dt=DT, total_time=RUN_T,
        n_snapshots=nframes, mask=mask, x_arr=x_arr, dx=dx, nx=nx,
        strategy=STRAT, n_seed_per_island=0, island_primary_Q=[Q_PRIM], seed=seed)
end

snap_at(snaps, hr) = argmin(abs(ii.t - hr * 60.0) for ii in snaps)

t0 = time()
fitjobs = [(ii, iii) for ii in CONDS for iii in FINAL_SEEDS]
fitprof = Vector{Dict{Int,Vector{Float64}}}(undef, length(fitjobs))
fitrec  = Vector{Dict{Int,Float64}}(undef, length(fitjobs))
Threads.@threads for ii in eachindex(fitjobs)
    c, sd = fitjobs[ii]
    img = zeros(NPOS, NPOS)
    snaps = run_one(TUMOR_R, c == "CD19" ? Q_CELL : 0.0, sd; nframes=N_SNAPS)
    p = Dict{Int,Vector{Float64}}(); n = Dict{Int,Float64}()
    for iii in FIT_HOURS
        f = snaps[snap_at(snaps, iii)]
        p[iii] = line_profile(render!(img, f.pos))
        n[iii] = f.n_rec[1]
    end
    fitprof[ii] = p; fitrec[ii] = n
end
SIMBIN = Dict{Tuple{String,Int},Vector{Float64}}()
NRECF  = Dict{Tuple{String,Int},Float64}()
for ii in CONDS
    idx = [iii for iii in eachindex(fitjobs) if fitjobs[iii][1] == ii]
    for iii in FIT_HOURS
        acc = sum(fitprof[iiii][iii] for iiii in idx) ./ length(idx)
        SIMBIN[(ii, iii)] = [binmean(acc, iiii) for iiii in RBINS]
        NRECF[(ii, iii)]  = mean(fitrec[iiii][iii] for iiii in idx)
    end
end
println("panel B: $(length(fitjobs)) runs in $(round(time() - t0; digits=1)) s")

fitted(k) = GAIN_A .* SIMBIN[k] .+ BKG_B
r2_one(k) = (o = OBSBIN[k]; f = fitted(k);
             1 - sum((o .- f) .^ 2) / sum((o .- mean(o)) .^ 2))
allobs = vcat((OBSBIN[ii] for ii in PANELS)...)
allfit = vcat((fitted(ii) for ii in PANELS)...)
R2 = 1 - sum((allobs .- allfit) .^ 2) / sum((allobs .- mean(allobs)) .^ 2)
println("  R² over the six profiles = $(round(R2; digits=4))  ",
        "(stored $(round(R2_STORED; digits=4)); these must agree)")
abs(R2 - R2_STORED) > 5e-3 && @warn "R² does not reproduce the stored value; " *
    "panel B is not the same run as the fit figure."
for ii in PANELS
    println("    $(rpad(ii[1], 5)) $(ii[2])h  R² = $(round(r2_one(ii); digits=3))   ",
            "recruited $(round(NRECF[ii]; digits=1))")
end

t1 = time()
tsruns = Vector{Any}(undef, length(SCEN))
Threads.@threads for ii in eachindex(SCEN)
    tsruns[ii] = run_one(TUMOR_R, SCEN[ii].Qc, GIF_SEED; nframes=N_SNAPS)
end
TS = Dict(SCEN[ii].key => tsruns[ii] for ii in eachindex(SCEN))
TSSNAP = Dict((ii.key, iii) => TS[ii.key][snap_at(TS[ii.key], iii)]
              for ii in SCEN for iii in TS_HOURS)
println("panel C: $(length(SCEN)) runs in $(round(time() - t1; digits=1)) s")

const N_CHECK = 200_000
q_check(Qc) = Qc * N_CELLS / N_CHECK
t1b = time()
ckruns = Vector{Any}(undef, length(SCEN))
Threads.@threads for ii in eachindex(SCEN)
    ckruns[ii] = run_one(TUMOR_R, q_check(SCEN[ii].Qc), GIF_SEED;
                         nframes=length(TS_HOURS) - 1, ncells=N_CHECK)
end
CK = Dict(SCEN[ii].key => ckruns[ii] for ii in eachindex(SCEN))
CKSNAP = Dict((ii.key, iii) => CK[ii.key][snap_at(CK[ii.key], iii)]
              for ii in SCEN for iii in TS_HOURS)
println("panel C check: $(length(SCEN)) runs at $(N_CHECK) particles in ",
        round((time() - t1b) / 60; digits=1), " min")
for ii in SCEN, iii in (24, 48, 72)
    println("   $(rpad(ii.short, 18)) $(lpad(iii, 2)) h  recruited fraction ",
            round(TSSNAP[(ii.key, iii)].n_rec[1] / N_CELLS; digits=4),
            " at $(N_CELLS) vs ",
            round(CKSNAP[(ii.key, iii)].n_rec[1] / N_CHECK; digits=4),
            " at $(N_CHECK)")
end

t2 = time()
otherR   = [ii for ii in RADII if ii != TUMOR_R]
sizejobs = [(ii, iii) for ii in otherR for iii in SCEN]
sizesnap = Vector{Any}(undef, length(sizejobs))
Threads.@threads for ii in eachindex(sizejobs)
    R, s = sizejobs[ii]
    sizesnap[ii] = run_one(R, s.Qc, GIF_SEED; nframes=1)[end]
end
SNAP = Dict((sizejobs[ii][1], sizejobs[ii][2].key) => sizesnap[ii]
            for ii in eachindex(sizejobs))
for ii in SCEN
    SNAP[(TUMOR_R, ii.key)] = TS[ii.key][end]
end
println("panel D: $(length(sizejobs)) runs in $(round(time() - t2; digits=1)) s")

t2b = time()
cksizejobs = [(ii, iii) for ii in otherR for iii in SCEN]
cksizesnap = Vector{Any}(undef, length(cksizejobs))
Threads.@threads for ii in eachindex(cksizejobs)
    R, s = cksizejobs[ii]
    cksizesnap[ii] = run_one(R, q_check(s.Qc), GIF_SEED; nframes=1, ncells=N_CHECK)[end]
end
CKSNAP_D = Dict((cksizejobs[ii][1], cksizejobs[ii][2].key) => cksizesnap[ii]
                for ii in eachindex(cksizejobs))
for ii in SCEN
    CKSNAP_D[(TUMOR_R, ii.key)] = CK[ii.key][end]
end
println("panel D check: $(length(cksizejobs)) runs at $(N_CHECK) particles in ",
        round((time() - t2b) / 60; digits=1), " min")
for ii in RADII, iii in SCEN
    println("   R = $(rpad(Int(ii), 5)) $(rpad(iii.short, 18)) recruited fraction ",
            round(SNAP[(ii, iii.key)].n_rec[1] / N_CELLS; digits=4), " at $(N_CELLS) vs ",
            round(CKSNAP_D[(ii, iii.key)].n_rec[1] / N_CHECK; digits=4), " at $(N_CHECK)")
end

n_producing(snap, R, secretes) =
    count(c -> snap.home[c] > 0 && secretes &&
          hypot(snap.pos[c,1], snap.pos[c,2]) <= R, 1:size(snap.pos, 1))
for ii in RADII, iii in SCEN
    sn = SNAP[(ii, iii.key)]
    println("    R = $(Int(ii)) µm  $(rpad(iii.lab, 36)) recruited ", sn.n_rec[1],
            " / ", N_CELLS, "  producing ", n_producing(sn, ii, iii.Qc > 0))
end

DENS = Dict{Tuple{Float64,String},Vector{Float64}}()
DMM  = Dict{Tuple{Float64,String},Vector{Float64}}()
open(resultpath("tumor_size_profiles.csv")) do io
    readline(io)
    for ii in eachline(io)
        f = split(strip(ii), ',')
        key = (parse(Float64, f[1]), String(f[2]))
        haskey(DENS, key) || (DENS[key] = Float64[]; DMM[key] = Float64[])
        push!(DMM[key],  parse(Float64, f[3]))
        push!(DENS[key], parse(Float64, f[4]))
    end
end

const SKIP = 10
println("panel E: read $(length(DENS)) curves from tumor_size_profiles.csv")

function density_map(pos, zoom)
    n = max(2, round(Int, 2zoom / BIN_UM))
    H = zeros(n, n)
    binmm2 = (2zoom / n / 1000)^2
    @inbounds for ii in axes(pos, 1)
        x = pos[ii, 1]; y = pos[ii, 2]
        (-zoom <= x < zoom && -zoom <= y < zoom) || continue
        ix = floor(Int, (x + zoom) / (2zoom) * n) + 1
        iy = floor(Int, (y + zoom) / (2zoom) * n) + 1
        H[iy, ix] += 1.0
    end
    c = [(-zoom + (ii - 0.5) * 2zoom / n) for ii in 1:n]
    (c, H ./ binmm2)
end

const DEMO_N     = 5000
const DEMO_SIGMA = 900.0
const DEMO_PSF   = 2.0
const DEMO_ANGLES = [0, 45, 90, 135]

function render_blobs(pos, sigma_px)
    img = zeros(NPOS, NPOS)
    c = CENTRE + 1
    rad = ceil(Int, 3 * sigma_px)
    nrm = 1 / (2π * sigma_px^2)
    @inbounds for ii in axes(pos, 1)
        cx = pos[ii, 1] / PX_UM + c; cy = pos[ii, 2] / PX_UM + c
        i0 = round(Int, cx); j0 = round(Int, cy)
        for iii in (i0 - rad):(i0 + rad), iiii in (j0 - rad):(j0 + rad)
            (1 <= iii <= NPOS && 1 <= iiii <= NPOS) || continue
            d2 = (iii - cx)^2 + (iiii - cy)^2
            img[iii, iiii] += nrm * exp(-d2 / (2 * sigma_px^2))
        end
    end
    img
end

demo_rng = MersenneTwister(20240901)
demo_pos = DEMO_SIGMA .* randn(demo_rng, DEMO_N, 2)
demo_img = render_blobs(demo_pos, DEMO_PSF)
demo_single = Dict(ii => line_at(demo_img, ii) ./ PX_AREA_MM2 for ii in DEMO_ANGLES)
demo_avg    = line_mean(demo_img) ./ PX_AREA_MM2
const DEMO_MM = [px_to_mm(ii - CENTRE) for ii in 0:(NPOS - 1)]
println("panel A: $(DEMO_N) synthetic cells, σ = $(Int(DEMO_SIGMA)) µm, ",
        "peak of the 180° average = $(round(maximum(demo_avg); digits=1)) cells/mm², ",
        "single-angle peak = $(round(maximum(demo_single[0]); digits=1))")

COND_LABEL = Dict("Par" => "Off-target (parental)", "CD19" => "Target (CD19)")
comma(n) = replace(string(n), r"(?<=\d)(?=(\d{3})+$)" => ",")

const DENSMAP = [PARULA_W[ii] for ii in range(0, 1; length = 256)]

const FIG_W = 2400
const WA    = 0.42
const WLAB  = 0.075
const WLEG  = 0.70

const SIDE_A = 310
const SIDE_C = 205
const SIDE_D = 268

const TITLE_H = 66
const ROWGAP  = 6
const RH = [70,
            350,
            350,
            100,
            250, 250,
            250, 250,
            250, 250,
            100,
            315, 315,
            315, 315,
            70,
            360]

const FIG_H = TITLE_H + sum(RH) + ROWGAP * length(RH) + 18

fig = CM.Figure(size = (FIG_W, FIG_H), backgroundcolor = :white,
                figure_padding = (14, 14, 8, 10))
top = fig.layout

function bare_axis(cell)
    ax = CM.Axis(cell; width = CM.Relative(1.0), height = CM.Relative(1.0),
                 xautolimitmargin = (0, 0), yautolimitmargin = (0, 0))
    CM.hidedecorations!(ax); CM.hidespines!(ax)
    ax
end

function block_head!(cell, letter, ttl, entries)
    g = CM.GridLayout(cell)
    CM.Label(g[1, 1], "$(letter)   $(ttl)"; fontsize = 15, font = :bold,
             halign = :left, tellwidth = false)
    ax = bare_axis(g[2, 1])
    CM.xlims!(ax, 0, 1); CM.ylims!(ax, 0, 1)
    n = length(entries)
    for (ii, (fillc, strokec, lab)) in enumerate(entries)
        x0 = (ii - 0.5) / n - min(0.10, 0.45 / n)
        CM.scatter!(ax, [x0], [0.5]; color = fillc, markersize = 15,
                    strokecolor = strokec, strokewidth = 1.2)
        CM.text!(ax, x0 + 0.012, 0.5; text = lab, align = (:left, :center),
                 fontsize = 12)
    end

    CM.rowsize!(g, 1, CM.Fixed(26))
    CM.rowsize!(g, 2, CM.Auto(false, 1))
    CM.rowgap!(g, 2)
    g
end

row_label!(cell, txt; fontsize = 13) =
    CM.Label(cell, txt; fontsize = fontsize, halign = :center, justification = :center)

function view_axis(cell, ttl, zoom, side; titlesize = 10)
    ax = CM.Axis(cell; title = ttl, titlesize = titlesize, titlegap = 4,
                 width = side, height = side,
                 xautolimitmargin = (0, 0), yautolimitmargin = (0, 0))
    CM.hidedecorations!(ax); CM.hidespines!(ax)
    CM.xlims!(ax, -1.02zoom, 1.02zoom); CM.ylims!(ax, -1.02zoom, 1.02zoom)
    ax
end

plot_axis(cell; kw...) =
    CM.Axis(cell; xgridvisible = false, ygridvisible = false,
            topspinevisible = false, rightspinevisible = false,
            xautolimitmargin = (0, 0), yautolimitmargin = (0, 0), kw...)

tumor_ring!(ax, R) = CM.lines!(ax, R .* cos.(range(0, 2π; length = 300)),
                               R .* sin.(range(0, 2π; length = 300));
                               color = :black, linewidth = 2)

function scale_bar!(ax, zoom)
    y0 = -0.90zoom
    CM.lines!(ax, [-1000.0, 0.0], [y0, y0]; color = :black, linewidth = 3)
    CM.text!(ax, 150.0, y0; text = "1 mm", align = (:left, :center), fontsize = 10)
end

function cell_panel!(cell, snap, R, ttl, zoom, side; ms = 2.6, scalebar = false,
                     titlesize = 10)
    ax = view_axis(cell, ttl, zoom, side; titlesize = titlesize)
    n = size(snap.pos, 1)
    naive = [snap.home[ii] == 0 for ii in 1:n]
    rec   = .!naive
    any(naive) && CM.scatter!(ax, snap.pos[naive, 1], snap.pos[naive, 2];
        color = (COL_NAIVE, DOT_ALPHA), markersize = ms, strokewidth = 0)
    any(rec) && CM.scatter!(ax, snap.pos[rec, 1], snap.pos[rec, 2];
        color = (COL_REC, DOT_ALPHA), markersize = ms + 0.4, strokewidth = 0)
    tumor_ring!(ax, R)
    scalebar && scale_bar!(ax, zoom)
    ax
end

function dens_panel!(cell, snap, R, ttl, cmax, zoom, side; titlesize = 10, scale = 1.0)
    ax = view_axis(cell, ttl, zoom, side; titlesize = titlesize)
    centers, H = density_map(snap.pos, zoom)
    n = length(centers)
    w = (centers[2] - centers[1]) * 1.02
    rects = CM.Rect2f[]; cols = eltype(DENSMAP)[]
    @inbounds for ii in 1:n, iii in 1:n
        v = H[iii, ii] * scale
        v > 0 || continue
        f = clamp(v / cmax, 0.0, 1.0)
        f < 1 / 256 && continue
        push!(rects, CM.Rect2f(centers[ii] - w / 2, centers[iii] - w / 2, w, w))
        push!(cols, DENSMAP[clamp(ceil(Int, f * 256), 1, 256)])
    end
    isempty(rects) || CM.poly!(ax, rects; color = cols, strokewidth = 0)
    tumor_ring!(ax, R)
    ax
end

CM.Label(top[1, 1],
    "ELISA-constrained fit and its predictions · $(comma(N_CELLS)) simulated T cells · " *
    "ΔFRO sensing · R² = $(round(R2; digits = 3)) over the six measured profiles";
    fontsize = 20, font = :bold, halign = :center)

CELL_LEGEND = [(COL_NAIVE, COL_NAIVE, "not yet recruited"),
               (COL_REC, COL_REC, "recruited (same colour in both conditions)"),
               (:white, :black, "tumor edge")]
g2 = CM.GridLayout(top[2, 1])
block_head!(g2[1, 1], "A", "How the T cell profile is measured · $(comma(DEMO_N)) synthetic cells",
    [(COL_NAIVE, COL_NAIVE, "synthetic cell"),
     (BP_INDIGO, BP_INDIGO, "the line being sampled"),
     (BP_GRAY_L, BP_GRAY_L, "the other rotations")])
block_head!(g2[1, 2], "B", "The fit: measured profiles and the fitted model",
    [(:black, :black, "measured mean (the fit target)"),
     (BP_GOLD, BP_NAVY, "ELISA-constrained model"),
     (:white, :black, "tumor edge, 1.5 mm")])
CM.colsize!(g2, 1, CM.Auto(false, WA)); CM.colsize!(g2, 2, CM.Auto(false, 1 - WA))

g3 = CM.GridLayout(top[3, 1])
g4 = CM.GridLayout(top[4, 1])

function demo_panel!(cell, ttl, angles)
    lim = px_to_mm(CENTRE) * 1000
    ax = view_axis(cell, ttl, lim / 1.02, SIDE_A; titlesize = 12)
    CM.scatter!(ax, demo_pos[:, 1], demo_pos[:, 2];
                color = (COL_NAIVE, DOT_ALPHA), markersize = 3.2, strokewidth = 0)
    single = length(angles) == 1
    for ii in angles
        a = ii * pi / 180
        CM.lines!(ax, [-lim * cos(a), lim * cos(a)], [-lim * sin(a), lim * sin(a)];
                  color = single ? BP_INDIGO : BP_GRAY_L,
                  linewidth = single ? 3.0 : 1.4)
    end
    CM.scatter!(ax, [0.0], [0.0]; color = :black, markersize = 10,
                strokecolor = :white, strokewidth = 1.5)
    ax
end

demo_panel!(g3[1, 1], "1 · one line through the centre", [20])
demo_panel!(g3[1, 2], "2 · rotate it in 1° steps", collect(0:15:165))

axA3 = plot_axis(g4[1, 1]; title = "3 · average every line profile", titlesize = 12,
    titlegap = 4, xlabel = "Distance from centre (mm)",
    ylabel = "T cell density (cells / mm²)",
    xlabelsize = 12, ylabelsize = 12, xticklabelsize = 11, yticklabelsize = 11)
for ii in DEMO_ANGLES
    CM.lines!(axA3, DEMO_MM, demo_single[ii]; color = BP_GRAY_L, linewidth = 0.9,
              label = "4 of the 180 single angles (clipped at the top)")
end
CM.lines!(axA3, DEMO_MM, demo_avg; color = BP_INDIGO, linewidth = 3.0,
          label = "the 180° average")
CM.xlims!(axA3, DEMO_MM[1], DEMO_MM[end])
CM.ylims!(axA3, 0, maximum(demo_avg) * 1.8)
CM.axislegend(axA3; position = :lt, labelsize = 11, framevisible = false,
              merge = true, padding = (4, 4, 2, 2))

ytopB = ceil(maximum(maximum(OBSHALF[ii]) for ii in PANELS) / 10) * 10 + 5
for ii in CONDS, (iii, iiii) in enumerate(FIT_HOURS)
    k = (ii, iiii)
    g = ii == CONDS[1] ? g3 : g4
    col = (ii == CONDS[1] ? 2 : 1) + iii
    ax = plot_axis(g[1, col]; title = "$(COND_LABEL[ii])  ·  $(iiii) h", titlesize = 12,
        titlegap = 4,
        xlabel = ii == CONDS[end] ? "Distance from tumor centre (mm)" : "",
        ylabel = iiii == FIT_HOURS[1] ? "Mean intensity (a.u.)" : "",
        xlabelsize = 12, ylabelsize = 12, xticklabelsize = 11, yticklabelsize = 11,
        xticks = 0:1:2, yticks = 0:10:Int(ytopB - 5))
    CM.vspan!(ax, 0, px_to_mm(CENTRE_SKIP); color = (CM.RGBf(0.85, 0.85, 0.85), 0.9))
    CM.lines!(ax, px_to_mm.(0:CENTRE), OBSHALF[k]; color = :black, linewidth = 2.0)
    CM.lines!(ax, px_to_mm.(RCENT), fitted(k); color = BP_GOLD, linewidth = 2.6)
    CM.scatter!(ax, px_to_mm.(RCENT), fitted(k); color = BP_GOLD, markersize = 7,
                strokecolor = BP_NAVY, strokewidth = 0.8)
    CM.vlines!(ax, [TUMOR_R / 1000]; color = :black, linewidth = 1.2, linestyle = :dot)
    CM.xlims!(ax, 0, px_to_mm(CENTRE)); CM.ylims!(ax, 0, ytopB)
    CM.text!(ax, 0.12, ytopB * 0.93; text = "R² = $(round(r2_one(k); digits = 3))",
             align = (:left, :center), fontsize = 11)
    CM.text!(ax, 0.12, ytopB * 0.83;
             text = "recruited $(comma(round(Int, NRECF[k]))) of $(comma(N_CELLS))",
             align = (:left, :center), fontsize = 11)
end
wfit = (1 - WA) / 3
CM.colsize!(g3, 1, CM.Auto(false, WA / 2)); CM.colsize!(g3, 2, CM.Auto(false, WA / 2))
for ii in 3:5; CM.colsize!(g3, ii, CM.Auto(false, wfit)); end
CM.colsize!(g4, 1, CM.Auto(false, WA))
for ii in 2:4; CM.colsize!(g4, ii, CM.Auto(false, wfit)); end

CMAX_C = maximum(maximum(density_map(TSSNAP[(ii.key, iii)].pos, ZOOM_C)[2])
                 for ii in SCEN for iii in TS_HOURS)
println("panel C density scale: 0 to $(round(CMAX_C; digits=0)) cells/mm² ",
        "($(Int(BIN_UM)) µm bins, shared by all $(2 * 2 * length(TS_HOURS)) density ",
        "panels: $(length(SCEN) * length(TS_HOURS)) at $(N_CELLS) particles and the ",
        "same number at $(N_CHECK))")
CMAX_D = maximum(maximum(density_map(SNAP[(ii, iii.key)].pos, ZOOM_D)[2])
                 for ii in RADII for iii in SCEN)
println("panel D density scale: 0 to $(round(CMAX_D; digits=0)) cells/mm² ",
        "($(Int(BIN_UM)) µm bins, shared by all $(2 * 2 * length(RADII)) panels)")

function density_bar!(cell, cmax; nstep = 256)
    ax = CM.Axis(cell; xlabel = "T cell density (cells / mm²)", xlabelsize = 12,
        xticklabelsize = 11, xlabelpadding = 2, xticksize = 4,
        yticksvisible = false, yticklabelsvisible = false, ygridvisible = false,
        xgridvisible = false, height = 20, width = CM.Relative(0.80),
        halign = :center, valign = :center,
        xautolimitmargin = (0, 0), yautolimitmargin = (0, 0))
    rects = [CM.Rect2f((ii - 1) / nstep * cmax, 0.0, cmax / nstep, 1.0) for ii in 1:nstep]
    cols  = [DENSMAP[clamp(ceil(Int, (ii - 0.5) / nstep * 256), 1, 256)] for ii in 1:nstep]
    CM.poly!(ax, rects; color = cols, strokewidth = 0)
    CM.xlims!(ax, 0, cmax); CM.ylims!(ax, 0, 1)
    ax
end

function head_with_bar!(row, letter, ttl, cmax)
    g = CM.GridLayout(top[row, 1])
    block_head!(g[1, 1], letter, ttl, CELL_LEGEND)
    density_bar!(g[1, 2], cmax)
    CM.colsize!(g, 1, CM.Auto(false, WLEG)); CM.colsize!(g, 2, CM.Auto(false, 1 - WLEG))

    CM.rowsize!(g, 1, CM.Relative(1.0))
    g
end

head_with_bar!(5, "C", "The population over 72 h at the fitted tumor size, " *
    "R = 1500 µm · rows 1–2 cells, rows 3–4 the same cells as density, " *
    "rows 5–6 the same run at one particle per real T cell", CMAX_C)

crows = [(SCEN[1], :cells), (SCEN[2], :cells), (SCEN[1], :dens), (SCEN[2], :dens)]
for (ii, (s, kind)) in enumerate(crows)
    g = CM.GridLayout(top[5 + ii, 1])
    row_label!(g[1, 1], (kind == :cells ? "Cells\n" : "Density\n") * s.rowlab)
    for (iii, iiii) in enumerate(TS_HOURS)
        snap = TSSNAP[(s.key, iiii)]
        if kind == :cells
            cell_panel!(g[1, iii + 1], snap, TUMOR_R,
                "$(iiii) h\nrecruited $(comma(snap.n_rec[1]))", ZOOM_C, SIDE_C;
                ms = 2.4, scalebar = (ii == 1 && iii == 1))
        else
            dens_panel!(g[1, iii + 1], snap, TUMOR_R, "$(iiii) h", CMAX_C, ZOOM_C, SIDE_C)
        end
    end
    CM.colsize!(g, 1, CM.Auto(false, WLAB))
    for iii in 2:(length(TS_HOURS) + 1)
        CM.colsize!(g, iii, CM.Auto(false, (1 - WLAB) / length(TS_HOURS)))
    end
    CM.colgap!(g, 4)
end

for (ii, s) in enumerate(SCEN)
    g = CM.GridLayout(top[9 + ii, 1])
    row_label!(g[1, 1], "Density, 200k check\n" * s.rowlab)
    for (iii, iiii) in enumerate(TS_HOURS)
        dens_panel!(g[1, iii + 1], CKSNAP[(s.key, iiii)], TUMOR_R, "$(iiii) h", CMAX_C,
                    ZOOM_C, SIDE_C; scale = N_CELLS / N_CHECK)
    end
    CM.colsize!(g, 1, CM.Auto(false, WLAB))
    for iii in 2:(length(TS_HOURS) + 1)
        CM.colsize!(g, iii, CM.Auto(false, (1 - WLAB) / length(TS_HOURS)))
    end
    CM.colgap!(g, 4)
end

head_with_bar!(12, "D", "Four tumor sizes at 72 h · left half cells, right half " *
    "the same cells as density · rows 3–4 repeat the density at one particle " *
    "per real T cell", CMAX_D)

for (ii, s) in enumerate(SCEN)
    g = CM.GridLayout(top[12 + ii, 1])
    row_label!(g[1, 1], s.rowlab)
    for (iii, iiii) in enumerate(RADII)
        snap = SNAP[(iiii, s.key)]
        cell_panel!(g[1, iii + 1], snap, iiii,
            "R = $(Int(iiii)) µm\ncells · recruited $(comma(snap.n_rec[1]))", ZOOM_D,
            SIDE_D; ms = 2.8, scalebar = (ii == 1 && iii == 1))
        dens_panel!(g[1, iii + 1 + length(RADII)], snap, iiii,
            "R = $(Int(iiii)) µm\ndensity", CMAX_D, ZOOM_D, SIDE_D)
    end
    CM.colsize!(g, 1, CM.Auto(false, WLAB))
    for iii in 2:(2 * length(RADII) + 1)
        CM.colsize!(g, iii, CM.Auto(false, (1 - WLAB) / (2 * length(RADII))))
    end
    CM.colgap!(g, 4)
end

for (ii, s) in enumerate(SCEN)
    g = CM.GridLayout(top[14 + ii, 1])
    row_label!(g[1, 1], "Density, 200k check\n" * s.rowlab)
    for (iii, iiii) in enumerate(RADII)
        dens_panel!(g[1, iii + 1 + length(RADII)], CKSNAP_D[(iiii, s.key)], iiii,
            "R = $(Int(iiii)) µm\ndensity, 200k particles", CMAX_D, ZOOM_D, SIDE_D;
            scale = N_CELLS / N_CHECK)
    end
    CM.colsize!(g, 1, CM.Auto(false, WLAB))
    for iii in 2:(2 * length(RADII) + 1)
        CM.colsize!(g, iii, CM.Auto(false, (1 - WLAB) / (2 * length(RADII))))
    end
    CM.colgap!(g, 4)
end

g13 = CM.GridLayout(top[17, 1])
block_head!(g13[1, 1], "E", "Predicted T cell density at 72 h against distance, " *
    "across those four sizes",
    [(SCEN[1].col, SCEN[1].col, SCEN[1].lab),
     (SCEN[2].col, SCEN[2].col, SCEN[2].lab),
     (:white, :black, "tumor edge")])

ytopE = ceil(maximum(maximum(DENS[(ii, iii.key)][SKIP+1:end]) for ii in RADII
                     for iii in SCEN) / 500) * 500
g14 = CM.GridLayout(top[18, 1])
for (ii, iii) in enumerate(RADII)
    ax = plot_axis(g14[1, ii];
        title = "Tumor radius $(Int(iii)) µm" * (iii == TUMOR_R ? "  (the fitted size)" : ""),
        titlesize = 13, titlegap = 4,
        xlabel = "Distance from tumor centre (mm)",
        ylabel = ii == 1 ? "T cell density (cells / mm²)" : "",
        xlabelsize = 13, ylabelsize = 13, xticklabelsize = 12, yticklabelsize = 12,
        xticks = 0:1:5, yticks = 0:500:Int(ytopE))
    for iiii in SCEN
        CM.lines!(ax, DMM[(iii, iiii.key)][SKIP+1:end], DENS[(iii, iiii.key)][SKIP+1:end];
                  color = iiii.col, linewidth = 2.6)
    end
    CM.vlines!(ax, [iii / 1000]; color = :black, linewidth = 1.2, linestyle = :dot)
    CM.xlims!(ax, 0, HALF_MM); CM.ylims!(ax, 0, ytopE)
end
for ii in 1:length(RADII); CM.colsize!(g14, ii, CM.Auto(false, 1)); end
CM.colgap!(g14, 26)

CM.colsize!(top, 1, CM.Fixed(FIG_W - 28))

CM.rowsize!(top, 1, CM.Fixed(TITLE_H))
for (ii, iii) in enumerate(RH); CM.rowsize!(top, ii + 1, CM.Fixed(iii)); end
CM.rowgap!(top, ROWGAP)

CM.resize_to_layout!(fig)

function save_both_cm(f, name)
    CM.save(figpath(name * ".png"), f; px_per_unit = 1)
    CM.save(figpath(name * ".pdf"), f; pt_per_unit = 1)
    (png = name * ".png", pdf = name * ".pdf")
end

save_both_cm(fig, "summary_figure")
let sz = size(fig.scene)
    println("SAVED summary_figure.png and summary_figure.pdf  ",
            "($(round(Int, sz[1])) × $(round(Int, sz[2])) px, ",
            "CairoMakie: embedded fonts, no rasters)")
end
