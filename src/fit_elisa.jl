include("model.jl")
include("params.jl")

const WELL_R   = well_radius_um("24-well")
const TUMOR_R  = 1500.0
const FIT_T    = 72 * 60.0
const STRAT    = :receptor

const N_CELLS  = parse(Int, get(ENV, "FIT_NCELLS", "10000"))
const PARAMS_CSV = "fit_parameters.csv"
const D_FIXED  = CK1_D
const K_FIXED  = CK_K_MEASURED
const LAMBDA   = sqrt(D_FIXED / K_FIXED)
const PX_UM    = 14000.0 / 1188.0
const CELL_R   = CELL_DIAMETER_UM / 2
px_to_mm(p) = p * PX_UM / 1000
const GRID_UM  = 50.0
const GRID_R   = 2600.0
const CENTRE   = 250
const NPOS     = 501
const CENTRE_SKIP = 10
const RBIN     = 5
const FIT_HOURS = [24, 47, 72]
const SNAP_MIN  = 60.0
const N_SNAPS   = round(Int, FIT_T / SNAP_MIN)

const MW_CXCL12 = 7963.0
const MW_CCL20  = 8000.0
const M_TUMOR = (1989.4816666666666 * 2e-12) / MW_CXCL12
const M_CELL  = (52.48407 * 0.1e-12)      / MW_CCL20
const N_TUMOR_WELL = 4.0e5;  const T_TUMOR_WELL = 48 * 60.0
const N_CELL_WELL  = 5.0e4;  const T_CELL_WELL  = 72 * 60.0
elisa_rate(M, N, T, k) = M * k / (N * (1 - exp(-k * T)))

const N_TUMOR_CELLS = 35_000
const SCALE = 200_000 / N_CELLS
const R_TUMOR = elisa_rate(M_TUMOR, N_TUMOR_WELL, T_TUMOR_WELL, K_FIXED)

const SIGMA_TUMOR = N_TUMOR_CELLS * R_TUMOR / (pi * TUMOR_R^2)
const MOL_PER_UM3_TO_NM = 1e24
const C_TUMOR_SURFACE = SIGMA_TUMOR * LAMBDA / (2 * D_FIXED) * MOL_PER_UM3_TO_NM
const Q_PRIM = C_TUMOR_SURFACE * ISLAND_REF_AREA * K_FIXED

const Q_PRIM_TOTAL = Q_PRIM / ISLAND_REF_AREA * (pi * TUMOR_R^2)

q_cell(k_cell) = Q_PRIM_TOTAL * SCALE *
                 elisa_rate(M_CELL, N_CELL_WELL, T_CELL_WELL, k_cell) /
                 (N_TUMOR_CELLS * R_TUMOR)

const NPOP  = parse(Int, get(ENV, "FIT_NPOP", "12"))
const NGEN  = parse(Int, get(ENV, "FIT_NGEN", "20"))
const FSEEDS = [11, 23, 37, 49]
const FINAL_SEEDS = [101, 113, 127, 139, 151, 163, 177, 191]

x_arr, dx, nx, mask = build_well_grid(GRID_R, GRID_UM)
center = reshape([0.0, 0.0], 1, 2)

println("TUMOR SIDE, derived from the CXCL12 ELISA:")
println("  per tumor cell     = $(round(R_TUMOR; sigdigits=4)) mol/min at k = ",
        "$(round(K_FIXED; sigdigits=4)) /min")
println("  sheet secretion σ  = $(round(SIGMA_TUMOR; sigdigits=4)) mol/(min·µm²) ",
        "over $(N_TUMOR_CELLS) cells")
println("  surface conc σλ/2D = $(round(C_TUMOR_SURFACE; digits=2)) nM")
println("  Q_prim             = $(round(Q_PRIM; sigdigits=6)) nM·µm²/min")
println("  peak CXCL12 deep inside the disk = ",
        round(Q_PRIM / ISLAND_REF_AREA / K_FIXED; digits=2), " nM (Kd = $(CK1_KD) nM)")
println("Q_cell is derived from k_cell, not fitted. Over the range the fit may try:")
for ii in (2.0, 26.5, 120.0, 720.0, 10080.0)
    kc = log(2) / ii
    println("   half-life $(rpad(ii, 8)) min  ->  k_cell = $(rpad(round(kc; sigdigits=3), 10)) ",
            "λ_cell = $(rpad(round(Int, sqrt(D_FIXED/kc)), 7)) µm  Q_cell = ",
            round(q_cell(kc); sigdigits=4))
end

rawmean = Dict{Tuple{String,Int}, Vector{Float64}}()
open(datapath("infiltration_measured_mean.csv")) do io
    readline(io)
    for ii in eachline(io)
        fields = split(strip(ii), ',')
        hr = parse(Int, fields[2]); hr in FIT_HOURS || continue
        key = (String(fields[1]), hr)
        haskey(rawmean, key) || (rawmean[key] = Float64[])
        push!(rawmean[key], parse(Float64, fields[4]))
    end
end
CONDS = ["Par", "CD19"]
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

function run_one(θ, Qcell, seed; nframes=N_SNAPS, ncells=N_CELLS)
    kc, kda, bwp, bwa, rd = θ
    run_multisource_recruitment(;
        well_radius=WELL_R, island_centers=center, island_radius=TUMOR_R,
        island_Q=[Qcell], n_cells=ncells, D=D_FIXED, k=K_FIXED, Kd=CK1_KD,
        ck_sensitivity=10.0, cell_speed=CELL_SPEED, rot_diffusion=rd,
        sensing_radius=SENSING_RADIUS, bias_weight=bwp, bias_auto=bwa,
        k_auto=kc, Kd_auto=kda, dt=DT, total_time=FIT_T, n_snapshots=nframes,
        mask=mask, x_arr=x_arr, dx=dx, nx=nx, strategy=STRAT,
        n_seed_per_island=0, island_primary_Q=[Q_PRIM], seed=seed)
end

snap_at(snaps, hr) = argmin(abs(ii.t - hr * 60.0) for ii in snaps)

function sim_profiles(θ; seeds=FSEEDS, ncells=N_CELLS)
    Qc = q_cell(θ[1])
    jobs = [(key, Q, iii) for (key, Q) in (("Par", 0.0), ("CD19", Qc)) for iii in seeds]
    prof = Vector{Dict{Int,Vector{Float64}}}(undef, length(jobs))
    nrec = Vector{Dict{Int,Float64}}(undef, length(jobs))
    Threads.@threads for ii in eachindex(jobs)
        key, Q, sd = jobs[ii]
        img = zeros(NPOS, NPOS)
        snaps = run_one(θ, Q, sd; ncells=ncells)
        p = Dict{Int,Vector{Float64}}(); n = Dict{Int,Float64}()
        for iii in FIT_HOURS
            f = snaps[snap_at(snaps, iii)]
            p[iii] = line_profile(render!(img, f.pos))
            n[iii] = f.n_rec[1]
        end
        prof[ii] = p; nrec[ii] = n
    end
    out  = Dict{Tuple{String,Int},Vector{Float64}}()
    half = Dict{Tuple{String,Int},Vector{Float64}}()
    rec  = Dict{Tuple{String,Int},Float64}()
    for (key, _) in (("Par", 0.0), ("CD19", Qc))
        idx = [iii for iii in eachindex(jobs) if jobs[iii][1] == key]
        for iii in FIT_HOURS
            acc = sum(prof[iiii][iii] for iiii in idx) ./ length(idx)
            half[(key, iii)] = acc
            out[(key, iii)]  = [binmean(acc, iiii) for iiii in RBINS]
            rec[(key, iii)]  = mean(nrec[iiii][iii] for iiii in idx)
        end
    end
    (out, rec, half)
end

const PVAR = Dict(ii => sum((OBSBIN[ii] .- mean(OBSBIN[ii])) .^ 2) for ii in PANELS)
const PWT  = vcat((fill(1 / PVAR[ii], length(OBSBIN[ii])) for ii in PANELS)...)
const SQWT = sqrt.(PWT)

function solve_ab(sim)
    d = vcat((sim[ii] for ii in PANELS)...)
    m = vcat((OBSBIN[ii] for ii in PANELS)...)
    X = hcat(d, ones(length(d)))
    ab = (SQWT .* X) \ (SQWT .* m)
    a, b = ab[1], ab[2]
    a < 0 && (a = 0.0; b = sum(PWT .* m) / sum(PWT))
    (a, b, sum(PWT .* (m .- (a .* d .+ b)) .^ 2))
end

nevals = Ref(0)
function objective(θ)
    nevals[] += 1
    solve_ab(sim_profiles(θ)[1])[3]
end

const PNAMES = ["k_cell (1/min)", "Kd_auto (nM)", "bias_prim (tumor chemokine)",
                "bias_auto (orthokine)", "rot_diffusion (rad²/min)"]
const LOGP = [true, true, true, true, false]
const LO = [log10(log(2) / 10080.0), log10(0.01), log10(0.05), log10(0.05), 0.005]
const HI = [log10(log(2) / 2.0),     log10(100.0), log10(500.0), log10(500.0), 1.00]
decode(z) = [LOGP[ii] ? 10.0^z[ii] : z[ii] for ii in eachindex(z)]
encode(θ) = [LOGP[ii] ? log10(θ[ii]) : θ[ii] for ii in eachindex(θ)]

function differential_evolution(f, lo, hi; npop, ngen, F=0.7, CR=0.9,
                                rng=MersenneTwister(7))
    n = length(lo)
    pop = [lo .+ (hi .- lo) .* rand(rng, n) for ii in 1:npop]
    fit = [f(ii) for ii in pop]
    println("  gen  0: best = $(round(minimum(fit); digits=3))")
    for ii in 1:ngen
        for iii in 1:npop
            d1, d2, d3 = iii, iii, iii
            while d1 == iii; d1 = rand(rng, 1:npop); end
            while d2 == iii || d2 == d1; d2 = rand(rng, 1:npop); end
            while d3 == iii || d3 == d1 || d3 == d2; d3 = rand(rng, 1:npop); end
            mut = pop[d1] .+ F .* (pop[d2] .- pop[d3])
            jr = rand(rng, 1:n)
            trial = [(rand(rng) < CR || iiii == jr) ? mut[iiii] : pop[iii][iiii]
                     for iiii in 1:n]
            for iiii in 1:n
                trial[iiii] < lo[iiii] && (trial[iiii] = lo[iiii] + (lo[iiii] - trial[iiii]))
                trial[iiii] > hi[iiii] && (trial[iiii] = hi[iiii] - (trial[iiii] - hi[iiii]))
                trial[iiii] = clamp(trial[iiii], lo[iiii], hi[iiii])
            end
            ft = f(trial)
            ft < fit[iii] && (pop[iii] = trial; fit[iii] = ft)
        end
        bi = argmin(fit)
        println("  gen $(lpad(ii, 2)): best = $(round(fit[bi]; digits=3))   ",
                "θ = ", round.(decode(pop[bi]); sigdigits=3))
    end
    bi = argmin(fit)
    (pop[bi], fit[bi])
end

function compass_polish(f, z0, lo, hi; step0=0.06, stepmin=0.004, maxeval=120)
    n = length(z0)
    u = (z0 .- lo) ./ (hi .- lo)
    tou(v) = lo .+ clamp.(v, 0.0, 1.0) .* (hi .- lo)
    best = f(tou(u)); ne = 1; step = step0
    while step >= stepmin && ne < maxeval
        improved = false
        for ii in 1:n, iii in (+1, -1)
            v = copy(u); v[ii] = clamp(v[ii] + iii * step, 0.0, 1.0)
            v[ii] == u[ii] && continue
            fv = f(tou(v)); ne += 1
            fv < best && (u = v; best = fv; improved = true)
            ne >= maxeval && break
        end
        improved || (step /= 2)
    end
    println("  polish: $(ne) evaluations, objective $(round(best; digits=3))")
    (tou(u), best)
end

t0 = time()
obj_z(z) = objective(decode(z))

const REUSE = get(ENV, "FIT_REUSE", "0") == "1" && isfile(resultpath(PARAMS_CSV))
local zbest, fbest
if REUSE
    vals = Dict{String,Float64}()
    open(resultpath(PARAMS_CSV)) do io
        readline(io)
        for ii in eachline(io)
            at = findlast(',', ii); at === nothing && continue
            vals[ii[1:at-1]] = parse(Float64, ii[at+1:end])
        end
    end
    zbest = encode([vals[ii] for ii in PNAMES]); fbest = obj_z(zbest)
    println("\n=== FIT_REUSE=1: reusing the stored parameters, no fit ===")
    println("  objective at the stored parameters: ", round(fbest; digits=3))
else
    println("\n=== differential evolution: $(NPOP) × $(NGEN), 5 free parameters, ",
            "$(length(FSEEDS)) seeds per evaluation ===")
    zde, fde = differential_evolution(obj_z, LO, HI; npop=NPOP, ngen=NGEN)
    println("\n=== local refinement ===")
    zbest, fbest = compass_polish(obj_z, zde, LO, HI)
    fbest > fde && (zbest = zde; fbest = fde)
end
θbest = decode(zbest)
println("\nfit done: $(nevals[]) evaluations in $(round((time()-t0)/60; digits=1)) min")

simF, recF, halfF = sim_profiles(θbest; seeds=FINAL_SEEDS)
aF, bF, _ = solve_ab(simF)
fitted(k) = aF .* simF[k] .+ bF
allobs = vcat((OBSBIN[ii] for ii in PANELS)...)
allfit = vcat((fitted(ii) for ii in PANELS)...)
R2 = 1 - sum((allobs .- allfit) .^ 2) / sum((allobs .- mean(allobs)) .^ 2)
r2_one(k) = (o = OBSBIN[k]; f = fitted(k);
             1 - sum((o .- f) .^ 2) / sum((o .- mean(o)) .^ 2))
function r2_cond(c)
    o = vcat((OBSBIN[(c, ii)] for ii in FIT_HOURS)...)
    f = vcat((fitted((c, ii)) for ii in FIT_HOURS)...)
    1 - sum((o .- f) .^ 2) / sum((o .- mean(o)) .^ 2)
end

const QC_BEST = q_cell(θbest[1])
println("\n=== best-fit parameters, 5 free, orthokine production ELISA-constrained ===")
for (ii, nm) in enumerate(PNAMES)
    println("  $(rpad(nm, 30)) = ", round(θbest[ii]; sigdigits=4))
end
println("  $(rpad("→ CCL20 half-life", 30)) = $(round(log(2)/θbest[1]; digits=1)) min")
println("  $(rpad("→ λ_cell = √(D/k_cell)", 30)) = ",
        "$(round(sqrt(D_FIXED/θbest[1]); digits=0)) µm")
println("  $(rpad("→ Q_cell from the ELISA", 30)) = $(round(QC_BEST; sigdigits=4))")
for (ii, nm) in enumerate(PNAMES)
    z = zbest[ii]; span = HI[ii] - LO[ii]
    atlo = z <= LO[ii] + 0.01span; athi = z >= HI[ii] - 0.01span
    (atlo || athi) && println("  ⚠ AT BOUND: $(nm) is at its ",
                              atlo ? "LOWER" : "UPPER", " bound; not an estimate")
end
println("  $(rpad("mean R² over the 6 profiles", 30)) = ",
        round(mean(r2_one(ii) for ii in PANELS); digits=4),
        "  (the objective is Σ(1−R²) = ",
        round(6 - sum(r2_one(ii) for ii in PANELS); digits=3), ")")
println("  $(rpad("intensity gain a", 30)) = $(round(aF; sigdigits=4))")
println("  $(rpad("intensity background b", 30)) = $(round(bF; digits=3))")
println("  R² over all $(length(PANELS)) profiles (held-out seeds) = $(round(R2; digits=4))")
for ii in CONDS
    println("    $(rpad(ii, 5)) pooled = $(round(r2_cond(ii); digits=3))   ",
            join(["$(iii)h $(round(r2_one((ii, iii)); digits=3))" for iii in FIT_HOURS],
                 "  "))
end
for ii in FIT_HOURS
    println("  recruited at $(ii) h: Par $(round(recF[("Par", ii)]; digits=1)), ",
            "CD19 $(round(recF[("CD19", ii)]; digits=1)) of $(N_CELLS)")
end

open(resultpath(PARAMS_CSV), "w") do io
    println(io, "parameter,value")
    for (ii, nm) in enumerate(PNAMES); println(io, nm, ",", θbest[ii]); end
    println(io, "CCL20 half-life (min) derived,", log(2) / θbest[1])
    println(io, "lambda_cell (um) derived,", sqrt(D_FIXED / θbest[1]))
    println(io, "Q_cell (nM um2/min) ELISA-CONSTRAINED,", QC_BEST)
    println(io, "Q_prim (nM um2/min) FIXED,", Q_PRIM)
    println(io, "Q_prim from CXCL12 ELISA surface conc (nM),", C_TUMOR_SURFACE)
    println(io, "peak CXCL12 deep in disk (nM),", Q_PRIM / ISLAND_REF_AREA / K_FIXED)
    println(io, "k tumor (1/min) FIXED,", K_FIXED)
    println(io, "Kd tumor (nM) FIXED,", CK1_KD)
    println(io, "D (um2/min) FIXED,", D_FIXED)
    println(io, "tumor cells assumed in the disk,", N_TUMOR_CELLS)
    println(io, "intensity gain a,", aF); println(io, "intensity background b,", bF)
    println(io, "R2 all profiles,", R2)
    for ii in CONDS
        println(io, "R2 $(ii) pooled,", r2_cond(ii))
        for iii in FIT_HOURS; println(io, "R2 $(ii) $(iii)h,", r2_one((ii, iii))); end
    end
    for ii in CONDS, iii in FIT_HOURS
        println(io, "recruited $(ii) $(iii)h,", recF[(ii, iii)])
    end
end
println("SAVED ", PARAMS_CSV)
