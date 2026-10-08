ENV["GKSwstype"] = "100"

using Plots, LinearAlgebra, Random, Statistics
gr()

const PARULA = cgrad([
    RGB(0.2078, 0.1647, 0.5294), RGB(0.0118, 0.3882, 0.8824),
    RGB(0.0784, 0.5216, 0.8314), RGB(0.0235, 0.6549, 0.7765),
    RGB(0.2196, 0.7255, 0.6196), RGB(0.5725, 0.7490, 0.4510),
    RGB(0.8510, 0.7294, 0.3373), RGB(0.9882, 0.8078, 0.1804),
    RGB(0.9765, 0.9843, 0.0549),
])
parula(f) = PARULA[clamp(float(f), 0.0, 1.0)]

const PC_DARK = parula(0.0)
const PC_TEAL = parula(0.45)
const PC_GOLD = parula(0.80)
const PC_YELLOW = parula(1.0)

const BP_INDIGO = RGB(81/255,  95/255, 170/255)
const BP_GOLD   = RGB(255/255, 195/255, 57/255)
const BP_JADE   = RGB(75/255, 189/255, 148/255)
const BP_CYAN   = RGB(11/255, 176/255, 213/255)
const BP_AZURE  = RGB(11/255, 142/255, 206/255)
const BP_LIME   = RGB(172/255, 200/255, 55/255)
const BP_NAVY   = RGB(0/255,   64/255, 113/255)
const BP_TEAL   = RGB(11/255, 128/255, 131/255)
const BP_GRAY   = RGB(127/255, 127/255, 127/255)
const BP_GRAY_L = RGB(153/255, 153/255, 153/255)

const PARULA_W = cgrad([RGB(1.0, 1.0, 1.0), parula(0.25), parula(0.45),
                        parula(0.62), parula(0.80), parula(1.0)])

const ISLAND_COLORS = [parula(0.0), parula(0.82), parula(0.5), parula(0.30), parula(1.0)]

function legend_strip(entries; width=900, height=46, msize=9, fontsize=9)
    p = plot(; framestyle=:none, legend=false, grid=false, ticks=false,
        xlim=(0, 1), ylim=(0, 1), size=(width, height), margin=0Plots.mm)
    n = length(entries)
    for (ii, (fill, stroke, lab)) in enumerate(entries)
        x0 = (ii - 0.5) / n - min(0.10, 0.45/n)
        scatter!(p, [x0], [0.5]; color=fill, ms=msize, markerstrokecolor=stroke,
            markerstrokewidth=1.0, label="")
        annotate!(p, x0 + 0.018, 0.5, text("  " * lab, fontsize, :left, :black))
    end
    p
end

const FIGDIR = abspath(joinpath(@__DIR__, "..", "figures"))
isdir(FIGDIR) || mkpath(FIGDIR)
figpath(name) = joinpath(FIGDIR, name)

const DATADIR = abspath(joinpath(@__DIR__, "..", "data"))
datapath(name) = joinpath(DATADIR, name)

const RESULTDIR = abspath(joinpath(@__DIR__, "..", "results"))
isdir(RESULTDIR) || mkpath(RESULTDIR)
resultpath(name) = joinpath(RESULTDIR, name)

function savefig_both(fig, name)
    base = endswith(name, ".png") ? name[1:end-4] : name
    savefig(fig, figpath(base * ".png"))
    savefig(fig, figpath(base * ".pdf"))
    (png = base * ".png", pdf = base * ".pdf")
end

const WELL_DIAMETERS_MM = Dict(
    "6-well" => 34.8, "12-well" => 22.1,
    "24-well" => 15.6, "48-well" => 11.0, "96-well" => 6.4,
)
well_radius_um(fmt) = (WELL_DIAMETERS_MM[fmt] / 2) * 1000.0

function build_well_grid(well_radius, grid_spacing)
    ghalf = ceil(Int, well_radius / grid_spacing)
    nx = 2 * ghalf + 1
    x = collect(range(-ghalf * grid_spacing, ghalf * grid_spacing; length=nx))
    dx = grid_spacing
    mask = falses(nx, nx)
    for ii in 1:nx, iii in 1:nx
        mask[iii, ii] = sqrt(x[iii]^2 + x[ii]^2) <= well_radius
    end
    return x, dx, nx, mask
end

function ar_diffstep!(C, Cn, mask, D, k, S, dx, dt, nx)
    coeff = D * dt / dx^2
    decay = k * dt
    @inbounds for ii in 1:nx, iii in 1:nx
        if !mask[iii,ii]; Cn[iii,ii] = 0.0; continue; end
        CL = (iii > 1  && mask[iii-1,ii]) ? C[iii-1,ii] : C[iii,ii]
        CR = (iii < nx && mask[iii+1,ii]) ? C[iii+1,ii] : C[iii,ii]
        CD = (ii > 1  && mask[iii,ii-1]) ? C[iii,ii-1] : C[iii,ii]
        CU = (ii < nx && mask[iii,ii+1]) ? C[iii,ii+1] : C[iii,ii]
        Cn[iii,ii] = C[iii,ii] + coeff * (CL + CR + CD + CU - 4 * C[iii,ii]) -
                     decay * C[iii,ii] + dt * S[iii,ii]
    end
    C .= Cn
end

function ar_interp(F, xa, dx, px, py, nx)
    ix_f = (px - xa[1]) / dx + 1
    iy_f = (py - xa[1]) / dx + 1
    ix = clamp(floor(Int, ix_f), 1, nx - 1)
    iy = clamp(floor(Int, iy_f), 1, nx - 1)
    fx = clamp(ix_f - ix, 0.0, 1.0)
    fy = clamp(iy_f - iy, 0.0, 1.0)
    @inbounds (1-fx)*(1-fy)*F[ix,iy] + fx*(1-fy)*F[ix+1,iy] +
              (1-fx)*fy*F[ix,iy+1] + fx*fy*F[ix+1,iy+1]
end

function ar_grad!(gx, gy, C, mask, dx, nx)
    inv2dx = 1.0 / (2 * dx)
    @inbounds for ii in 1:nx, iii in 1:nx
        if !mask[iii,ii] || iii <= 1 || iii >= nx || ii <= 1 || ii >= nx
            gx[iii,ii] = 0.0; gy[iii,ii] = 0.0; continue
        end
        gx[iii,ii] = (C[iii+1,ii] - C[iii-1,ii]) * inv2dx
        gy[iii,ii] = (C[iii,ii+1] - C[iii,ii-1]) * inv2dx
    end
end

# mirificus biologicus ludus

const ISLAND_REF_AREA = π * 1000.0^2

function deposit_bilinear!(S, x_arr, dx, px, py, q, nx, mask)
    fx = (px - x_arr[1]) / dx + 1
    fy = (py - x_arr[1]) / dx + 1
    ix = clamp(floor(Int, fx), 1, nx - 1); iy = clamp(floor(Int, fy), 1, nx - 1)
    a = clamp(fx - ix, 0.0, 1.0); b = clamp(fy - iy, 0.0, 1.0)
    s = q / dx^2
    @inbounds begin
        mask[ix,   iy]   && (S[ix,   iy]   += s * (1 - a) * (1 - b))
        mask[ix+1, iy]   && (S[ix+1, iy]   += s * a * (1 - b))
        mask[ix,   iy+1] && (S[ix,   iy+1] += s * (1 - a) * b)
        mask[ix+1, iy+1] && (S[ix+1, iy+1] += s * a * b)
    end
end

@inline function in_island(px, py, island_centers, island_radius, m)
    m > 0 || return false
    (px - island_centers[m,1])^2 + (py - island_centers[m,2])^2 <= island_radius^2
end

@inline function sense_field(Cf, gxf, gyf, x_arr, dx, nx, px, py, sensing_radius, Kd,
                             strategy, ck_sensitivity, Est, Ist, ci, dt, τ_E, τ_I)
    gxi = ar_interp(gxf, x_arr, dx, px, py, nx)
    gyi = ar_interp(gyf, x_arr, dx, px, py, nx)
    gm = sqrt(gxi^2 + gyi^2) + 1e-10
    twx, twy = -gxi / gm, -gyi / gm
    cc = max(ar_interp(Cf, x_arr, dx, px, py, nx), 0.0)
    cf = max(ar_interp(Cf, x_arr, dx, px + sensing_radius * twx,
                                     py + sensing_radius * twy, nx), 0.0)
    cb = max(ar_interp(Cf, x_arr, dx, px - sensing_radius * twx,
                                     py - sensing_radius * twy, nx), 0.0)
    local st
    if strategy == :fractional
        st = (cf - cb) / (cc + Kd)
    elseif strategy == :receptor
        st = cf / (cf + Kd) - cb / (cb + Kd)
    elseif strategy == :legi
        fro = cc / (cc + Kd)
        Est[ci] += dt * (fro - Est[ci]) / τ_E
        Ist[ci] += dt * (fro - Ist[ci]) / τ_I
        st = -max(Est[ci] - Ist[ci], 0.0)
    else
        st = 0.0
    end
    (ck_sensitivity * st * twx, ck_sensitivity * st * twy)
end

function run_multisource_recruitment(;
        well_radius, island_centers::Matrix{Float64}, island_radius::Float64,
        island_Q::Vector{Float64},
        n_cells, D, k, Kd, ck_sensitivity,
        cell_speed, rot_diffusion, sensing_radius, bias_weight,
        dt, total_time, n_snapshots,
        mask, x_arr, dx, nx,
        strategy::Symbol=:fractional,
        n_seed_per_island::Int=0,
        island_primary_Q::Vector{Float64}=zeros(size(island_centers, 1)),
        seed_radius=nothing, bias_auto=nothing, k_auto=nothing, Kd_auto=nothing,
        τ_E=0.5, τ_I=10.0, seed=42)

    Random.seed!(seed)
    M = size(island_centers, 1)

    split_k = k_auto !== nothing
    k_a = split_k ? Float64(k_auto) : Float64(k)
    split_k && bias_auto === nothing &&
        error("k_auto needs two-field sensing: pass bias_auto as well")

    Kd_a = Kd_auto === nothing ? Float64(Kd) : Float64(Kd_auto)
    Kd_auto !== nothing && bias_auto === nothing &&
        error("Kd_auto needs two-field sensing: pass bias_auto as well")

    two_field = bias_auto !== nothing

    r_seed = seed_radius === nothing ? well_radius : seed_radius

    pos  = zeros(n_cells, 2)
    home = zeros(Int, n_cells)

    idx = 0
    for ii in 1:M, iii in 1:n_seed_per_island
        idx >= n_cells && break
        idx += 1
        r = island_radius * sqrt(rand()); θ = 2π * rand()
        pos[idx, :] = island_centers[ii, :] .+ [r * cos(θ), r * sin(θ)]
        home[idx] = ii
    end

    while idx < n_cells
        r = r_seed * sqrt(rand())
        θ = 2π * rand()
        x, y = r * cos(θ), r * sin(θ)
        ok = true
        for ii in 1:M
            if (x - island_centers[ii,1])^2 + (y - island_centers[ii,2])^2 <=
               (island_radius * 1.2)^2
                ok = false; break
            end
        end
        if ok
            idx += 1
            pos[idx, :] = [x, y]
        end
    end
    dirs = 2π .* rand(n_cells)

    primary_nodes = [Tuple{Int,Int}[] for ii in 1:M]
    for ii in 1:M
        island_primary_Q[ii] > 0 || continue
        for iii in 1:nx, iiii in 1:nx
            if mask[iiii,iii] &&
               (x_arr[iiii]-island_centers[ii,1])^2 +
               (x_arr[iii]-island_centers[ii,2])^2 <= island_radius^2
                push!(primary_nodes[ii], (iiii,iii))
            end
        end
    end

    C  = zeros(nx, nx); Cn = zeros(nx, nx)
    gx = zeros(nx, nx); gy = zeros(nx, nx); S = zeros(nx, nx)
    E_st = zeros(n_cells); I_st = zeros(n_cells)

    Cauto = two_field ? zeros(nx, nx) : C

    Cauton = split_k ? zeros(nx, nx) : Cn
    gxa   = two_field ? zeros(nx, nx) : gx
    gya   = two_field ? zeros(nx, nx) : gy
    E_st2 = two_field ? zeros(n_cells) : E_st
    I_st2 = two_field ? zeros(n_cells) : I_st

    Sprim = zeros(nx, nx)
    for ii in 1:M
        island_primary_Q[ii] > 0 || continue
        dens = island_primary_Q[ii] / ISLAND_REF_AREA
        for iii in primary_nodes[ii]; Sprim[iii[1], iii[2]] += dens; end
    end
    Cprim = zeros(nx, nx); Cprimn = zeros(nx, nx)

    n_steps = round(Int, total_time / dt)
    snap_int = max(1, n_steps ÷ n_snapshots)
    σ_θ = sqrt(2 * rot_diffusion * dt)

    counts() = Int[count(==(ii), home) for ii in 1:M]
    snaps = [(t=0.0, pos=copy(pos), home=copy(home), n_rec=counts(),
              C_grid=copy(C), C_prim=copy(Cprim), C_auto=copy(split_k ? Cauto : C .- Cprim))]

    for ii in 1:n_steps

        for iii in 1:n_cells
            if home[iii] == 0
                for iiii in 1:M
                    if (pos[iii,1]-island_centers[iiii,1])^2 +
                       (pos[iii,2]-island_centers[iiii,2])^2 <= island_radius^2
                        home[iii] = iiii; break
                    end
                end
            end
        end

        fill!(S, 0.0)
        for iii in 1:n_cells
            home_island = home[iii]
            (home_island > 0 && island_Q[home_island] > 0) || continue
            in_island(pos[iii,1], pos[iii,2], island_centers, island_radius,
                      home_island) || continue
            deposit_bilinear!(S, x_arr, dx, pos[iii,1], pos[iii,2],
                              island_Q[home_island], nx, mask)
        end

        if split_k

            ar_diffstep!(Cauto, Cauton, mask, D, k_a,  S,     dx, dt, nx)
            ar_diffstep!(Cprim, Cprimn, mask, D, k,    Sprim, dx, dt, nx)
            C .= Cprim .+ Cauto
            ar_grad!(gx,  gy,  Cprim, mask, dx, nx)
            ar_grad!(gxa, gya, Cauto, mask, dx, nx)
        else

            S .+= Sprim

            ar_diffstep!(C, Cn, mask, D, k, S, dx, dt, nx)
            ar_diffstep!(Cprim, Cprimn, mask, D, k, Sprim, dx, dt, nx)
            if two_field
                Cauto .= C .- Cprim
                ar_grad!(gx,  gy,  Cprim, mask, dx, nx)
                ar_grad!(gxa, gya, Cauto, mask, dx, nx)
            else
                ar_grad!(gx, gy, C, mask, dx, nx)
            end
        end

        for iii in 1:n_cells
            px, py = pos[iii,1], pos[iii,2]

            local bwx, bwy
            if two_field
                bpx, bpy = sense_field(Cprim, gx, gy, x_arr, dx, nx, px, py,
                    sensing_radius, Kd, strategy, ck_sensitivity, E_st, I_st,
                    iii, dt, τ_E, τ_I)
                bax, bay = sense_field(Cauto, gxa, gya, x_arr, dx, nx, px, py,
                    sensing_radius, Kd_a, strategy, ck_sensitivity, E_st2, I_st2,
                    iii, dt, τ_E, τ_I)
                bwx = bias_weight * bpx + bias_auto * bax
                bwy = bias_weight * bpy + bias_auto * bay
            else
                bx, by = sense_field(C, gx, gy, x_arr, dx, nx, px, py,
                    sensing_radius, Kd, strategy, ck_sensitivity, E_st, I_st,
                    iii, dt, τ_E, τ_I)
                bwx = bias_weight * bx
                bwy = bias_weight * by
            end

            ra = dirs[iii] + σ_θ * randn()
            cvx = cos(ra) + bwx
            cvy = sin(ra) + bwy
            cm = sqrt(cvx^2 + cvy^2)
            if cm < 1e-10; cm = 1.0; end
            cvx /= cm; cvy /= cm

            npx = px + cell_speed * dt * cvx
            npy = py + cell_speed * dt * cvy
            nd = atan(cvy, cvx)

            rn = sqrt(npx^2 + npy^2)
            if rn > well_radius
                ov = rn - well_radius
                npx -= 2 * ov * npx / rn
                npy -= 2 * ov * npy / rn
                nd = π - nd + 2 * atan(npy, npx)
            end

            pos[iii,1] = npx; pos[iii,2] = npy; dirs[iii] = nd
        end

        if ii % snap_int == 0 || ii == n_steps
            push!(snaps, (t=ii*dt, pos=copy(pos), home=copy(home),
                          n_rec=counts(), C_grid=copy(C),
                          C_prim=copy(Cprim), C_auto=copy(split_k ? Cauto : C .- Cprim)))
        end
    end
    return snaps
end
