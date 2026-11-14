using StaticArrays
using Plots
using SplitApplyCombine

import KEEP.PointMass4 as PM4
import KEEP.PointMassPara as PMP

include("state_funcs.jl")

p = PMP.build_vbpara()

α0, τ0, dα0, dτ0 = 0, 1, 0, 0
tf = 100
sol_plus = PM4.integrate(SA[α0, 1, dα0, dτ0, 0], tf, p; save_everystep=true)
sol_minus = PM4.integrate(SA[α0, -1, dα0, dτ0, 0], tf, p; save_everystep=true)

tf = 100
t = range(tf - 5, tf; step=0.01)
# ...[Keep your simulation code up to `t = range(...)`] ...

## This is very cool: torus plot
xyz_torus = (t, sol, sym=false) -> begin
    α_ = α(t, sol)
    τ_ = τ(t, sol)
    if sym
        α_ = -α_
        τ_ = -τ_
    end
    x = (2 + cos(τ_)) * cos(α_)
    y = (2 + cos(τ_)) * sin(α_)
    z = 1.5 * sin(τ_)
    return SA[x, y, z]
end

lw = 4 # slightly thicker lines for better visibility
plot_sym = false

# 1. Increase surface resolution to hide the blocky rendering mesh
u_surf = v_surf = range(-π, π; length=150)
X = [(2 + cos(v)) * cos(u) for u in u_surf, v in v_surf]
Y = [(2 + cos(v)) * sin(u) for u in u_surf, v in v_surf]
Z = [1.5 * sin(v) for u in u_surf, v in v_surf]

# Pre-calculate trajectory paths for cleanliness
xyz_plus_coords = invert(xyz_torus.(t, Ref(sol_plus)))
xyz_minus_coords = invert(xyz_torus.(t, Ref(sol_minus)))

if plot_sym
    xyz_plus_sym_coords = invert(xyz_torus.(t, Ref(sol_plus), true))
    xyz_minus_sym_coords = invert(xyz_torus.(t, Ref(sol_minus), true))
end

# 2. Create a helper function to generate clean, floating subplots
function create_view(cam_angle)
    # framestyle=:none removes the bounding box, grid, and axis ticks completely
    p = surface(
        X,
        Y,
        Z;
        c=:bone,              # A smooth, neutral color gradient
        alpha=0.6,            # Surface transparency
        colorbar=:none,       # Hide colorbar
        framestyle=:none,     # Hide all axes/grids
        camera=cam_angle,
        legend=false,
    )

    # 3. Use vibrant, contrasting colors for the limit cycles
    plot!(p, xyz_plus_coords...; lw=lw, color=:dodgerblue)
    plot!(p, xyz_minus_coords...; lw=lw, color=:crimson)

    if plot_sym
        plot!(p, xyz_plus_sym_coords...; lw=lw - 1, color=:dodgerblue, alpha=0.5, ls=:dash)
        plot!(p, xyz_minus_sym_coords...; lw=lw - 1, color=:crimson, alpha=0.5, ls=:dash)
    end

    return p
end

# 4. Generate 4 distinct views (azimuth, elevation)
f1 = create_view((0, 90))   # Top View
f2 = create_view((0, 0))    # Front View
f3 = create_view((90, 0))   # Side View
f4 = create_view((45, 30))  # Perspective/Isometric View

# 5. Combine into a clean layout
fig = plot(
    f1,
    f2,
    f3,
    f4;
    layout=(2, 2),
    size=(1000, 1000),
    background_color=:white,
    title=["Top View" "Front View" "Side View" "Perspective View"],
    titlefont=10,
)

display(fig)

## Makie
import CairoMakie as cm

fig = cm.Figure(; size=(1000, 1000), backgroundcolor=:white)

# Extract points (Makie uses Point3f arrays)
pts_plus = [cm.Point3f(xyz_torus(ti, sol_plus)...) for ti in t]
pts_minus = [cm.Point3f(xyz_torus(ti, sol_minus)...) for ti in t]

cameras = [(0.0, π / 2), (0.0, 0.0), (π / 2, 0.0), (π / 4, π / 6)]

for (i, cam) in enumerate(cameras)
    row, col = divrem(i - 1, 2) .+ 1

    # Perspectiveness = 0 creates perfectly orthographic 2D views for the sides/top
    ax = cm.Axis3(
        fig[row, col];
        azimuth=cam[1],
        elevation=cam[2],
        aspect=:data,
        perspectiveness=(i == 4 ? 0.4 : 0.0),
    )

    # Renders the torus as a solid object with smooth lighting
    cm.surface!(ax, X, Y, Z; color=:gainsboro, alpha=0.1, shading=cm.MultiLightShading)

    cm.lines!(ax, pts_plus; color=:dodgerblue, linewidth=4)
    cm.lines!(ax, pts_minus; color=:crimson, linewidth=4)

    cm.hidedecorations!(ax)
    cm.hidespines!(ax)
end

display(fig)
