# Section_Modelling.jl — compute-only library
#
# Numerical results and plot-data for the spindle-torus and τ-parametrisation figures.
# All plotting, styling and notes live in paper.jl; this file only computes.

using LinearAlgebra: norm

"""
    compute_modelling()

Compute every numerical result and plot-data array needed by the modelling figures.

# Returned fields
- `a`: arm length (2 m).
- `l`: tether length (3 m).
- `β0`: example line elevation, deg2rad(30).
- `kite_pos`: example kite position [a + l*cos(β0), 0, l*sin(β0)].
- `xs`, `ys`, `zs`: mesh coordinates of the spindle torus as flat, NaN-separated line segments.
- `θc`: 200-sample angle grid for the arm-tip circle and the example lines.
- `τ_φ`: example τ for the Δφ marker (90°).
- `τ_θ`: example τ for the Δθ marker (225°).
- `dx`, `dy`: arrow components of the two τ̇ quivers.
- `up_col`, `down_col`: palette indices for the two quiver directions.
- `τ`: 501-sample angle grid for the (φ, θ) lemniscate curve.
- `θ`: coordinate function θ(τ) = sind(2τ).
- `φ`: coordinate function φ(τ) = sind(τ).
- `phi_vals`: φ evaluated on `τ` (x coordinates of the lemniscate curve).
- `theta_vals`: θ evaluated on `τ` (y coordinates of the lemniscate curve).
"""
function compute_modelling()
    # Spindle torus: set of possible kite positions
    #
    #   arm tip  A(α) = a * (cos α, sin α, 0)
    #   kite     K(α, β) = A(α) + l * (cos β * (cos α, sin α, 0) + sin β * ẑ)
    #
    # The line of length l points anywhere in the radial-vertical plane. Revolving
    # the sphere of radius l about the arm axis gives a torus; l > a => spindle.
    a, l = 2, 3

    # sanity check: a sampled point is exactly l from its arm tip
    β0 = deg2rad(30)  # line elevation
    kite_pos = [a + l * cos(β0), 0, l * sin(β0)]
    @assert isapprox(norm(kite_pos - [a, 0, 0]), l; atol=1e-8)

    # mesh as flat line segments (GR's wireframe chokes on matrix grids); a mesh
    # rather than a surface so the self-intersecting interior stays visible
    lg = range(0, 2π; length=120)  # points along each grid line
    lk = range(0, 2π; length=40)   # grid lines
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    for βj in lk
        for αi in lg
            push!(xs, (a + l * cos(βj)) * cos(αi))
            push!(ys, (a + l * cos(βj)) * sin(αi))
            push!(zs, l * sin(βj))
        end
        push!(xs, NaN)
        push!(ys, NaN)
        push!(zs, NaN)
    end
    for αj in lk
        for βi in lg
            push!(xs, (a + l * cos(βi)) * cos(αj))
            push!(ys, (a + l * cos(βi)) * sin(αj))
            push!(zs, l * sin(βi))
        end
        push!(xs, NaN)
        push!(ys, NaN)
        push!(zs, NaN)
    end

    # set of possible arm-tip positions: circle of radius a in the xy-plane
    θc = range(0, 2π; length=200)

    ## Parametrisation in tau
    # TODO [use the plot from applications/ecc2026]

    τ_φ = 90.0
    τ_θ = 225.0
    dx, dy = 1 / 4, 1 / 2
    up_col, down_col = 1, 2

    τ = range(0, 360; length=501)
    θ = τ -> sind(2 * τ)
    φ = τ -> sind(τ)
    phi_vals = φ.(τ)
    theta_vals = θ.(τ)

    # θ and φ are intentionally not returned: local closures cannot be
    # serialised across processes. paper.jl uses sind(2τ) / sind(τ) directly.
    return (;
        a,
        l,
        β0,
        kite_pos,
        xs,
        ys,
        zs,
        θc,
        τ_φ,
        τ_θ,
        dx,
        dy,
        up_col,
        down_col,
        τ,
        phi_vals,
        theta_vals,
    )
end
