module lab
import Base.@kwdef
using Plots
using DataFrames
# Assume SI units unless otherwise specified!

export g, Δt, Bob, Pendulum, simulate, extremes, periods

const g = 9.81
const Δt = 1 / 10000

@kwdef struct Bob
    mass::Number
    diameter::Number
    γ::Number
end

struct Pendulum
    pivot::@NamedTuple{x::Number, y::Number}
    radius::Number
    θ::Number
    bob::Bob
    Pendulum(radius::Number, θ::Number, bob::Bob) = new((x = 0, y = 0), radius, θ, bob)
end

function simulate(pendulum::Pendulum, duration::Number; resistance=true)
    (; radius, θ, bob) = pendulum
    (; mass, diameter, γ) = bob

    if resistance == false
        γ = 0
    end

    # change in acceleration in a vacuum
    ddθ(θ::Number) = -(sin(θ) * g) / radius
    drag(ω::Number) = γ * (radius) * (diameter ^ 2) * (ω ^ 2) / mass

    df = DataFrame(time=Float64[], θ=Float64[], ω=Float64[], α=Float64[])
    ω = 0
    for t in range(0, duration; step=Δt)
        α = ddθ(θ) - copysign(drag(ω), ω)
        push!(df, [t, θ, ω, α])
        ω += α * Δt
        θ += ω * Δt
    end
    df
end

function extremes(df::AbstractDataFrame)::AbstractDataFrame
    extremes = DataFrame(first(df))
    for i in 3:nrow(df)
        if sign(df[i-1, :ω]) ≠ sign(df[i, :ω])
            push!(extremes, df[i-1, :])
        end
    end
    extremes
end

function periods(extremes::AbstractDataFrame)::Vector{Number}
    acc = []
    for i in range(1, nrow(extremes) - 2; step=2)
        push!(acc, extremes[i + 2, :time] - extremes[i, :time])
    end
    acc
end

function main()
    b = Bob(0.01, .1, 0.25)
    p = Pendulum(1., deg2rad(50), b)
    df = simulate(p, 30)
    @time df = simulate(p, 30)
    # e = extremes(df)
    df
end

end
