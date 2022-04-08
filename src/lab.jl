# Assume SI units unless otherwise specified!

const g = 9.81
const Δt = 1 / 1000

struct Bob
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

function simulate(pendulum::Pendulum, duration::Number)
    (; radius, θ, bob) = pendulum
    (; mass, diameter, γ) = bob

    # change in acceleration in a vacuum
    ddθ(θ::Number) = -(sin(θ) * g) / radius
    drag(ω::Number) = γ * (radius) * (diameter ^ 2) * (ω ^ 2) / mass

    ts = range(0, duration; step=Δt)
    pos::AbstractVector{Number} = []
    ω = 0
    for t in ts:
        # println("t = $t:\t(θ = $(rad2deg(θ))°;\tω = $ω; α = $α)")
        push!(pos, rad2deg(θ))
        α = ddθ(θ) - copysign(drag(ω), ω)
        ω += α * Δt
        θ += ω * Δt
    end
    println("Finish")
    ts, pos
end

function main()
    b = Bob(0.01, .1, .25)
    p = Pendulum(1., deg2rad(15), b)
    ts, θs = simulate(p, 30)
end

clearconsole()
main()
