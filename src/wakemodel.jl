

abstract type AbstractWakeModel end
abstract type AbstractWakeModel2d <: AbstractWakeModel end

"""

"""
struct WakeModel2d <: AbstractWakeModel2d
    Cd::Float64
    D::Float64
    θ::Float64
    ηlim::Float64
    α::Float64
    β::Float64
    x₀::Float64
    function WakeModel2d(Cd, D, ηlim=2.0, α=0.7231, β=0.1859)
        θ = Cd * D / 2
        x₀ = D*(0.5 - 1.0 / (β*Cd*ηlim^2))
        new(Cd, D, θ, ηlim, α, β, x₀)
    end
end
Base.Broadcast.broadcastable(w::WakeModel2d) = Ref(w)


xmfun(w::WakeModel2d, x) = (x - w.x₀) / (2*w.θ)
u0fun1(w::WakeModel2d, xm) = 1.0 / sqrt(w.α * xm)
u0fun(w::WakeModel2d, x) = u0fun1(w, xmfun(w, x))
L0fun1(w::WakeModel2d, xm) = w.θ * sqrt(w.β * xm)
L0fun(w::WakeModel2d, x) = L0fun1(w, xmfun(w, x))

profile(η) = exp(-log(2)*η^2)

function ufun(wake::WakeModel2d, x, y, Uoo=1.0)
    xm = xmfun(wake, x)
    L0 = L0fun1(wake, xm)
    u0 = u0fun1(wake, xm)
    Uoo * (1.0 - u0*profile(y/L0))
end

function vfun(wake::WakeModel2d, x, y, Uoo=1.0)
    xm = xmfun(wake, x)
    L0 = L0fun1(wake, xm)
    u0 = u0fun1(wake, xm)
    γ = log(2)
    η = y/L0
    return Uoo * (-u0*u0/4*sqrt(wake.α*wake.β) * η * profile(η))
end


function vorticity(wake::WakeModel2d, x, y, Uoo=1.0)

    xm = xmfun(wake, x)
    L0 = L0fun1(wake, xm)
    u0 = u0fun1(wake, xm)
    η = y/L0
    f = profile(η)
    γ = log(2)
    dudy = 2u0*γ*η*f / L0
    dvdx = -y/4 * wake.θ * sqrt(wake.β/wake.α) * sqrt(2/(wake.θ*wake.β)) *
        (x - wake.x₀)^(-2.5) * f * (2γ*η^2 - 3)

    return dvdx - dudy
end

