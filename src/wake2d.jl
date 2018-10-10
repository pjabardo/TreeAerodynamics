
abstract type AbstractWake end
abstract type AbstractWake2d end

struct Wake2d
    nx::Int
    ny::Int
    η::Vector{Float64}
    xw::Vector{Float64}
    yv::Matrix{Float64}
    Γ::Matrix{Float64}
    σ::Matrix{Float64}
    A::Matrix{Float64}
    ω::Matrix{Float64}
end
Base.Broadcast.broadcastable(w::Wake2d) = Ref(w)

function wake2d(w::WakeModel2d, x₁, η₁)
    
    nx₁ = length(x₁)
    ny₁ = length(η₁)

    nx = nx₁ - 1
    ny = ny₁ - 1

    yv = zeros(ny, nx)
    Γ  = zeros(ny, nx)
    σ  = zeros(ny, nx)

    xw = copy(x₁)
    γ = log(2)
    η₀ = copy(η₁) 
    #for k = 1:ny
    #    η₀[k] = (η₁[k+1] + η₁[k]) / 2
    #end
    a = zeros(ny, nx)
    ww = zeros(ny, nx)
    for i in 1:nx
        xmean = (xw[i] + xw[i+1]) / 2
        δx = (xw[i+1] -  xw[i])
        xm = xmfun(w, xmean)
        u₀ = u0fun1(w, xm)
        L₀ = L0fun1(w, xm)
        for k = 1:ny
            η = (η₀[k+1] + η₀[k])/2 
            yv[k,i] = η * L₀
            δy = (η₀[k+1] - η₀[k]) * L₀
            f = profile(η)
            dudy = 2u₀*γ*η*f / L₀
            dvdx = -yv[k,i]/4 * w.θ * sqrt(w.β/w.α) * sqrt(2/(w.θ*w.β)) *
                (xmean - w.x₀)^(-2.5) * f * (2γ*η^2 - 3)
            #ω = dvdx - dudy
            ω = -dudy
            A = δx * δy
            a[k,i] = A
            ww[k,i] = ω
            Γ[k,i] = ω*A
            σ[k,i] = min(δx, δy) / 2
        end
    end

    return Wake2d(nx, ny, η₀, xw, yv, Γ, σ, a, ww)
end


        
    


function updatewake2d!(w::WakeModel2d, wk::Wake2d)

    nx = wk.nx
    ny = wk.ny

    η₀ = wk.η
    xw = wk.xw
    yv = wk.yv
    Γ  = wk.Γ
    σ  = wk.σ
    
    γ = log(2)

    for i in 1:nx
        xmean = (xw[i] + xw[i+1]) / 2
        δx = (xw[i+1] - xw[i])
        xm = xmfun(w, xmean)
        u₀ = u0fun1(w, xm)
        L₀ = L0fun1(w, xm)
        for k = 1:ny
            η = (η₀[k+1] + η₀[k])/2 
            y[k,i] = η * L₀
            δy = (η₀[k+1] - η₀[k]) * L₀
            f = profile(η)
            dudy = 2u₀*γ*η*f / L₀
            dvdx = -y[k,i]/4 * w.θ * sqrt(w.β/w.α) * sqrt(2/(w.θ*w.β)) *
                (xmean - w.x₀)^(-2.5) * f * (2γ*η^2 - 3)
            ω = dvdx - dudy
            A = δx * δy
            Γ[k,i] = ω*A
            σ[k,i] = min(δx, δy) / 12
        end
    end

end

