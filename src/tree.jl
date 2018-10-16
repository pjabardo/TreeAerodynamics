
struct HexagonTree{T <: DragType}
    xc::Float64
    yc::Float64
    D::Float64
    Cd::T
    level::Int
    branches::Array{Branch2d{T},1}
end



function HexagonTree(branch::Branch2d{T}) where {T <: DragType}

    D = diameter(branch)
    Cd = branch.Cd
    xc = xcenter(branch)
    yc = ycenter(branch)

    level = 1
    
    b = Array{Branch2d{T},1}(undef, 7)

    b[1] = branch
    θ = 0.0
    r = D
    for i = 1:6
        b[i+1] = Branch2d(D, Cd, xc + r * cosd(θ), yc + r * sind(θ))
        θ += 60.0
    end

    return HexagonTree{T}(xc, yc, D, Cd, level, b)

end

function HexagonTree(h::HexagonTree{T}) where {T <: DragType}
    D = h.D
    xc = h.xc
    yc = h.yc
    Cd = h.Cd
    level = h.level
    branches = h.branches
    nb1 = length(branches)

    b = Array{Branch2d{T},1}(undef, nb1*7)

    
    r = D * 3^level
    
    for k = 1:nb1
        b[k] = Branch2d{T}(D, Cd, branches[k].xc, branches[k].yc)
    end
    count = nb1 + 1
    θ = 0.0
    for i = 1:6
        dx = xc + r * cosd(θ)
        dy = yc + r * sind(θ)
        for k = 1:nb1
            b[count] = Branch2d{T}(D, Cd, branches[k].xc + dx, branches[k].yc + dy)
            count += 1
        end
        θ += 60.0
    end

    return HexagonTree{T}(xc, yc, D, Cd, level+1, b)

end

diameter(h::HexagonTree) = h.D * 3^(h.level)

function force(branches, Ux, Uy, ρ=1.2)

    Um = hypot.(Ux, Uy)

    Cd = dragcoeff.(branches, Um)
    α = atan(Uy, Ux)
    cα = cos.(α)
    sα = sin.(α)
    Fx = 0.0
    Fy = 0.0

    nb = length(Um)

    for i = 1:nb
        D = diameter(branches)
        fi = 0.5 * ρ * Cd * D * Um[i]^2
        
        Fx += fi*cα[i]
        Fy += fi*sα[i]
    end

    return Fx, Fy
end

function forcecoeff(branches, Ux, Uy, Uoo, Lref)


end


        
    
