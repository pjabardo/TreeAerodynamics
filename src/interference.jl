



function inducevel(w::Wake2d, xw, yw, Uoo, x, y)

    nx = w.nx
    ny = w.ny

    u = 0.0
    v = 0.0

    for i = 1:nx
        
        xwi = xw[i]
        ywi = yw[i]
        nrmx, nrmy = nrmdir(xw[i+1]-xwi, yw[i+1]-ywi)

        for k = 1:ny

            y0 = w.yv[k,i]
            x₁ = xwi + y0*nrmx
            x₂ = xwi - y0*nrmx
            y₁ = ywi + y0*nrmy
            y₂ = ywi - y0*nrmy

            σ = w.σ[k,i]
            Γ = Uoo * w.Γ[k,i]

            δx₁ = x - x₁
            δy₁ = y - y₁

            δx₂ = x - x₂
            δy₂ = y - y₂

            r₁ = hypot(δx₁, δy₁)
            r₂ = hypot(δx₂, δy₂)

            nx₁ = -δy₁/r₁; ny₁ = δx₁/r₁
            nx₂ = -δy₂/r₂; ny₂ = δx₂/r₂

            σ² = σ*σ
            uθ₁ =  Γ/2π * r₁/(r₁*r₁ + σ²)
            uθ₂ = -Γ/2π * r₂/(r₂*r₂ + σ²)

            dx1 = uθ₁*nx₁
            dx2 = uθ₂*nx₂
            dy1 = uθ₁*ny₁
            dy2 = uθ₂*ny₂
            u +=  dx1 + dx2
            v +=  dy1 + dy2
        end
    end
    
#    println(u, " ", v)        
    return u, v        
    

end



function branchvel!(wakes, branches, Ux, Uy, xw, yw, Uxo, Uyo, Uoo=1.0)
    nb = length(wakes)

    Uxo .= Uoo
    Uyo .= 0.0
    
    for i = 1:nb
        uui = hypot(Ux[i], Uy[i])
        for k = i+1:nb
            du, dv = inducevel(wakes[i], xw[i], yw[i], uui, branches[k].xc, branches[k].yc)
            Uxo[k] += du
            Uyo[k] += dv

            uuk = hypot(Ux[k], Uy[k])
            du, dv = inducevel(wakes[k], xw[k], yw[k], uuk, branches[i].xc, branches[i].yc)
            Uxo[i] += du
            Uyo[i] += dv
        end
    end
    
end
using PyPlot
function wakedispl!(wakes, branches, Uoo, Ux, Uy, xw, yw, xwo, ywo)

    nb = length(wakes)
    for i = 1:nb
        nw = length(xw[i])
        xwi = xw[i]
        ywi = yw[i]
        xwoi = xwo[i]
        ywoi = ywo[i]
        sx, sy = versor(Ux[i], Uy[i])
        xwoi[1] = branches[i].xc + sx * branches[i].D/2
        ywoi[1] = branches[i].yc + sy * branches[i].D/2
        
        for j = 1:nw-1
            xx = (xwi[j+1] + xwi[j]) * 0.5
            yy = (ywi[j+1] + ywi[j]) * 0.5
            dx = (xwi[j+1] - xwi[j]) 
            dy = (ywi[j+1] - ywi[j])
            rr = hypot(dx, dy)
            uu = Uoo
            vv = 0.0
            for k = 1:nb
                if k != i
                    uuk = hypot(Ux[k], Uy[k])
                    du, dv = inducevel(wakes[k], xw[k], yw[k], uuk, xx, yy)
                    uu += du
                    vv += dv
                end
            end

            sx, sy = versor(uu, vv)
            xwoi[j+1] = xwoi[j] + sx * rr
            ywoi[j+1] = ywoi[j] + sy * rr
            
        end
    end
    
        
end

#function inducevel(w::Wake2d, xw, yw, Uoo, x, y)

function wakeinterference(wakes, branches; maxiter=20, err=1e-6, rlx=0.2, rlx2=0.05)

    

    nb = length(wakes)

    xw = [wakes[i].xw .+ branches[i].xc for i in 1:nb]
    yw = [fill(branches[i].yc, length(wakes[i].xw)) for i in 1:nb]

    xwo = [wakes[i].xw .+ branches[i].xc for i in 1:nb]
    ywo = [fill(branches[i].yc, length(wakes[i].xw)) for i in 1:nb]

    Ux = ones(nb)
    Uy = zeros(nb)
    Uxo = zeros(nb)
    Uyo = zeros(nb)
    niter = 0
    for iter = 1:maxiter
        niter = iter
        # Compute the velocity at each branch
        branchvel!(wakes, branches, Ux, Uy, xw, yw, Uxo, Uyo)
        
        errmax = max(maximum(abs, Ux-Uxo), maximum(abs, Uy-Uyo))

        for i in 1:nb
            Ux[i] = Ux[i] + rlx * (Uxo[i] - Ux[i])
            Uy[i] = Uy[i] + rlx * (Uyo[i] - Uy[i])
        end
        wakedispl!(wakes, branches, 1.0, Ux, Uy, xw, yw, xwo, ywo)

        errmax2 = 0.0
        for i in 1:nb
            errmax2 = max(errmax2, abs(xw[i][8] - xwo[i][8]), abs(yw[i][8] - ywo[i][8]))
        end
        
        println(iter, " - ", errmax, " - ", errmax2)
        idx = 1:15
        for i in 1:nb
            xwi = xw[i]
            ywi = yw[i]
            xwoi= xwo[i]
            ywoi= ywo[i]
            #plot(branches[i].xc, branches[i].yc, "go")
            
            #plot(xwi[idx], ywi[idx], "r-")
            #plot(xwoi[idx], ywoi[idx], "b--")
            
            for k = 1:length(xw[i])
                xwi[k] = xwi[k] + rlx2 * (xwoi[k] - xwi[k])
                ywi[k] = ywi[k] + rlx2 * (ywoi[k] - ywi[k])
            end
        end
        #sleep(2)
        #cla()
        
        if errmax2 < err
            break
        end
        
        
    end

    return Ux, Uy, xw, yw, niter
    
end
