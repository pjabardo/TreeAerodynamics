



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


function branchvel(k, wakes, branches, Ux, Uy, xw, yw, uoofun)
    nb = length(wakes)

    

    D = branches[k].D
    xc = branches[k].xc
    yc = branches[k].yc

    ξ = 3*D/8
    w₁ = 0.111111111111
    w₂ = 0.888888888889/4
    
    Uxo, Uyo = uoofun(xc, yc)
    for i = 1:nb
        if i != k
            uui = hypot(Ux[i], Uy[i])
        
            du0, dv0 = inducevel(wakes[i], xw[i], yw[i], uui, xc, yc)
            du1, dv1 = inducevel(wakes[i], xw[i], yw[i], uui, xc+ξ, yc)
            du2, dv2 = inducevel(wakes[i], xw[i], yw[i], uui, xc-ξ, yc)
            du3, dv3 = inducevel(wakes[i], xw[i], yw[i], uui, xc, yc+ξ)
            du4, dv4 = inducevel(wakes[i], xw[i], yw[i], uui, xc, yc-ξ)
            Uxo += w₁*du0 + w₂*(du1 + du2 + du3 + du4)
            Uyo += w₁*dv0 + w₂*(dv1 + dv2 + dv3 + dv4)

        end
    end

    return Uxo, Uyo
end


function wakedispl!(i, jinit, branch, wakes, Ux, Uy, xw, yw, xwo, ywo, uoofun)

    
    nb = length(wakes)
     
    nw = length(xw[i])
    xwi = xw[i]
    ywi = yw[i]
    uui = hypot(Ux[i], Uy[i])

    if jinit==1
        sx, sy = versor(Ux[i], Uy[i])
        xwo[1] = branch.xc + sx * branch.D/2
        ywo[1] = branch.yc + sy * branch.D/2
    else
        for j in 1:jinit-1
            xwo[j] = xwi[j]
            ywo[j] = ywi[j]
        end
    end
        
    for j = jinit:nw-1
        xx = (xwi[j+1] + xwi[j]) * 0.5
        yy = (ywi[j+1] + ywi[j]) * 0.5
        dx0 = wakes[i].xw[j+1] - wakes[i].xw[j]
        
        uu, vv = uoofun(xx, yy)
        for k = 1:nb
            if k != i
                uuk = hypot(Ux[k], Uy[k])
                du, dv = inducevel(wakes[k], xw[k], yw[k], uuk, xx, yy)
                uu += du
                vv += dv
                end
            end

        
        sx, sy = versor(uu, vv)
        
        # Manter a distância entre pontos fixos
        #xwo[j+1] = xwo[j] + sx * rr
        #ywo[j+1] = ywo[j] + sy * rr
        # Compensar a differença entre velocidades
        uuinduced = hypot(uu, vv)
        rr = dx0 * (1.0 + (uuinduced-uui) / uui) # Esta equação não é exata!
        xwo[j+1] = xwo[j] + sx * rr 
        ywo[j+1] = ywo[j] + sy * rr
        
        
    end
end
    
        




function wakeinterference(wakes, branches, uoofun; maxiter=20, err=1e-6, rlx=0.02, Lerr=10.0)

    

    nb = length(wakes)

    xw = [wakes[i].xw .+ branches[i].xc for i in 1:nb]
    yw = [fill(branches[i].yc, length(wakes[i].xw)) for i in 1:nb]

    xwo = [wakes[i].xw .+ branches[i].xc for i in 1:nb]
    ywo = [fill(branches[i].yc, length(wakes[i].xw)) for i in 1:nb]

    Ux = zeros(nb)
    Uy = zeros(nb)

    for i = 1:nb
        Ux[i], Uy[i] = uoofun(branches[i].xc, branches[i].yc)
    end
    
    Uxo = zeros(nb)
    Uyo = zeros(nb)
    niter = 0
    jinit = 1
    for iter = 1:maxiter
        niter = iter
        # Compute the velocity at each branch
        errmax = 0.0
        for k = 1:nb
            ux, uy = branchvel(k, wakes, branches, Ux, Uy, xw, yw, uoofun)
            erru = max(maximum(abs, Ux[k]-ux), maximum(abs, Uy[k]-uy))
            if erru > errmax
                errmax = erru
            end
            
            Uxo[k] = ux
            Uyo[k] = uy
        end
        #return (Uxo, Uyo)
        
        #for i in 1:nb
        #    Ux[i] = Ux[i] + rlx * (Uxo[i] - Ux[i])
        #    Uy[i] = Uy[i] + rlx * (Uyo[i] - Uy[i])
        #end

        # Calculate the displacement of the wake:

        for i = 1:nb
            wakedispl!(i, jinit, branches[i], wakes, Ux, Uy, xw, yw, xwo[i], ywo[i], uoofun)
        end

        
        Ux .= Uxo
        Uy .= Uyo

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
                xwi[k] = xwi[k] + rlx * (xwoi[k] - xwi[k])
                ywi[k] = ywi[k] + rlx * (ywoi[k] - ywi[k])
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
