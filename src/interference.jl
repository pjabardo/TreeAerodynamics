



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

function meanvelcontrib(branch, wake, uu, xw, yw)

    D = branch.D
    xc = branch.xc
    yc = branch.yc

    ξ = 3*D/8
    w₁ = 0.111111111111
    w₂ = 0.222222222222

    du0, dv0 = inducevel(wake, xw, yw, uu, xc, yc)
    du1, dv1 = inducevel(wake, xw, yw, uu, xc+ξ, yc)
    du2, dv2 = inducevel(wake, xw, yw, uu, xc-ξ, yc)
    du3, dv3 = inducevel(wake, xw, yw, uu, xc, yc+ξ)
    du4, dv4 = inducevel(wake, xw, yw, uu, xc, yc-ξ)

    du = w₁*du0 + w₂*(du1 + du2 + du3 + du4)
    dv = w₁*dv0 + w₂*(dv1 + dv2 + dv3 + dv4)

    return du, dv

end


function branchvel(k, wakes, branches, Ux, Uy, xw, yw, uoofun)
    nb = length(wakes)

    Uxo, Uyo = uoofun(branches[k].xc, branches[k].yc)
    for i = 1:nb
        if i != k && Ux[i] > 0  # Does not deal with recirculation...
            
            uui = hypot(Ux[i], Uy[i])

            du, dv = meanvelcontrib(branches[k], wakes[i], uui, xw[i], yw[i])
            Uxo += du
            Uyo += dv

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
        xwo[1] = branch.xc + sx * wakes[i].xw[1]
        ywo[1] = branch.yc + sy * wakes[i].xw[1]
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
            #if k != i
                uuk = hypot(Ux[k], Uy[k])
                du, dv = inducevel(wakes[k], xw[k], yw[k], uuk, xx, yy)
                uu += du
                vv += dv
            #    end
        end

        
        sx, sy = versor(uu, vv)
        
        # Manter a distância entre pontos fixos
        #xwo[j+1] = xwo[j] + sx * rr
        #ywo[j+1] = ywo[j] + sy * rr
        # Compensar a differença entre velocidades
        uuinduced = hypot(uu, vv)
        #rr = dx0 * (1.0 + (uuinduced-uui) / uui) # Esta equação não é exata!
        rr = dx0
        if j >= nw-1
            xwo[j+1] = xwo[j] +  rr 
            ywo[j+1] = ywo[j] 
        else
            xwo[j+1] = xwo[j] + sx * rr 
            ywo[j+1] = ywo[j] + sy * rr
        end
        
        
    end
end
    
        

function lengthidx(xw, L)

    n = length(xw)
    idx = 1
    for i = 1:n
        idx = i
        if xw[i] > L
            break
        end
    end

    return idx
end


function wakeinterference!(nstart, nwakes, wakes, branches, xw, yw, xwo, ywo, Ux, Uy, uoofun; maxiter=30000, err=1e-3, rlx=0.002, rlxu=0.1, Lerr=10.0)

    
    

    nb = nwakes

    wakes1 = wakes[1:nb]
                   
    maxxw = maximum((w.xw[end] for w in wakes1))
    idxerr = [lengthidx(w.xw, Lerr) for w in wakes1]
    
    for i = 1:nb
        Ux[i], Uy[i] = uoofun(branches[i].xc, branches[i].yc)
    end
    
    niter = 0
    jinit = ones(Int, nb)
    
    for iter = 1:maxiter
        niter = iter
        # Compute the velocity at each branch
        println("Chegou")
        errmax = 0.0
        for k = 1:nb
            ux, uy = branchvel(k, wakes1, branches, Ux, Uy, xw, yw, uoofun)
            erru = max(maximum(abs, Ux[k]-ux), maximum(abs, Uy[k]-uy))
            if erru > errmax
                errmax = erru
            end
            
            Ux[k] = Ux[k] + rlxu * (ux - Ux[k])
            Uy[k] = Uy[k] + rlxu * (uy - Uy[k])
        end
        #return (Uxo, Uyo)
        
        #for i in 1:nb
        #    Ux[i] = Ux[i] + rlx * (Uxo[i] - Ux[i])
        #    Uy[i] = Uy[i] + rlx * (Uyo[i] - Uy[i])
        #end

        # Calculate the displacement of the wake:

        for i = nstart:nb
            wakedispl!(i, jinit[i], branches[i], wakes1, Ux, Uy, xw, yw, xwo[i], ywo[i], uoofun)
        end


        errmax2 = 0.0
        for i in 1:nb
            errmax2 = max(errmax2, abs(xw[i][idxerr[i]] - xwo[i][idxerr[i]]),
                          abs(yw[i][idxerr[i]] - ywo[i][idxerr[i]]))
        end

        if mod(iter, 200)==0
            println(iter, " - ", errmax, " - ", errmax2, " - ", idxerr[1])
        end 
            
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
        
        if errmax2 < err * Lerr
            if Lerr > maxxw
                break
            end
            Lerr = 1.5*Lerr
            jinit .= idxerr
            idxerr .= [lengthidx(w.xw, Lerr) for w in wakes1]
            println("==================================")
            println(Lerr, " - ", jinit[1], " - ", idxerr[1])
            
        end
        
        
    end

    return niter
    
end


function wakeinterference(wakes, branches, uoofun; maxiter=30000, err=1e-6, rlx=0.2, rlxu=0.3, Lerr=10.0)

    
    

    nb = length(wakes)

    xw = [wakes[i].xw .+ branches[i].xc for i in 1:nb]
    yw = [fill(branches[i].yc, length(wakes[i].xw)) for i in 1:nb]

    xwo = [wakes[i].xw .+ branches[i].xc for i in 1:nb]
    ywo = [fill(branches[i].yc, length(wakes[i].xw)) for i in 1:nb]

    Ux = zeros(nb)
    Uy = zeros(nb)
    


    for n = 2:nb
        
        wakeinterference!(n, wakes, branches, xw, yw, xwo, ywo, Ux, Uy, uoofun;
                          maxiter=maxiter, err=err, rlx=rlx, rlxu=0.1, Lerr=10.0)
        rlx = rlx * 0.8
    end
    
    return Ux, Uy, xw, yw
    
end


function fixedwakeinterference(wakes, branches, uoofun; maxiter=30000, err=1e-5, rlx=0.2)

    nb = length(wakes)

    xw = [wakes[i].xw .+ branches[i].xc for i in 1:nb]
    yw = [fill(branches[i].yc, length(wakes[i].xw)) for i in 1:nb]


    Ux = zeros(nb)
    Uy = zeros(nb)
    Uxo = zeros(nb)
    Uyo = zeros(nb)

    for i = 1:nb
        Ux[i], Uy[i] = uoofun(branches[i].xc, branches[i].yc)
    end
    
    niter = fixedwakeinterference!(wakes, branches, uoofun, Ux, Uy, Uxo, Uyo, xw, yw;
                                  maxiter=maxiter, err=err, rlx=rlx)
    
   
    return Ux, Uy, niter
    
end

function fixedwakeinterference!(wakes, branches, uoofun, Ux, Uy, Uxo, Uyo, xw, yw;
                                maxiter=30000, err=1e-5, rlx=0.2)

    nb = length(wakes)

    uxmax = maximum(abs, Ux)
    uymax = maximum(abs, Uy)

    umax = hypot(uxmax, uymax)
    
    niter = 0
    
    for iter = 1:maxiter
        
        niter = iter
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
        if mod(iter,10)==0
            println(iter, "; ", errmax)
        end
        
        @. Ux = Ux + rlx * (Uxo - Ux)
        @. Uy = Uy + rlx * (Uyo - Uy)
        
        if errmax < err*umax
            println("Converged: ", iter, "; ", errmax/umax)
            break
        end
    end
    
    return niter
    
end








function velinterference(branches::Vector{Branch2d{T}}, uoofun, x, η;
                          maxiter=30000, err=1e-5, rlx=0.2,
                          Uxinit=nothing, Uyinit=nothing) where {T <: DragType}

    nb = length(branches)
    nw = length(x)
    xw = [x .+ branches[i].xc for i in 1:nb]
    yw = [fill(branches[i].yc, nw) for i in 1:nb]
    D = [diameter(b) for b in branches]

    Cd = zeros(nb)
    Cd2 = zeros(nb)
    
    Ux = zeros(nb)
    Uy = zeros(nb)

    Um = zeros(nb)
    Um2 = zeros(nb)

    for i = 1:nb
        Ux[i], Uy[i] = uoofun(branches[i].xc, branches[i].yc)
    end
    
    if Uxinit != nothing
        Ux .= Uxinit
    end
    if Uyinit != nothing
        Uy .= Uyinit
    end

    for i = 1:nb
        Um[i] = hypot(Ux[i], Uy[i])
    end
        

    
    Uxo = zeros(nb)
    Uyo = zeros(nb)

    wm = Array{WakeModel2d,1}(undef, nb)
    wakes = Array{Wake2d,1}(undef, nb)
    maxiter1 = maxiter

    if T==Float64
        maxiter1 = 1
    end

    for i = 1:nb
        Cd2[i] = dragcoeff(branches[i], Um[i])
    end
    
    lastiter = false
    for iter = 1:maxiter1
        println("============================================")
        println("CD iter: ", iter)
        
        for i = 1:nb
            
            Cd[i] = Cd2[i]
            wm[i] = WakeModel2d(Cd[i], D[i])
            wakes[i] = wake2d(wm[i], x, η)

        end

        err1 = lastiter ? err : err*100
        niter = fixedwakeinterference!(wakes, branches, uoofun, Ux, Uy, Uxo, Uyo, xw, yw;
                                       maxiter=maxiter, err=err1, rlx=rlx)
        @. Um2 = hypot(Ux, Uy)
        @. Cd2 = dragcoeff(branches, Um2)
         
        #display(hcat(Um, Cd, Um2, Cd2))
        errmax = maximum(abs, (Cd2[i] - Cd[i] for i = 1:nb))
        if lastiter
            break
        elseif errmax < err
            lastiter = true
        end
        
        println("CD Iteration: ", iter, "; ", errmax)
    end   
    #@. Cd = dragcoeff(branches, Um)
    
    return Ux, Uy, Cd
    
end

#function maketree(L, D, Cd

function simulatevels(branches, Uoo, x, η; err=1e-5, rlx=0.2)
    
    nu = length(Uoo)
    nb = length(branches)
    drag = []
    
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    println("Vel: ", Uoo[1])
    push!(drag, velinterference(branches, (x,y)->(Uoo[1], 0.0), x, η; err=err, rlx=rlx))
    
    for i = 2:nu
    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    println("Vel: ", Uoo[i])
        fun = (x,y) -> (Uoo[i], 0.0)
        Ux = drag[i-1][1] .* Uoo[i]/Uoo[i-1]
        Uy = drag[i-1][1] .* Uoo[i]/Uoo[i-1]
        rlx = min(0.5, rlx*1.5)
        push!(drag, velinterference(branches, fun, x, η, err=err, rlx=rlx, Uxinit=Ux, Uyinit=Uy))
    end
    
    return drag
    
    
end
