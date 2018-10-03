    using PyPlot
include("TreeAerodynamics.jl")

myrand(n, xmin=0.0, xmax=1.0) = xmin .+ rand(n) * (xmax-xmin)


nb = 10
D0 = 1.0
D = fill(D0, nb)
Cd = fill(0.35, nb)

xc = myrand(nb, 1, 10)
yc = myrand(nb, 1, 10)

nb = length(D)
η1 = [0.0, 0.1, 0.25, 0.40, 0.55, 0.7, 0.849, 1.0, 1.15, 1.30, 1.45, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5];
x1 = x = TA.geomseq(D0/2, 100*D0, 30, 1.12); 
nw = length(x1)


branches = [TA.Branch2d(D[i], Cd[i], xc[i], yc[i]) for i in 1:nb]
wm = [TA.WakeModel2d(TA.dragcoeff(b), b.D) for b in branches]

wakes = [TA.wake2d(w, x1, η1) for w in wm];

ux, uy, xw, yw, niter= TA.wakeinterference(wakes, branches, maxiter=2000, rlx=0.5, rlx2=0.1)
