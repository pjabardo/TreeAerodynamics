

function geomseq(x1, x2, n, r=1.1)
    if r==1.0
        return collect(range(x1, stop=x2, length=n))
    end
    L = x2-x1
    ns = n-1
    Sn = (r^ns-1)/(r-1)
    dx = L/Sn
    x = zeros(n)
    x[1] = x1
    x[end] = x2

    for i = 2:ns
        x[i] = x[i-1]+dx
        dx *= r
    end
    return x
end


function nrmdir(δx, δy)
    r = hypot(δx, δy)
    return -δy/r, δx/r
end


function versor(x,y)
    r = hypot(x,y)
    return x/r, y/r
end
