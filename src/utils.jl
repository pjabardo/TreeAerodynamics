
"""
    geomseq(x1, x2, n[, r=1.1])

Generate a sequence of `n` numbers between `x1` and `x2`
where the ratio `(x[i+1]-x[i]) / (x[i]-x[i-1]) == r`

The sequence of numbers include both ends `x1` and `x2`.

# Examples
```julia-repl
julia> geomseq(1,2,10, 1.2)
10-element Array{Float64,1}:
 1.0               
 1.0480794616724993
 1.1057748156794986
 1.1750092404878978
 1.2580905502579767
 1.3577881219820713
 1.477425208050985 
 1.6209897113336815
 1.7932671152729172
 2.0               
```
"""
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

"""
    nrmdir(δx, δy)

Generate an unit vector normal to the vector (δx, δy) moving counter-clockwise

# Examples
```julia-repl
julia> nrmdir(1,0)
(0.0, 1.0)

julia> nrmdir(2,0)
(0.0, 1.0)

julia> nrmdir(0.5, 0.5)
(-0.7071067811865475, 0.7071067811865475)

```
"""
function nrmdir(δx, δy)
    r = hypot(δx, δy)
    return -δy/r, δx/r
end


"""
    versor(x,y)

Return a unit vector with the same direction as the vector (x,y)

# Examples
```julia-repl
julia> versor(5,0)
(1.0, 0.0)

julia> versor(1,1)
(0.7071067811865475, 0.7071067811865475)
```

"""
function versor(x,y)
    r = hypot(x,y)
    return x/r, y/r
end
