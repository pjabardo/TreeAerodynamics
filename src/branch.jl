
DragType = Union{Float64,DragCoeff}


"""
    Branch2d(D, Cd, xc, yc)

Represents the main characteristics of a tree branch (or any other body actually)

"""

struct Branch2d{CDTYPE <: DragType}
    "Diameter of the branch"
    D::Float64
    "Drag coefficient of the branch"
    Cd::CDTYPE
    "X center of the branch"
    xc::Float64
    "Y center of the branch"
    yc::Float64
    #Branch2d(D, Cd, xc, yc) = new{CDTYPE}(D, Cd, xc, yc)
end


#Branch2d(D, Cd::Float64, xc, yc) = Branch2d{Float64}(D, Cd, xc, yc)
#Branch2d(D, Cd::T, xc, yc) where T<:AbstractDrag = Branch2d{T}(D, Cd, xc, yc)

"Diameter of the branch"
diameter(b::Branch2d) = b.D

"Default drag coefficient of the branch"
dragcoeff(b::Branch2d) = dragcoeff(b.Cd)
"Drag coefficient for a given wind speed"
dragcoeff(b::Branch2d, U) = dragcoeff(b.Cd, U)
"Coordinates of the branch"
center(b::Branch2d) = (b.xc, b.yc)
"X coordinate of the branch"
xcenter(b::Branch2d) = b.xc
"Y coordinate of the branch"
ycenter(b::Branch2d) = b.yc
