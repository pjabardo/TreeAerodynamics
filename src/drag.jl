

abstract type AbstractDrag end


"""
    DragCoeff(cd0, cdmin, a)

Models a variable drag coefficient. Uses an exponential to calculate
a drag that varies with wind velocity:

``
C_D = C_{D,0}⋅exp(-α⋅U) + C_{D,min}
``
"""
struct DragCoeff <: AbstractDrag
    cd0::Float64
    cdmin::Float64
    α::Float64
end
DragCoeff(cd0, cdmin, α) = new(cd0, cdmin, α)
Base.Broadcast.broadcastable(cd::DragCoeff) = Ref(cd)

"Calculate the drag coefficient"
dragcoeff(c::DragCoeff) = c.cd0+c.cdmin
"Calculate the drag coefficient at velocity U"
dragcoeff(c::DragCoeff, U) = c.cd0*exp(-c.α*U) + c.cdmin

dragcoeff(c::Number) = c
dragcoeff(c::Number, U) = c


struct DragPowerLaw <: AbstractDrag
    Cd1::Float64
    p::Float64
    Umin::Float64
end
DragPowerLaw(cd, p, umin=0.0) = DragPowerLaw(cd, p, umin)
Base.Broadcast.broadcastable(cd::DragPowerLaw) = Ref(cd)

function dragcoeff(c::DragPowerLaw, U=5.0)
    uu = U > c.Umin ? U : c.Umin
    return c.Cd1 * uu ^ c.p
end

    
