#!../julia/julia

module GOE_Finite_Level_Densities
using Winston
#using Calculus
export demoGOE, demoGUE

### To be added to Base:
import Base: float, promote_rule, ^
float(x::BigInt) = BigFloat(x)
promote_rule{T<:FloatingPoint}(::Type{BigInt},::Type{T}) = BigFloat
^(x::Float32, y::BigInt) = BigFloat(x)^y
^(x::Float64, y::BigInt) = BigFloat(x)^y
###

HermiteH(n::Integer, x::FloatingPoint) = factorial(n)*sum((i)->(-1)^i*(2x)^(n-2i)/factorial(i)/factorial(n-2i),0:div(n,2))
HermiteH_D(n::Integer, x::FloatingPoint) = (n!=0 ? 2n*HermiteH(n-1, x) : 0)
HermiteH_I(n::Integer, x::FloatingPoint) = HermiteH(n-1, 0.0) - exp(-x^2)*HermiteH(n-1, x)

phi(n::Integer, x::FloatingPoint) = 1/sqrt(2^n*factorial(n)*sqrt(pi))*exp(-x^2/2)*HermiteH(n, x)
phi_D(n::Integer, x::FloatingPoint) = 1/sqrt(2^n*factorial(n)*sqrt(pi))*(-x*exp(-x^2/2)*HermiteH(n, x) + exp(-x^2/2)*HermiteH_D(n, x))
phi_I(n::Integer, x::FloatingPoint) = (delta=x*.001; (x == 0 ? 0 : sum((x)->delta*phi(n,x),0:delta:x))) # Euler numerical integration
#phi_I(n::Integer, x::FloatingPoint) = integrate((x)->phi(n,float(x)),0.0,x) # numerical integration
#phi_I(n::Integer, x::FloatingPoint) = 1/sqrt(2^n*factorial(n)*sqrt(pi))* ???

eigDensityGUE(n::Integer, x::FloatingPoint) = sum((k)->phi(k,x)^2, 0:(n-1))
#eigDensityGUE(n::Integer, x::FloatingPoint) = n*phi(n,x)^2-sqrt(n*(n+1))*phi(n-1,x)*phi(n+1,x)
plot_eigDensityGUE(n::Integer) = (xs=-6:.1:6; plot(xs, [eigDensityGUE(n, float(x)) for x=xs])) #note: N can be a regular Int, or a BigInt!

eigDensityGOE(n::Integer, x::FloatingPoint) = sum((k)->phi(2k,x)^2-phi_D(2k,x)*phi_I(2k,x), 0:(div(n,2)-1))
plot_eigDensityGOE(n::Integer) = (xs=-6:.1:6; plot(xs, [eigDensityGOE(n, float(x)) for x=xs]))

pause() = (print("Press <enter> to continue"); readline(STDIN); nothing)

demoGUE() = demoGUE(nothing, Any[1,4,16,BigInt(24),BigInt(32),BigInt(48),BigInt(64)])
function demoGUE(filename::Union(String,Nothing), range)
    println("Type Ctrl-C to cancel at any time")
    for i = range
        println("plot_eigDensityGUE($i)")
        p=plot_eigDensityGUE(i)
        setattr(p, "title", "")
        setattr(p, "title", "GUE Eig Density - N=$i")
        Winston.display(p)
        if isa(filename,String)
            file(p, "$(filename)_$(i).png")
        else
            pause()
        end
    end
end

demoGOE() = demoGOE(nothing, 1:20)
function demoGOE(filename::Union(String,Nothing), range)
    println("Type Ctrl-C to cancel at any time")
    for i = range*2
        local p
        println("plot_eigDensityGOE($i)")
        try
            p=plot_eigDensityGOE(i)
        catch e
            if isa(e, DomainError)
                p=plot_eigDensityGOE(BigInt(i))
            else
                rethrow()
            end
        end
        setattr(p, "title", "GOE Eig Density - N=$i")
        Winston.display(p)
        if isa(filename,String)
            file(p, "$(filename)_$(i).png")
        else
            pause()
        end
    end
end
    
end

