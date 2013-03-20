module convergence
    
using Winston
pause() = (print("Press <enter> to continue"); readline(STDIN); nothing)

H(n::Integer, x::FloatingPoint) = factorial(n)*sum((i)->(-1)^i*(2x)^(n-2i)/factorial(i)/factorial(n-2i),0:div(n,2))
pi_j(j::Integer, x::FloatingPoint) = H(j,x)/sqrt(sqrt(pi)*factorial(j)*2^j)
#Kn(n::Integer, x::FloatingPoint, y::FloatingPoint) = exp((-x^2-y^2)/2)*sum([pi_j(j,x)*pi_j(j,y) for j=0:(n-1)])
Kn(n::Integer, x::FloatingPoint, y::FloatingPoint) = exp((-x^2-y^2)/2)*sqrt(n-1)*(pi_j(n-1,x)*pi_j(n,y)-pi_j(n,x)*pi_j(n-1,y))/(y-x)

#sine_kernel(m, x, y) = 2*sqrt(pi)*m*sin(x-y)/(x-y)
sine_kernel(m, x, y) = sqrt(2*m-1)/sqrt(m)/pi*sin(2*sqrt(int(m))*(x-y))/(x-y)

demo2() = demo2(nothing)
function demo2(filename::Union(Nothing,String))
    x = 2.005
    ys = -10:.1:10
    ms = Integer[1, 2, 4, 6, 9, BigInt(12), BigInt(15), BigInt(18), BigInt(21)]
    for m in ms
        p = plot(ys, [sine_kernel(m, x, y) for y in ys], "b--",
             ys, [Kn(2m, x, y) for y in ys], "r-")
        l1,l2 = p.content1.components
        setattr(l1, "label", "Sine Kernel")
        setattr(l2, "label", "K_n")
        add(p, Legend(0.1,0.9, Any[l1,l2]))
        setattr(p, "title", "K_n vs. Sine Kernel at X=$x - M=$m")
        setattr(p, "xlabel", "Y")
        setattr(p, "ylabel", "Z")
        Winston.display(p)
        if isa(filename,String)
            file(p, "$(filename)_$(m).png")
        end
        pause()
    end
end

demo3() = demo3(nothing)
function demo3(filename::Union(Nothing,String))
    x = 0.005
    ys = -30:.1:10
    ns = Integer[1, 2, 3, 4, 5, 10, 15, 20, 40, 60]
    for n in ns
        p = plot(ys, [(airy(x)*airyprime(y)-airyprime(x)*airy(y))/(x-y) for y in ys], "b--",
             ys, [1/sqrt(2)/n^(1/6)*Kn(n>20?BigInt(n):n, sqrt(2n)+x/sqrt(2)/n^(1/6), sqrt(2n)+y/sqrt(2)/n^(1/6)) for y in ys], "r-")
        l1,l2 = p.content1.components
        setattr(l1, "label", "Theoretical Limit")
        setattr(l2, "label", "Airy Process")
        add(p, Legend(0.1,0.9, Any[l1,l2]))
        setattr(p, "title", "Airy Process Convergence at X=$x - N=$n")
        setattr(p, "xlabel", "Y")
        setattr(p, "ylabel", "Z")
        Winston.display(p)
        if isa(filename,String)
            file(p, "$(filename)_$(n).png")
        end
        pause()
    end

end

function demo()
    demo2()
    demo3()
end

end
