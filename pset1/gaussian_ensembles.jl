module Gaussian_Ensembles

using Winston
using Distributions

function sample_eigs(n,beta)
  d = randn(n)*sqrt(2)
  s = float64([sqrt(rand(Chisq(beta*(n-i)))) for i=1:(n-1)])
  e = eigvals(SymTridiagonal(d,s))/sqrt(n)
end

function sample_kth_moment(n,beta,k)
    mean(sample_eigs(n,beta).^k)
end


function samples(n,beta,k,t)
    v = zeros(t);
    for i = 1:t
        v[i] = sample_kth_moment(n,beta,k)
    end
    v
end

pause() = (print("Press <enter> to continue"); readline(STDIN); nothing)

demo1() = demo1(1,nothing)
function demo1(k,filename::Union(String,Nothing))
    println()
    n = 10
    beta = 1
    v = 0.0
    step = 1000
    xmax = 12000
    runs = 40
    xs = step:step:xmax
    ys = zeros(div(xmax,step),runs)
    for run = 1:runs
        for i = xs
            ys[i/step,run] = mean(samples(n,beta,k,i))
        end
        print("run $run or $runs\r")
    end
    println()
    xs = [xs]*ones(1,runs)
    p = plot(xs,ys,"*")
    setattr(p, "title", "Demo1: mean k=$k")
    Winston.display(p)
    if isa(filename,String)
        file(p, "$(filename).png")
    end
end

function demo()
    demo1()
    pause()

end

end
