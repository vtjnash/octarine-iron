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

demo2() = demo2_hist(2,nothing)
function demo2_hist(n,filename::Union(String,Nothing))
    nbins = 50
    beta = 1
    k = 2
    t = 40000
    v = samples(n,beta,k,t)
    step = (max(v)-min(v))/nbins
    b = min(v) : step : max(v)
    a = hist(v,b);
    b += step/2 # center the bin
    xx=(0:.01:1)*max(b)
    j=n*(n-1)*beta/2+n
    x=xx*(n^2/2)
    p = plot( xx, ((n^2/2)*x.^(j/2-1).*exp(-x/2)/2^(j/2)/gamma(j/2)), "r-",
        b, a/sum(a)/(b[2]-b[1]), "b-")
    setattr(p, "title", "Demo2 (hist): n=$n, j=$j")
    Winston.display(p)
    if isa(filename,String)
        file(p, "$(filename)_n=$n.png")
    end
end
function demo2_kde(n,filename::Union(String,Nothing))
    npoints = 50
    beta = 1
    k = 2
    t = 40000
    v = samples(n,beta,k,t)
    a, b = kde(v,npoints);
    xx=(0:.01:1)*max(b)
    j=n*(n-1)*beta/2+n
    x=xx*(n^2/2)
    p = plot( xx, ((n^2/2)*x.^(j/2-1).*exp(-x/2)/2^(j/2)/gamma(j/2)), "r-",
        b, a, "b-")
    setattr(p, "title", "Demo2 (kde): n=$n, j=$j")
    Winston.display(p)
    if isa(filename,String)
        file(p, "$(filename)_n=$n.png")
    end
end


function demo()
    demo1()
    pause()
    demo2()
end

end
