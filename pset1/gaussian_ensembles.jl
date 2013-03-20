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
    npoints = 500
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

demo5() = for n = [5:5:25; 50:50:500]
        demo5_kde(n,nothing)
        pause()
    end
function demo5_hist(n,filename::Union(String,Nothing))
    nbins = sqrt(n)
    beta = 1
    e1 = sample_eigs(n,beta)
    e2 = sample_eigs(n,beta)
    e3 = sample_eigs(n,beta)
    e4 = sample_eigs(n,beta)
    step1 = (max(e1)-min(e1))/nbins
    step2 = (max(e2)-min(e2))/nbins
    step3 = (max(e3)-min(e3))/nbins
    step4 = (max(e4)-min(e4))/nbins
    b1 = min(e1) : step1 : max(e1)
    b2 = min(e2) : step2 : max(e2)
    b3 = min(e3) : step3 : max(e3)
    b4 = min(e4) : step4 : max(e4)
    b1 += step1/2
    b2 += step2/2
    b3 += step3/2
    b4 += step4/2
    a1 = hist(e1,b1);
    a2 = hist(e2,b2);
    a3 = hist(e3,b3);
    a4 = hist(e4,b4);
    x = [-2:0.01:2]
    p = plot(x,sqrt(4-x.^2)/(2*pi),"r-",
        b1, a1/sum(a1)/(b1[2]-b1[1]),"b-",
        b2, a2/sum(a2)/(b2[2]-b2[1]),"g-",
        b3, a3/sum(a3)/(b3[2]-b3[1]),"y-",
        b4, a4/sum(a4)/(b4[2]-b4[1]),"-")
    setattr(p, "title", "Demo5 (hist): n=$n")
    Winston.display(p)
    if isa(filename,String)
        file(p, "$(filename)_n=$n.png")
    end
end
function demo5_kde(n,filename::Union(String,Nothing))
    npoints = 2000
    beta = 1
    e1 = sample_eigs(n,beta)
    e2 = sample_eigs(n,beta)
    e3 = sample_eigs(n,beta)
    e4 = sample_eigs(n,beta)
    a1,b1 = kde(e1,npoints);
    a2,b2 = kde(e2,npoints);
    a3,b3 = kde(e3,npoints);
    a4,b4 = kde(e4,npoints);
    x = [-2:0.01:2]
    p = plot(x,sqrt(4-x.^2)/(2*pi),"r-",
        [b1], a1,"b-",
        [b2], a2,"g-",
        [b3], a3,"y-",
        [b4], a4,"-")
    setattr(p, "title", "Demo5 (kde): n=$n")
    Winston.display(p)
    if isa(filename,String)
        file(p, "$(filename)_n=$n.png")
    end
end

h(x) = (x = x.*(abs(x).<2)+2.*(abs(x).>=2); sqrt(4-x.^2)/(2*pi))
function demo5diff_kde(n,filename::Union(String,Nothing))
    npoints = 2000
    beta = 1
    e1 = sample_eigs(n,beta)
    e2 = sample_eigs(n,beta)
    e3 = sample_eigs(n,beta)
    e4 = sample_eigs(n,beta)
    a1,b1 = kde(e1,npoints);
    a2,b2 = kde(e2,npoints);
    a3,b3 = kde(e3,npoints);
    a4,b4 = kde(e4,npoints);
    x = [-2:0.01:2]
    p = plot([b1], a1-h([b1]),"b-",
        [b2], a2-h([b2]),"g-",
        [b3], a3-h([b3]),"y-",
        [b4], a4-h([b4]),"-")
    println("mean $(mean((a1-h([b1])).^2))")
    setattr(p, "title", "Demo5 error term (kde): n=$n")
    Winston.display(p)
    if isa(filename,String)
        file(p, "$(filename)_n=$n.png")
    end
end

function demo()
    demo1()
    pause()
    demo2()
    pause()
    demo5()
    pause()
end

end
