module WeirdDetector

using Interpolations
using DataFrames

export Point, periodogram, aliases

struct Point 
    t :: Float32
    modt :: Float32
    F :: Float32
    sigmaF :: Float32
    smoothedF :: Float32
end
Point(t::Real, F::Real, sigmaF::Real) = Point(Float32(t), NaN, Float32(F), Float32(sigmaF), NaN)
Point(t::Real, F::Real) = Point(t, F, 5f-5)

function fold!(period::Float32, data::Array{Point})
    for (i,p) in enumerate(data)
        data[i] = Point(p.t, p.t%period, p.F, p.sigmaF, p.smoothedF)
    end
    sort!(data, by=(x -> x.modt))
end

"helper function for smoothing"
function access(data::Vector{A}, i::Int) :: A where A
    const npoints = length(data)
    j = mod(i,npoints) == 0 ? npoints : mod(i,npoints)
    data[j] 
end

"helper function for smoothing"
function wrappedtime(data::Vector{Point}, i::Int, P::Float32) :: Float32
    const npoints = length(data)
    modt = access(data, i).modt
    if i < 1
        return modt - P
    elseif i > npoints
        return modt + P
    else 
        return modt
    end
end

function smooth2!(data::Vector{Point}, P::Float32, kw::Float32)
    kw = Int(round(kw * length(data)))
    s = map(collect(-kw:kw)) do i
        p = access(data, i)
        (p.F / p.sigmaF, 1/p.sigmaF)
    end
    sumFoverSigma = sum(first.(s))
    sumInverseSigma = sum(last.(s))
    lb = -kw
    ub = kw 
    for i = 1:length(data)
        point = access(data, lb)
        sumFoverSigma -= point.F / point.sigmaF
        sumInverseSigma -= 1e0 / point.sigmaF
        lb += 1

        ub += 1
        point = access(data, ub)
        sumFoverSigma += point.F / point.sigmaF
        sumInverseSigma += 1e0 / point.sigmaF

        data[i] = Point(data[i].t, data[i].modt, data[i].F, data[i].sigmaF, 
                        Float32(sumFoverSigma/sumInverseSigma), kw)
    end 
    kw
end

"smooth the data with a rolling mean"
function smooth!(data::Vector{Point}, P::Float32, kw::Float32)
    kw = kw * P
    sumFoverSigma = 0e0 
    sumInverseSigma = 0e0 
    lb = 2
    ub = 1 
    while wrappedtime(data, lb, P) > data[1].modt - kw
        lb -= 1
        point = access(data, lb)
        sumFoverSigma += point.F / point.sigmaF
        sumInverseSigma += 1e0 / point.sigmaF
    end
    for i = 1:length(data)
        while wrappedtime(data, lb, P) < data[i].modt - kw
            point = access(data, lb)
            sumFoverSigma -= point.F / point.sigmaF
            sumInverseSigma -= 1e0 / point.sigmaF
            lb += 1
        end 
        while wrappedtime(data, ub, P) < data[i].modt + kw
            ub += 1
            point = access(data, ub)
            sumFoverSigma += point.F / point.sigmaF
            sumInverseSigma += 1e0 / point.sigmaF
        end 
        data[i] = Point(data[i].t, data[i].modt, data[i].F, data[i].sigmaF, 
                        Float32(sumFoverSigma/sumInverseSigma), ub-lb+1)
    end 
    mean((x->x.kw).(data)) 
end

"calculate reduced chi squared on data that has already been folded and smoothed"
function chi2(data::Vector{Point}) :: Float32
    chi2 = 0e0
    for p in data
        chi2 += ((p.F - p.smoothedF)/p.sigmaF)^2
    end
    Float32(chi2)
end

"fold, smooth, and calculate reduced chi squared"
function chi2(data::Vector{Point}, period::Float32; kw::Float32=0.002f0)
    fold!(period, data)
    smooth!(data, period, kw)
    chi2(data)
end

function kurtosis(data::Vector{Point}) :: Float32
    fs = (x->x.smoothedF).(data)
    mu = mean(fs)
    sigma = std(fs)
    mean((fs .- mu).^4)/sigma^4
end

"returns dataframe containing chi-squared and kurtosis by period"
function periodogram(data::Vector{Point}, periods::Vector{Float32}, kw=0.002f0; 
                     parallel=true, datakw=false)
    df = DataFrame(chi2=Float32[], kurtosis=Float32[])
    const s = div(length(periods), nworkers()*3)
    results = pmap(periods, distributed=parallel, batch_size=s) do p
        fold!(p, data)
        row = Vector{Any}()
        if datakw
            smooth2!(data, p, kw)
        else
            smooth!(data, p, kw)
        end
        push!(row, chi2(data))
        push!(row, kurtosis(data))
        push!(df, row)
    end
    df
end

function aliases(downto::Int=5; upperBound=nothing) :: Vector{Rational}
    lines = Set{Rational}()
    for n in 1:downto
        for m in 1:n
            if upperBound == nothing
                push!(lines, m//n)
            else
                i = 0
                while i + m//n < upperBound 
                    push!(lines, i + m//n)
                    i += 1
                end
            end
        end
    end
    collect(lines)
end

function flatten(periods::Vector{Float32}, chi2s::Vector{Float32}, stepwidth::Float32=10.21f0;
                 preflipped=false)
    nchi2s = similar(chi2s)
    steps = Int(ceil(periods[end]/stepwidth))
    for n in 1:steps
        lb = first(searchsorted(periods, (n-1)*stepwidth))
        ub = last(searchsorted(periods, n*stepwidth))
        line(x, p) = p[1]*x + p[2]
        c2(p) = sum((line(periods[lb:ub], p) .- chi2s[lb:ub]).^2)
        p0 = [0e0, 0e0]
        p0[1] = (chi2s[ub] - chi2s[lb])/(periods[ub] - periods[lb])
        p0[2] = chi2s[lb] - p0[1]*periods[lb]
        res = Optim.optimize(c2, p0, Newton())
        p = Optim.minimizer(res)
        nchi2s[lb:ub] .= Float32.(line(periods[lb:ub], p) - chi2s[lb:ub])
    end
    if preflipped
        -nchi2s
    else
        nchi2s
    end
end

function movingstd(npowers::Vector{Float32}, kw=200) :: Vector{Float32}
    sigmas = similar(npowers)
    for i in 1:length(npowers)
        lb = i - kw < 1 ? 1 : i - kw
        ub = i + kw > length(npowers) ? length(npowers) : i + kw
        sigmas[i] = std(npowers[lb:ub])
    end
    sigmas
end

function scrambled_periodogram(df::DataFrame, periods::Vector{Float32}; kwargs...)
    df = deepcopy(df)
    dropmissing!(df)
    df[:F] = df[randperm(size(df)[1]), :F]
    df[:sigmaF] = df[randperm(size(df)[1]), :sigmaF]
    data = pointsify(df)
    periodogram(data, periods; kwargs...)
end

function optimal_periods(pmin=0.25f0, pmax=5f1; n=5)
    pmin = Float32(pmin)
    pmax = Float32(pmax)
   exp.((log(pmin) : (0.001f0/n) : log(pmax))) 
end

end
