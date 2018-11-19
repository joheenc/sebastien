This code implements the algorithm in [Wheeler & Kipping 201?]().  Find posterior samples from the paper at [github.com/ajwheeler/WDdata](https://github.com/ajwheeler/WDdata).

This code is tested with Julia 1.0, but it should be compatible with 0.6 as well.  It requires the `Interpolations`, `DataFrames`, and `FITSIO` packages.  To install them run these commands on the Julia REPL

    Pkg.add("Interpolations")
    Pkg.add("DataFrames")
    Pkg.add("FITSIO")

Here is a brief usage example

    push!(LOAD_PATH, "/directory/containing/this/code")
    using WeirdDetector
    
    #here's how you can get a Kepler light curve detrended like we did in the paper
    getFITS(8462852, fitsdir="./") #download FITS files for Boyajian's star to the current directory
    df = loadFITS(8462852, fitsdir="./") #load all quarters into single data frame, perform outlier rejection and detrending
    data = pointsify(df) #convert to data type taken by periodogram()

    #Here's how to prep arbitary data
    t = ...      # array of epochs
    F = ...      # array of flux values
    sigmaF = ... # array of 1-sigma uncertainties
    data = Point.(t, F, sigmaF) #this constructs an array of Points with Julia's broadcast syntax

    #construct array of periods used in paper.  Any array of Float32's can be used instead.
    periods = optimal_periods(2.0, 50.0)
    #construct a DataFrame containing the chi-squared and kurtosis values for each period
    #if you want to parallelize, start Julia with "julia -n <number of cores>"
    output = periodogram(data, periods, parallel=true) 
    
    #do some postprocessing 
    output[:delt_chi2] = flatten(output[:chi2], periods) #calculate \Delta \chi^2 (see equation N)

    null_output = scrambled_periodogram(data, periods)
    null_output[:delt_chi2] = flatten(null_output[:chi2], periods)
    sigma = movingstd((null_output[:kurtosis] .- 3) .* (null_output[:delt_chi2]))

    output[:zeta] = (output[:kurtosis] .- 3) .* output[:delt_chi2] ./ sigma

    #plot the periodogram
    using PyPlot
    plot(periods, output[:zeta])
