This code implements the algorithm in [Wheeler & Kipping 201?]().  Find posterior samples from the paper at [github link tk]().

This code is tested with Julia 1.0, but it should be compatible with 0.6 as well.  It requires the `Interpolations`, `DataFrames`, and `FITSIO` packages.  To install them run these commands on the Julia REPL

    Pkg.add("Interpolations")
    Pkg.add("DataFrames")
    Pkg.add("FITSIO")

Here is a brief usage example

    include("WeirdDetector")
    using WeirdDetector #make key functions accessible without being qualified
    
    #here's how you can get a Kepler light curve detrended like we did in the paper
    getFITS(8462852, fitsdir="./") #download FITS files for Boyajian's star
    df = loadFITS(8462852, fitsdir="./") #load all quarters into single data frame, perform ourlier rejection and detrending
    data = pointsify(df) #convert to data type taken by periodogram()

    #Here's how to prep arbitary data
    t = ...      # array of epochs
    F = ...      # array of flux values
    sigmaF = ... # array of 1-sigma uncertainties
    data = Point.(t, F, sigmaF) #this constructs an array of Points with Julia's broadcast syntax

    periods = optimal_periods(pmin=0.25, pmax=50) #construct array of periods used in paper.  Any array of Float32's will do.
    output = periodogram(data, periods, parallel=true) #construct a DataFrame containing the chi-squared and kurtosis values for each period
                                                       #if you want to parallelize, start Julia with "julia -n <number of cores>"
    
    #do some postprocessing 
    
