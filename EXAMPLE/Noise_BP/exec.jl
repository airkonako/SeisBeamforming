#==================================================#
# Input Parameters


freq_range = [0.01, 1.0] # frequency range for evaluating beam power: used to choose omega from FFTfreq [1/s]
vel_range  = [1000, 5000]      # evaluating velocity range [cmin cmax] [m/s]

Δvel = 500 # grid space of velocity [m/s]
Δθ   = 10 # grid space of azimuth of trial slowness vector [deg]



finame      = "./dataset/BPnetwork.jld2"
fodir   = "./dataset"
foname  = "BP_beamform" # data is saved at ./dataset/$foname.jld2
#==================================================#

mkpath(fodir)
fopath=joinpath(fodir, foname*".jld2")

InputDictionary = Dict([
        "freq_range"=> freq_range,
        "vel_range" => vel_range,
        "Δvel"   => Δvel,
        "Δθ"     => Δθ,
        "finame"=> DL_time_unit,
        "fopath"          => fopath
    ])

# mass request with input Dictionary
seisdownload(InputDictionary, MAX_MEM_PER_CPU=float(MAX_MEM_PER_CPU))
