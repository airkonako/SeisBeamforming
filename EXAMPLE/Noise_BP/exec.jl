
#include("./src/SeisBeamforming.jl")
include("/Users/kurama/Documents/kurama/research/SCECInternship2019/Beamforming/SeisBeamforming/src/SeisBeamforming.jl")
using .SeisBeamforming
using SeisIO, Noise, JLD2, FileIO, PlotlyJS, Geodesy, Statistics, LinearAlgebra, Dates

#==================================================#
# Input Parameters

channel = ["BP1"] # channel name to do beam forming; ["BP1", "HHZ" is available]

bp_len = 60 * 60 # time window to conduct beam forming
bp_step = 60 * 30 # overlap of time window
freq_min = 0.1 # low freq at band pass filter before FFT
freq_max = 5.0 # high freq at band pass filter before FFT
sampling_freq = 20.0 # sampling frequency used to downsampling at FFT
# examined parameters
freq_range = [0.2, 3.0] # frequency range for evaluating beam power: used to choose omega from FFTfreq [1/s]
fexamin_num = 10# num of examined frequency

vel_range  = [1000, 5000]      # evaluating velocity range [cmin cmax] [m/s]
Δvel = 200 # grid space of velocity [m/s]

Δθ   = 10 # grid space of azimuth of trial slowness vector [deg]



finame      = "./dataset/BPnetwork.jld2"
fodir   = "./dataset"
foname  = "BP_beamform" # data is saved at ./dataset/$foname.jld2

IsPlotFigure = false # plot figures to debug
#==================================================#

mkpath(fodir)
fopath=joinpath(fodir, foname*".jld2")

InputDictionary = Dict([
        "channel"  	=> channel,
	"bp_len"	=> bp_len,
	"bp_step"	=> bp_step,
	"freq_min"	=> freq_min,
	"freq_max"	=> freq_max,
	"sampling_freq"	=> sampling_freq,
        "freq_range"	=> freq_range,
	"fexamin_num"	=> fexamin_num,
        "vel_range" 	=> vel_range,
        "Δvel"   	=> Δvel,
        "Δθ"     	=> Δθ,
        "finame"	=> finame,
        "fopath"        => fopath,
	"IsPlotFigure"  => IsPlotFigure
    ])

# mass request with input Dictionary
seisnoisebeamform(InputDictionary)
