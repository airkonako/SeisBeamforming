
#include("./src/SeisBeamforming.jl")
include("/Users/kurama/Documents/kurama/research/SCECInternship2019/Beamforming/SeisBeamforming/src/SeisBeamforming.jl")
using .SeisBeamforming
using SeisIO, Noise, JLD2, FileIO, Dates, Printf, HDF5
#==================================================#
# Input Parameters

finame      = "./dataset/BP_beamform.jld2"
fodir   = "./dataset"
foname  = "BP_polarplot" # data is saved at ./dataset/$foname.jld2

#==================================================#

mkpath(fodir)
fopath=joinpath(fodir, foname*".h5")

t = jldopen(finame)

DLtimestamplist = t["info/DLtimestamplist"]



h5open(fopath, "w") do file
    write(file, "info/DLtimestamplist", string.(DLtimestamplist))  # alternatively, say "@write file A"
end

for tstamp = DLtimestamplist

    Beam = t[joinpath("BeamData", tstamp)]

    # obtain examined frequency vector
    fvec = Float64[]
    fvecstamplist = String[]
    for freqid = 1:length(Beam)
        push!(fvec, Beam[freqid].f[1])
        push!(fvecstamplist, @sprintf("f:%04.2fHz", Beam[freqid].f[1]))
    end


    group1 = joinpath("BeamPower", tstamp)
    h5open(fopath, "r+") do file
        write(file, group1*"/fvec",  fvec)
        write(file, group1*"/fvecstamplist",  fvecstamplist)
    end


    for freqid = 1:length(Beam)
        B = Beam[freqid]
        cvec = B.c
        θvec = round.(rad2deg.(B.θ), digits=3)
        freq = B.f
        time = string.(u2d.(B.t))

        group1 = joinpath("BeamPower", tstamp, fvecstamplist[freqid])
        h5open(fopath, "r+") do file
            write(file, group1*"/cvec",  cvec)
            write(file, group1*"/thetavec",  θvec)
            write(file, group1*"/f",  fvec[freqid])
            write(file, group1*"/time",  time)
        end


        for tid = 1:length(time)
            bp = B.bp[:, : , tid]

            group1 = joinpath("BeamPower", tstamp, fvecstamplist[freqid], time[tid])
            h5open(fopath, "r+") do file
                write(file, group1*"/bp",  abs.(bp))
            end
        end
    end

end
