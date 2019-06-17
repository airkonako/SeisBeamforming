module SeisBeamforming

using SeisIO, Noise, JLD2

export seisnoisebeamform

include("Types/BeamData.jl")
include("Types/show.jl")

function seisnoisebeamform()

    BP = BeamData()


end

end # module
