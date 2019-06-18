import Base:in, +, -, *, ==, convert, isempty, isequal, length, push!, sizeof, append!

export BeamData

const beamfields = [:name, :bp_len,:bp_step, :f, :c, :θ]

"""
   BeamData
A structure for Beam power of ambient noise data.
## Fields: BeamData
| **Field** | **Description** |
|:--------------------|:--------------------|
| :name       | name of input array |
| :misc         | misc dictionary. |
| :notes        | notes |
| :bp_len         | time window length to compute beam power |
| :bp_step       | time window step to compute beam power  |
| :f    | examined frequency [1/s] |
| :c    | examined phase velocity [m/s] |
| :θ    | examined aximuth of slowness vector [rad] |
| :t    | time stamp of each time window |
| :bp   | beam power b(f, t, θ, c)|
| :tstack  | time stamp of stacked data. |
| :tstack_len       | time window length of stacked beam power. |
| :bpstack      | stacked beam power. |
"""
mutable struct BeamData

  name     ::String                          # name of input file
  misc     ::Dict{String,Any}                # misc
  notes    ::Array{String,1}                 # notes
  bp_len   ::Int                             # beam power time window length [s]
  bp_step  ::Int                             # beam power time window step   [s]
  f::Array{Float64,1}                         # examined frequency [1/s]
  c::Array{Float64,1}                        # examined wave velocity
  θ::Array{Float64,1}                        # examined azimuth
  t::Array{Float64,1}                         # time
  bp::Array{<:Union{Complex{Float32},Complex{Float64}},3}  # beam power data
  tstack::Float64                             #stacked time
  tstack_len::Int                             #stacked time length
  bpstack::Array{Complex{Float32},2}          #stacked beam power

  function BeamData(
      name     ::String,
      misc     ::Dict{String,Any},
      notes    ::Array{String,1},
      bp_len   ::Int,
      bp_step  ::Int,
      f        ::Array{Float64,1},
      c        ::Array{Float64,1},
      θ        ::Array{Float64,1},
      t        ::Array{Float64,1},
      bp       ::Array{<:Union{Complex{Float32},Complex{Float64}},3},
      tstack   ::Float64,
      tstack_len::Int,
      bpstack  ::Array{Complex{Float32},2}
      )

      return new(name, misc, notes, bp_len, bp_step, f, c, θ, t, bp, tstack, tstack_len, bpstack)
    end
end

BeamData(;
          name     ::String                    = "",
          misc     ::Dict{String,Any}          = Dict{String,Any}(),
          notes    ::Array{String,1}           = Array{String,1}(undef, 0),
          bp_len   ::Int                       = zero(Int),
          bp_step  ::Int                       = zero(Int),
          f        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          c        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          θ        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          t        ::Array{Float64,1}          = Array{Float64,1}(undef, 0),
          bp       ::Array{<:Union{Complex{Float32},Complex{Float64}},3} =
                     Array{Complex{Float32},3}(undef, 0, 0, 0),
          tstack    ::Float64                   = zero(Float64),
          tstack_len::Int                       = zero(Int),
          bpstack  ::Array{Complex{Float32},2}  = Array{ComplexF32,2}(undef, 0, 2)

          ) = BeamData(name, misc, notes, bp_len, bp_step, f, c, θ, t, bp, tstack, tstack_len, bpstack)


in(s::String, BP::BeamData) = BP.name==s

isempty(B::BeamData) = min(B.name=="",
                           B.bp_len == zero(typeof(B.bp_len)),
                           B.bp_step == zero(typeof(B.bp_step)),
                           B.tstack == zero(typeof(B.tstack)),
                           B.tstack_len == zero(typeof(B.tstack_len)),
                           isempty(B.f),isempty(B.c),
                           isempty(B.θ),isempty(B.t),
                           isempty(B.bp),isempty(B.bpstack),
                           isempty(B.misc),isempty(B.notes))

isequal(B::BeamData, U::BeamData) = minimum([hash(getfield(F,i))==hash(getfield(U,i)) for i in beamfields]::Array{Bool,1})
==(B::BeamData, U::BeamData) = isequal(B,U)::Bool

append!(B::BeamData, U::BeamData)  = (if B == U;
  setfield!(B, :t, vcat(getfield(B,:t), getfield(U,:t)));
  setfield!(B, :bp, hcat(getfield(B,:bp), getfield(U,:bp)))
  elseif isempty(B);
  [setfield!(B, i, getfield(U,i)) for i in beamfields];
  [setfield!(B, i, getfield(U,i)) for i in [:t,:bp]];
  end;
  return B)

+(B::BeamData, U::BeamData) = (T = deepcopy(B); return append!(T, U))

function sizeof(B::BeamData)
s = sum([sizeof(getfield(B,f)) for f in beamfields])
if !isempty(B.notes)
  s += sum([sizeof(i) for i in B.notes])
end
if !isempty(B.misc)
  s += sum([sizeof(i) for i in values(B.misc)])
end
return s
end
