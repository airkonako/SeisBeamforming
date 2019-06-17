import Base:show, size, summary

const show_os = 12
const unindexed_fields = (:v, :n)

summary(B::BeamData) = string(typeof(B), " with ", length(B.t), " beampower",
                             size(B.bp,2) == 1 ? "" : "s")

function show(io::IO, B::BeamData)
  W = max(80, displaysize(io)[2]) - show_os
  w = min(W, 35)
  nc = size(B.bp,1) == 0 ? 0 : 1
  N = min(nc, div(W-1, w))
  M = min(N+1, nc)
  bplen,numbp = size(B.bp)
  println(io, "BeamData with ", length(B.t), " beampowers")
  Field = fieldnames(BeamData)
  for f in Field
    if (f in unindexed_fields) == false
      targ = getfield(B, f)
      t = typeof(targ)
      fstr = uppercase(String(f))
      print(io, lpad(fstr, show_os-2), ": ")
      if f == :notes || f == :misc
        SeisIO.show_str(io, String[string(length(targ), " entries")], w, N)
      elseif (t <: AbstractFloat || t <: InstrumentPosition || t<: InstrumentResponse)
        println(io, repr("text/plain", targ, context=:compact=>true))
      elseif f == :t
        if length(B.t) > 0
          SeisIO.show_str(io, String[timestamp(B.t[1]), " (", string(length(B.t)), " BeamPowers)"], w, N)
        else
          SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
        end
      elseif f == :bp
        SeisIO.show_str(io, String[summary(targ)], w, N)
      else
        SeisIO.show_str(io, String[repr("text/plain", targ, context=:compact=>true)], w, N)
      end
    end
  end
  return nothing
end
show(B::BeamData) = show(stdout, B)
