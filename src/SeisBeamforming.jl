module SeisBeamforming

using SeisIO, Noise, JLD2, FileIO, PlotlyJS, Geodesy, Statistics, LinearAlgebra

export seisnoisebeamform

include("./src/Types/BeamData.jl")
include("./src/Types/show.jl")

function seisnoisebeamform(InputDict::Dict)

    InputDictionary = Dict([
            "freq_range"=> freq_range,
            "vel_range" => vel_range,
            "Δvel"   => Δvel,
            "Δθ"     => Δθ,
            "finame"=> DL_time_unit,
            "fopath"          => fopath
        ])

    #finame = InputDict["finame"]
    finame = "./EXAMPLE/Noise_BP/dataset/BPnetwork.jld2"
    t = load(finame, "info")
    stationlist = load(finame, "info/stationlist")
    DLtimestamplist = load(finame, "info/DLtimestamplist")
    starttime = load(finame, "info/starttime")
    #store station location
    lat = Float64[]
    lon = Float64[]
    x   = Float64[]
    y   = Float64[]

    channel = "BP1"

    for i = 1:length(stationlist)
        Stemp = load(finame, joinpath(DLtimestamplist[1], stationlist[i]))
        if occursin(channel, stationlist[i])
            push!(lat, Stemp.loc.lat)
            push!(lon, Stemp.loc.lon)
            # ignore elevation
        end
    end

	if isempty(lat); error("no available station to compute beampower"); end

    origin_lla = LLA(mean(lat), mean(lon))
    trans = ENUfromLLA(origin_lla, wgs84)

    for i = 1:length(stationlist)
        point_lla = LLA(lat[i], lon[i])
        point_enu = trans(point_lla)
        push!(x, point_enu.e)
        push!(y, point_enu.n)
    end
    #check with plotting
    trace1 = PlotlyJS.scatter(;x=x./1e3, y=y./1e3,
                          mode="markers", name="Stations",
						  textposition="top center",
                          text=stationlist,
                          marker_size=12)
    trace2 = PlotlyJS.scatter(;x=[0.0], y=[0.0],
						mode="markers",
 						marker=attr(symbol=:star),
						name="Origin",
						marker_size=12)
    layout = Layout(;width=500, height=500,
				 showlegend=true,
                 xaxis_range=[-15, 15],
                 yaxis_range=[-15, 15],
                 xaxis=attr(title="x [km]"),
				 yaxis=attr(title="y [km]"))
    PlotlyJS.plot([trace1, trace2], layout)


	for i = 1:length(DLtimestamplist) # this becomes unit duration of bpstack
		bp_len = 60*60.0 # hourly beamforming
		bp_step = 60*30 # 30 minutes overlap
		freq_min = 0.1
		freq_max = 5.0
		fexamin_num = 30 # num of examined freqency

		# compute fft for all stations
		FFT1 = Array{FFTData, 1}(undef, length(stationlist))

		# compute examined f, c, θ
		freq_range = [0.2, 3.0]
		vel_range = [1000, 4000]
		Δvel   = 500.0
		Δθ    = 10

		#assume all station have same sampling frequency
		L = FFT1[1].fs * bp_len
		fvec0 = FFT1[1].fs*collect(0:L/2)/L
		fexam = exp10.(range(log10(freq_range[1]),stop=log10(freq_range[2]),length=fexamin_num))

		exfreq = []

		for j = 1:length(fexam)

			id1 = findfirst(x -> x>=fexam[j], fvec0)

			push!(exfreq, (Int(id1), fvec0[id1]))
		end

		exvel = collect(vel_range[1]:Δvel:vel_range[2])
		exθ = deg2rad.(collect(0:Δθ:360-Δθ))

		for j = 1:length(stationlist)
			Stemp = load(finame, joinpath(DLtimestamplist[i], stationlist[j]))
			Stemp.t[1,2]  = round(Stemp.t[1,2], sigdigits=13) #to avoid t_len skip
			FFT1[j] = compute_fft(SeisData(Stemp), freq_min, freq_max,
			 Stemp.fs, Int(bp_step), Int(bp_len))
		end

		kcount = 1
		for k = 1:length(FFT1[1].t) #loop every bp_len

			for l = 1:length(exfreq) #loop every examined frequency

				# loop below compute one polar coord figure
				for ev = exvel #loop examined phase velocity

					for eθ in exθ

						omega_id = exfreq[m][1]
						kv = (exfreq[m][2] / ev) .* [cos(eθ), sin(eθ)]

						# compute size of e and Y at this loop
						Y_stationid = Int[]
						for stid = 1:length(stationlist)
							try
								FFT1[tid].t[k]
								push!(Y_stationid, tid)
							catch
								continue;
							end
						end

						size_Y = length(Y_stationid)

						#compute the array plane wave response [exp(-ikv.xN)]^*T
						xv = Float64[]
						ee = ComplexF64[]
						for stid = Y_stationid
							xv= [x[stid]; y[stid]]
							push!(ee, exp(-im*dot(kv, xv)))
						end

						# compute spatial correlation matrix YY
						Y = ComplexF64[]
						for stid = Y_stationid
							f1 = FFT1[stid].fft[exfreq[l][1]]
							push!(Y, f1)
						end

						#reshape e and Y following Brooks et al. (2009)

						ee = transpose(conj.(ee))
						Y = transpose(conj.(Y))

						# compute YY'
						YY = Y*transpose(conj.(Y))

						#add to BeamData





    BP = BeamData()




end

end # module
