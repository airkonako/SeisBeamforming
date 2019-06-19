module SeisBeamforming

using SeisIO, Noise, JLD2, FileIO, PlotlyJS, Geodesy, Statistics, LinearAlgebra, Dates

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
	DL_time_unit= load(finame, "info/DL_time_unit")
	starttime = load(finame, "info/starttime")
	endtime = load(finame, "info/endtime")

	bp_len = 60*60.0 # hourly beamforming
	bp_step = 60*30 # 30 minutes overlap
	freq_min = 0.1
	freq_max = 5.0
	fexamin_num = 30 # num of examined freqency


    # store station location
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


	# save info into jld2
	jldopen(fopath, "w") do file
		file["info/DLtimestamplist"]  = DLtimestamplist;
		file["info/stationlist"]  = stationlist;
		file["info/starttime"]   = starttime;
		file["info/endtime"]     = endtime;
	end

	for i = 1:length(DLtimestamplist) # this becomes unit duration of bpstack

		# make beampower timestamp

		i = 2
		t1 = parse.(Int64, split(DLtimestamplist[i], ".")[1:2])
		hm = split(DLtimestamplist[1], ".")[3]
		m1, d1=SeisIO.j2md(t1[1], t1[2])
		st = DateTime(join([t1[1], "-", m1, "-",  d1, hm]))

		bptimestamp = []
		ttemp = deepcopy(st)
		while ttemp <= st + Second(DL_time_unit) - Second(bp_len)
			push!(bptimestamp, ttemp)
			global ttemp += Second(bp_step)
		end

		# compute fft for all stations
		FFT1 = Array{FFTData, 1}(undef, length(stationlist))

		for j = 1:length(stationlist)
			Stemp = load(finame, joinpath(DLtimestamplist[i], stationlist[j]))
			Stemp.t[1,2]  = round(Stemp.t[1,2], sigdigits=13) #to avoid t_len skip
			FFT1[j] = compute_fft(SeisData(Stemp), freq_min, freq_max,
			 Stemp.fs, Int(bp_step), Int(bp_len))
		end

		# compute examined f, c, θ
		freq_range = [0.2, 3.0]
		vel_range = [1000, 4000]
		Δvel   = 500.0
		Δθ    = 10

		# assume all station have same sampling frequency
		L = FFT1[1].fs * bp_len
		fvec0 = FFT1[1].fs*collect(0:L/2)/L
		fexam = exp10.(range(log10(freq_range[1]),stop=log10(freq_range[2]),length=fexamin_num))

		exfreq = []
		@doc "exfreq contains (id, frequency) in the FFT.fft array to compute Y array." exfreq

		for j = 1:length(fexam)

			id1 = findfirst(x -> x>=fexam[j], fvec0)

			push!(exfreq, (Int(id1), fvec0[id1]))
		end

		exvel = collect(vel_range[1]:Δvel:vel_range[2])
		exθ = deg2rad.(collect(0:Δθ:360-Δθ))

		# assign size of BeamData.bp
		size_v = length(exvel)
		size_θ = length(exθ)
		size_t = length(bptimestamp)

		Beam = Array{BeamData, 1}(undef, length(exfreq))

		for l = 1:length(exfreq) #loop every examined frequency
			l = 1

			Btemp = BeamData()
			Btemp.name = channel
			Btemp.bp_len = bp_len
			Btemp.bp_step = bp_step
			Btemp.f = [exfreq[l][2]]
			Btemp.c = exvel
			Btemp.θ = exθ
			Btemp.t = d2u.(bptimestamp)
			Btemp.bp = Array{Complex{Float32},3}(undef, size_θ, size_v, size_t)
			Btemp.tstack = d2u(bptimestamp[round(Int64,length(bptimestamp)/2)+1]) #mean time of stack duration
			Btemp.tstack_len = parse(Int64, DL_time_unit)
			Btemp.bpstack = zeros(ComplexF32, size_θ, size_v)

			stack_count = 0
			for k = 1:length(bptimestamp) #loop every bp_len
				k = 2
				# compute size of e and Y at this loop: i.e. gather available stations at this timewindow
				Y_station_and_timestampid = []
				@doc "Y_station_and_timestampid contains the station id and timestamp id corresponding to bptimestamp[k]
				 to compute Y array." Y_station_and_timestampid

				for stid = 1:length(stationlist)
					for timestampid = 1:length(FFT1[stid].t)
						if u2d(FFT1[stid].t[timestampid]) == bptimestamp[k]
							#this station is used as one station in seismic array
							push!(Y_station_and_timestampid, (stid, timestampid))
						end
					end
				end

				if Y_station_and_timestampid >= 2

					# loop below is to compute one polar coord figure
					global stack_count += 1

					for eθid in 1:length(exθ)
						eθ = exθ[eθid]

						for evid = 1:length(exvel) #loop examined phase velocity
							ev = exvel[evid]

							omega_id = exfreq[l][1]
							kv = (exfreq[l][2] / ev) .* [cos(eθ), sin(eθ)]

							size_Y = length(Y_station_and_timestampid)

							#compute the array plane wave response [exp(-ikv.xN)]^*T
							xv = Float64[]
							ee = ComplexF64[]

							for stid = 1:length(Y_station_and_timestampid)
								stationid = Y_station_and_timestampid[stid][1]
								xv= [x[stationid]; y[stationid]]
								push!(ee, exp(-im*dot(kv, xv)))
							end

							# compute spatial correlation matrix YY
							Y = ComplexF64[]

							for stid = 1:length(Y_station_and_timestampid)
								stationid = Y_station_and_timestampid[stid][1]
								timestampid = Y_station_and_timestampid[stid][2]
								f1 = FFT1[stationid].fft[exfreq[l][1], timestampid]
								push!(Y, f1)
							end

							# compute e'YY'e
							YY = Y*transpose(conj.(Y))

							Btemp.bp[eθid, evid, k] = transpose(conj.(ee)) * YY * ee

							#add into stack
							Btemp.bpstack[eθid, evid] += transpose(conj.(ee)) * YY * ee
						end
					end


				else
					println("$(bptimestamp[k]): num of station in array is less than two; thus skip this beam forming.")
					Btemp.bp[:, :, k] = zeros(ComplexF32, size_θ, size_v)
				end
			end

			# normalize stacked beam power by number of stacking
			Btemp.bpstack = Btemp.bpstack ./ stack_count

			# store into beam data struct
			Beam[l] = Btemp
		end

		# save into jld2 file

		varname = DLtimestamplist[i]
		save(fopath, Dict(varname => Beam))

		# save info into jld2
		jldopen(fopath, "r+") do file
			varname = joinpath("BeamData", DLtimestamplist[i])
			file[varname]  = Beam;
		end

end


end # module
