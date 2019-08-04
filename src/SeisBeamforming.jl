module SeisBeamforming

using SeisIO, SeisNoise, JLD2, FileIO, PlotlyJS, Geodesy, Statistics, LinearAlgebra, Dates

export seisnoisebeamform

#include("./src/Types/BeamData.jl")
#include("./src/Types/show.jl")

include("./Types/BeamData.jl")
include("./Types/show.jl")

function seisnoisebeamform(InputDict::Dict)

    finame 				= InputDict["finame"]
	fopath        		= InputDict["fopath"]
	IsPlotFigure = InputDict["IsPlotFigure"]

    stationlist 		= load(finame, "info/stationlist")
	DLtimestamplist 	= load(finame, "info/DLtimestamplist")
	DL_time_unit		= load(finame, "info/DL_time_unit")
	starttime 			= load(finame, "info/starttime")
	endtime 			= load(finame, "info/endtime")

	channel			= InputDict["channel"]
	bp_len 			= InputDict["bp_len"] #hourly beamforming
	bp_step 		= InputDict["bp_step"] # 30 minutes overlap
	freq_min 		= InputDict["freq_min"]
	freq_max 		= InputDict["freq_max"]
	sampling_freq 	= InputDict["sampling_freq"]
	freq_range 		= InputDict["freq_range"]
	fexamin_num 	= InputDict["fexamin_num"] # num of examined freqency
	vel_range 		= InputDict["vel_range"]
	Δvel   			= InputDict["Δvel"]
	Δθ    			= InputDict["Δθ"]

    # store station location
    lat = Float64[]; lon = Float64[]; x = Float64[]; y = Float64[];

    for i = 1:length(stationlist)
        Stemp = load(finame, joinpath(DLtimestamplist[1], stationlist[i]))
        if any(occursin.(channel, stationlist[i][end-2:end]))
            push!(lat, Stemp.loc.lat)
            push!(lon, Stemp.loc.lon)
            #- ignore elevation -#
        end
    end

	if isempty(lat); error("no available station to compute beampower in this array. abort."); end

	# project lat and lon to x-y coordinate
    origin_lla = LLA(mean(lat), mean(lon))
    trans = ENUfromLLA(origin_lla, wgs84)

    for i = 1:length(stationlist)
        point_lla = LLA(lat[i], lon[i])
        point_enu = trans(point_lla)
        push!(x, point_enu.e)
        push!(y, point_enu.n)
    end

    # check with plotting
	if IsPlotFigure
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
	                 #xaxis_range=[-15, 15],
	                 #yaxis_range=[-15, 15],
	                 xaxis=attr(title="x [km]"),
					 yaxis=attr(title="y [km]"))
	    p = PlotlyJS.plot([trace1, trace2], layout)
		display(p)
	end


	# save info into jld2
	jldopen(fopath, "w") do file
		file["info/DLtimestamplist"]  = DLtimestamplist;
		file["info/stationlist"]  = stationlist;
		file["info/starttime"]   = starttime;
		file["info/endtime"]     = endtime;
	end

	# Loop with DL time sample; this becomes unit duration of bpstack
	for i = 1:length(DLtimestamplist)
		println(DLtimestamplist[i])

		# make beampower timestamp: bptimestamp
		t1 = parse.(Int64, split(DLtimestamplist[i], ".")[1:2])
		hm = split(DLtimestamplist[1], ".")[3]
		m1, d1 = SeisIO.j2md(t1[1], t1[2])
		st = DateTime(join([t1[1], "-", m1, "-",  d1, hm]))

		bptimestamp = []
		ttemp = deepcopy(st)

		while ttemp <= st + Second(DL_time_unit) - Second(bp_len)
			push!(bptimestamp, ttemp)
			ttemp += Second(bp_step)
		end

		# compute fft for all stations
		FFT1 = Array{FFTData, 1}(undef, length(stationlist))

		for j = 1:length(stationlist)
			Stemp = load(finame, joinpath(DLtimestamplist[i], stationlist[j]))
			Stemp.t[1,2]  = round(Stemp.t[1,2], sigdigits=13) #to avoid t_len skip
			try
				FFT1[j] = compute_fft(SeisData(Stemp), freq_min, freq_max,
			 	sampling_freq, Int(bp_step), Int(bp_len))
			catch # unknow error with FFT1

				FFT1[j] = FFTData(t = [0.0])
				println("skip station $(stationlist[j])")
			end
		end

		# assume all station have same sampling frequency
		L = sampling_freq * bp_len
		fvec0 = sampling_freq * collect(0:L/2)/L

		# compute examined frequencies
		fexam = exp10.(range(log10(freq_range[1]),stop=log10(freq_range[2]),length=fexamin_num))
		exfreq = []
		#@doc "exfreq contains (id, frequency) in the FFT.fft array to compute Y array." exfreq
		for j = 1:length(fexam)
			id1 = findfirst(x -> x>=fexam[j], fvec0)
			push!(exfreq, (Int(id1), fvec0[id1]))
		end

		# compute examined phase velocity and θ
		exvel = collect(vel_range[1]:Δvel:vel_range[2])
		exθ = deg2rad.(collect(0:Δθ:360))

		# assign size of BeamData.bp
		size_v = length(exvel)
		size_θ = length(exθ)
		size_t = length(bptimestamp)

		Beam = Array{BeamData, 1}(undef, length(exfreq))

		for l = 1:length(exfreq) #loop every examined frequency

			Btemp = BeamData()
			Btemp.name = join(channel, ",")
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

				# compute size of e and Y at this loop: i.e. gather available stations at this timewindow
				Y_station_and_timestampid = []
				#@doc "Y_station_and_timestampid contains the station id and timestamp id corresponding to bptimestamp[k]
				# to compute Y array." Y_station_and_timestampid

				for stid = 1:length(stationlist)
					for timestampid = 1:length(FFT1[stid].t)
						if u2d(FFT1[stid].t[timestampid]) == bptimestamp[k]
							#this station is used as one station in seismic array
							push!(Y_station_and_timestampid, (stid, timestampid))
						end
					end
				end

				if length(Y_station_and_timestampid) >= 2 # if stations are more than 2 in array

					# loop below is to compute one polar coord figure
					stack_count += 1

					for eθid in 1:length(exθ)
						eθ = exθ[eθid]

						for evid = 1:length(exvel) #loop examined phase velocity
							ev = exvel[evid]

							# examined frequency
							ω_id = exfreq[l][1]
							ω	 = 2 * pi * exfreq[l][2] # ω = angular frequency
							kv   = (ω / ev) .* [cos(eθ), sin(eθ)]

							#compute the array plane wave response [exp(-ikv.xN)]^*T
							xv = Float32[]
							ee = ComplexF32[]

							for stid = 1:length(Y_station_and_timestampid)
								stationid = Y_station_and_timestampid[stid][1]
								xv= [x[stationid]; y[stationid]]
								push!(ee, exp(-im*dot(kv, xv)))
							end

							# compute spatial correlation matrix YY
							Y = ComplexF32[]

							for stid = 1:length(Y_station_and_timestampid)
								stationid = Y_station_and_timestampid[stid][1]
								timestampid = Y_station_and_timestampid[stid][2]
								f1 = FFT1[stationid].fft[ω_id, timestampid]
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

		# save info into jld2
		jldopen(fopath, "r+") do file
			varname = joinpath("BeamData", DLtimestamplist[i])
			file[varname]  = Beam;
		end
	end

	print("Beam forming is successfully done at: ")
	println(now())

end


end # module
