using SeisIO, PlotlyJS, JLD2, Dates, Noise, Printf, ORCA

fi = "./dataset/BPnetwork.jld2"

t = jldopen(fi)

DLtimestamplist = t["info/DLtimestamplist"]
stationlist = t["info/stationlist"]

tid  = 1
for tid = 1:length(DLtimestamplist)
    plotspan=1000
    # make beampower timestamp: bptimestamp
    t1 = parse.(Int64, split(DLtimestamplist[tid], ".")[1:2])
    hm = split(DLtimestamplist[1], ".")[3]
    m1, d1 = SeisIO.j2md(t1[1], t1[2])
    st = DateTime(join([t1[1], "-", m1, "-",  d1, hm]))
    tstamp = Noise.timestamp(st)


    layout = Layout(width=800, height=600,
                    xaxis=attr(title="Time [hours]", tickvals =[0, 6, 12, 18, 24]),
                    yaxis=attr(tickvals = collect(1:length(stationlist)),
                    ticktext=stationlist, automargin=true),
                    #font =attr(size=12),
                    showlegend=false,
                    xaxis_range=[0, 24],
                    #yaxis_range=[],
                    title = @sprintf("Time:%s", tstamp)
                    )


    p = PlotlyJS.plot([NaN], layout)

    amp = 1e4
    #plot trace
    for i = 1:length(stationlist)
        station = stationlist[i]
        S = t[joinpath(DLtimestamplist[tid], station)]
        tvec = collect(0:length(S.x)-1)/S.fs / 60 / 60
        trace1 = scatter(;x=tvec[1:plotspan:end], y=amp*S.x[1:plotspan:end].+float(i), mode="lines")
        PlotlyJS.addtraces!(p, trace1)
    end

    PlotlyJS.deletetraces(p,1)

    #display(p)

    mkpath("./fig_dawdata")
    figname = @sprintf("./fig_dawdata/%s.png",tstamp)
    savefig(p, figname)
end
