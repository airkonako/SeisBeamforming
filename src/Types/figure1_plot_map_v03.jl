using Printf, GMT, JLD2

west = -122
east = -119.5
south = 35
north = 37

basemapregion ="$west/$south/$east/$(north)r"
frame=(axes=:NSEW, annot=1, ticks="30m",grid="15m")
center = [-120.0 36.0]
#proj = (name=:Cassini, center=center)
proj = (name=:Mercator)
Ispace = "0.1m"
color="./GMT_relief.cpt"

figsize = 14

#make large xyz file
#grdlist = ["grdn36w120", "grdn36w121"]
grdlist = ["grdn36w120", "grdn36w121", "grdn36w122", "grdn37w120", "grdn37w121", "grdn37w122"]
#
# for gg = grdlist
#    println(gg)
#    #convert it to grd file
#    fi = "./Arcgrid/xyzfiles/$gg"
#    G1 = xyz2grd(fi, region=basemapregion, proj=proj, I=Ispace)
#
#    #---process grd data---#
#    #remove inf in z
#    replace!(G1.z, -Inf=>NaN)
#    #replace!(G1.z, 0=>-2000) this cause bug
#
#    #======================#
#
#    #save into grd file
#    gmtwrite("./$gg.grd", G1)
# end

#
# for gg = grdlist
#
#    #convert it to grd file
#    fi = "./$gg.grd"
#    G1 = gmtread(fi, grd=true)
#
#    #save into grd file
#    gmtwrite("./$(gg)_lim.grd", G1)
# end

# plot

icount = 0
for gg = grdlist
   global icount += 1
   fi = "./$(gg).grd"
   println(fi)
   if icount == 1
      grdimage(fi, region=basemapregion, proj=proj, frame=frame,figsize=figsize,  color=color, coast=true, shade=(auto=true,), nan_alpha=true)
   elseif icount == length(grdlist)
      grdimage!(fi, region=basemapregion, proj=proj, frame=frame,figsize=figsize,  color=color, coast=true, shade=(auto=true,), nan_alpha=true)
   else
      grdimage!(fi, region=basemapregion, proj=proj, frame=frame,figsize=figsize,  color=color, coast=true, shade=(auto=true,), nan_alpha=true)
   end
end


# plot fault
t = jldopen("sanandreasfault.jld2")
fault = t["fault"]


for i = 1:length(fault)
#for i = 1:10

   f1 = fault[i]

   x = Float64[]
   y = Float64[]
   for j = 1:length(f1)
      push!(x, f1[j][1])
      push!(y, f1[j][2])
   end

   # select faults
   if  all(x -> x >= west, x) && all(x -> x <= east, x) && all(x -> x >= south, y) && all(x -> x <= north, y)
      println(x, y)
      xy = mat2ds(y, x=x)
      lines!(xy, pen=(1.5,:red)) # you can ignore error with this line
   end
end

#colorbar!(pos=(width=10, justify=:TC, outside=true, offset=(-6, 11), horizontal=true),
#   frame=(annot=1000, ticks=500),  truncate=(0, 3000), cmap=color, show=1)

# plot stations
station = jldopen("stationcoords.jld2")

# NC
NC = station["NC"]

for i = 1:length(NC)
   GMT.plot!([NC[i][1] NC[i][2]], marker="d", markersize=0.35, pen=(0.6),
   markeredgecolor=:black, markerfacecolor="135/206/250")
end

# BP
BP = station["BP"]

for i = 1:length(BP)
   GMT.plot!([BP[i][1] BP[i][2]], marker="v", markersize=0.35, pen=(0.8),
   markeredgecolor=:black, markerfacecolor="251/224/0")
end


#OFFICIAL WEB SITE SAYS INCOLECT VARIABLE NOM: scale_lat -> scale_at_lat
mapscale=(map=true, anchor=[-121.5, 35.2], scale_at_lat=36,  length=50, fancy=true, label="[km]")
p = coast!(region=basemapregion, shore=(pen=(0.8,:black)), L=mapscale,
 frame=frame, S="0/141/203", stamp="", show=true)
