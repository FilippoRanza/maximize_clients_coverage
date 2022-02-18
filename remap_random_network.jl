using HDF5
using JSON


json_data = JSON.parsefile("example-network.json")
points = json_data["points"]

mat_points = zeros(length(points), 2)
for (i, pt) in enumerate(points)
    mat_points[i, 1] = pt[1]
    mat_points[i, 2] = pt[2]
end


mat_points .+= 2500
mat_points .*= (150 / 5000)

h5open("stations.hdf5", "w") do file
    file["network/1"] = mat_points
end





