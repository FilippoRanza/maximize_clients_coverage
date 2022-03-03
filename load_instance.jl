using HDF5

struct Instance
    clients
    stations
    Γ
    ϕ
    Σ
    Π
    k
    c
    T
end

function load_istance(file_name, instance)
    h5open(file_name) do file
        inst = file[instance]
        clients = read(inst["clients"])
        stations = read(inst["stations"])
        Γ = read(inst["gamma"])
        ϕ = read(inst["phi"])
        Σ = read(inst["sigma"])
        Π = read(inst["pi"])
        k = read(inst["k"])
        c = read(inst["c"])
        T = read(inst["T"])
        Instance(clients, stations, Γ, ϕ, Σ, Π, k, c, T)
    end
end


