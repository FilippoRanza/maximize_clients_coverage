

function make_instance(client_matrix, station_list; threashold=1e-3, scale=250)
    A, idx = init_output(client_matrix, station_list)
    scale /= log(2)
    for (j, (sx, sy)) in enumerate(eachrow(station_list))
        p = (x, y) -> exp(-dist(x, y, sx, sy) / scale)
        for (i, (cx, cy, _)) in enumerate(idx) 
            A[i, j] = p(cx, cy) 
        end
    end
    A[A .<= threashold] .= 0
    A, idx
end

function init_output(client_matrix, station_list)
    client_count = length(client_matrix[client_matrix .> 0])
    station_count = size(station_list, 1)
    A = zeros(client_count, station_count)
    idx = non_zero_list(client_matrix)
    A, idx
end

function non_zero_list(M)
    x, y = size(M)
    [(i, j, M[i, j]) for i ∈ 1:x for j ∈ 1:y if M[i, j] > 0]
end


dist = (x₁, y₁, x₂, y₂) -> √((x₁ - x₂)^2 + (y₁ - y₂)^2)

