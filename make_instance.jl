
function make_instance(client_matrix, station_matrix; p_dist=2, threashold=1e-3)
    A, Λ = init_alpha_matrix(client_matrix, station_matrix)
    
    for (j, (x, y, v, d_csr, d_cli)) in enumerate(eachrow(client_matrix))
        p_csr, p_cli = make_probability(x, y, d_csr, d_cli, p_dist)
        for (i, (sx, sy)) in enumerate(eachrow(station_matrix))
            A[i, j] = p_cli(sx, sy) * v 
            Λ[i] += p_csr(sx, sy) * v
        end
    end
    A[A .<= threashold] .= 0
    A, Λ 
end

function make_probability(i, j, d_csr, d_client, p_dist) 
    d_csr /= log(p_dist)
    d_client /= log(p_dist)
    p_csr = (x, y) -> exp(-dist(i, j, x, y) / d_csr)
    p_client = (x, y) -> exp(-dist(i, j, x, y) / d_client)
    p_csr, p_client
end


function init_alpha_matrix(cm, sm) 
    x = size(cm, 1)
    y = size(sm, 1)
    zeros(y, x), zeros(y)
end

dist = (x₁, y₁, x₂, y₂) -> √((x₁ - x₂)^2 + (y₁ - y₂)^2)
