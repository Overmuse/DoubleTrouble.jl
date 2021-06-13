function index_to_pair(idx)
    j = (1+isqrt(8idx-7))÷2+1
    i = idx - ((j-2)*(j-1))÷2
    (i, j)
end

function spreads(data)
    m, n = size(data)
    out = zeros(m, n*(n-1)÷2)
    x = 1
    for i in 2:size(data, 2), j in 1:i-1
        out[:, x] = log.(data[:, j] ./ data[1, j]) .- log.(data[:, i] ./ data[1, i])
        x += 1
    end
    out
end 

