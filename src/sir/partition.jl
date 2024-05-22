

function make_partition(N::Int, n_blocks::Int)  # will results in possibly one more block than n_blocks
    if n_blocks > N
        throw(ArgumentError("n must be smaller than N"))
    end
    k = Int(floor(N/n_blocks))
    res = [i*k+1:(i+1)*k for i in 0:n_blocks-1]
    if N> n_blocks*k
        push!(res, n_blocks*k+1 : N)
    end
    res
end

#make_partition(12,1)
