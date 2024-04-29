function partition_into_blocks_close(N::Int, n::Int)
    # Check if n is smaller than N
    if n >= N
        throw(ArgumentError("n must be smaller than N"))
    end
    
    # Initialize partitions array
    partitions = Array{Vector{Int}, 1}(undef, n)
    
    # Initialize blocks
    for i in 1:n
        partitions[i] = Int[]
    end
    
    # Calculate the number of elements per block
    elements_per_block = div(N, n)
    
    # Distribute the elements evenly
    current_block = 1
    current_element = 1
    for _ in 1:n
        while length(partitions[current_block]) < elements_per_block && current_element <= N
            push!(partitions[current_block], current_element)
            current_element += 1
        end
        current_block += 1
    end
    
    # Distribute the remaining elements
    while current_element <= N
        push!(partitions[mod(current_element, n) == 0 ? n : mod(current_element, n)], current_element)
        current_element += 1
    end
    
    return partitions
end

