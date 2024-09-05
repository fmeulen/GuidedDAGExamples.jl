
"""
    transpose_state_matrix(M::Matrix{State})

    transpose matrix of type State
"""
function transpose_state_matrix(M::Matrix{State})
    # Get the size of the matrix
    rows, cols = size(M)
    
    # Create an empty matrix with transposed dimensions
    transposed_matrix = Matrix{State}(undef, cols, rows)
    
    # Loop over the elements and assign them to the transposed positions
    for i in 1:rows
        for j in 1:cols
            transposed_matrix[j, i] = M[i, j]
        end
    end
    
    return transposed_matrix
end

"""
    plotpath(X;  name="path") 

    Plot path of simulated interacting particles.
    
    X should be Vector{Vector{State}} (alias for Array{Array{State, 1}, 1})
    
    Example 
    X = [
        [_S_, _S_, _S_, _I_],
        [_L_, _R_, _I_, _R_],
        [_R_, _L_, _S_, _I_]];
    
    plotpath(X)
"""
function plotpath(X;  name="path") 
    M = [X[j][i] for j in eachindex(X), i in eachindex(X[1])]
    M = transpose_state_matrix(M)

    # Create a dictionary to map states to numeric values
    color_map = Dict(
        _S_ => 1,  # Yellow
        _I_ => 2,  # Red
        _R_ => 3,  # Green
        _L_ => 4   # Black
    )
    
    # Convert the matrix to a numeric matrix using the color map
    numeric_matrix = [color_map[state] for state in M]

    # Define the color palette corresponding to the numeric values
    # The colors are assigned in the order corresponding to 1, 2, 3, and 4
    color_palette = [:yellow, :red, :green, :black]

    # Plot the heatmap using the numeric matrix and the custom color palette
    p = heatmap(numeric_matrix, color=color_palette,  clim=(1, 4),title=name,  xlabel="time", ylabel="individual")
    p 
end
 

