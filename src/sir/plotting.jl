function plotpath(X; name="path") 
    n_particles = length(X[1])
    Xc = copy(X)
    push!(Xc, vcat([_S_, _I_, _R_], fill(_L_, n_particles-3)))
    # construct observation ColorPalette
#    defaultpalette = palette(cgrad(:default, categorical=true), 3)
    # white = RGBA{Float64}(255, 255, 255)
    # white = RGBA{Float64}(16, 59, 223, 0.12)
    white = RGBA(52, 162, 231, 0.23)
#    observationcolors = vec(hcat(white, defaultpalette.colors.colors...))
#    observationpalette = ColorPalette(typeof(defaultpalette.colors)(observationcolors, "", ""))

    p = heatmap(obs2matrix(Xc)', xlabel="time", ylabel="individual", 
                colorbar=true, 
                #color=observationpalette, 
                dps=600, title=name, background_color_subplot=white)
    return p
end
