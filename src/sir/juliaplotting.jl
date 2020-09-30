

if PLOT
    # traceplots
    trλ = Plots.plot(1:ITER, map(x->x.λ, θθ ))
    trμ = Plots.plot(1:ITER, map(x->x.μ, θθ ))
    trν = Plots.plot(1:ITER, map(x->x.ν, θθ ))
    Plots.plot(trλ, trμ, trν)#, layout = @layout [a;b;c])



    if PLOT
    l = @layout [a b ; c d]
    pheatforw = Plots.heatmap(obs2matrix(Xtrue),title="forward")
    pheatguidedinit = Plots.heatmap(obs2matrix(Xinit),title="guided init")
    pheatguidedfinal = Plots.heatmap(-obs2matrix(Xfinal),title="guided final")

    u = -1.0*obs2matrix(Xtrue)
    for j in axes(u, 1), k in axes(u, 2)
        if j ∈ 2:J:J*nobs || j ∈ 3:J:J*nobs || j ∈ 4:J:J*nobs
            u[j,k] = u[j-1,k]
        elseif j ∉ 1:J:J*nobs
            u[j,k] = NaN
        end
    end
    pobs = heatmap(u, title="observed")

    pp = Plots.plot(pobs, pheatguidedinit, pheatforw, pheatguidedfinal, layout=l, size=(1000,1200))
    png(pp, "heatmaps.png")





















    u = 1.0*obs2matrix(Xtrue)
for j in axes(u, 1), k in axes(u, 2)
    if j ∉ 1:J:J*nobs
        u[j,k] = NaN
    end
    if j ∈ 2:J:J*nobs || j ∈ 3:J:J*nobs
        u[j,k] = u[j-1,k]
    end
end
pobs = heatmap(u, title="observed")

    # configuration plots
    out = []
    yax = []
    u = obs2matrix(Xobs)
    u0 = 0.0 * u[1,:] .- .5
    for i in 1:size(u)[1]
        push!(out, u[i,:])
        push!(yax, i+(i-1)*J)
        push!(out, u0, u0, u0,u0)

        push!(yax, i+(i-0.8)*J, i+(i-0.6)*J, i+(i-0.4)*J, i+(i-0.2)*J)
    end
    ttt= [out[j][i] for j in eachindex(out), i in eachindex(out[1])]
    pheatobs  = Plots.heatmap(1:n, yax, ttt,title="observed")
    Xfinal= collect(Iterators.flatten(X))
    l = @layout [a b ; c d]
    pheatforw = Plots.heatmap(obs2matrix(Xtrue[1:(nobs-1)*J]),title="forward")
    pheatguidedinit = Plots.heatmap(obs2matrix(Xinit),title="guided init")
    pheatguidedfinal = Plots.heatmap(obs2matrix(Xfinal),title="guided final")

    pp = Plots.plot(pheatforw, pheatobs, pheatguidedinit, pheatguidedfinal)#, layout=@layout [a b ; c d])
    png(pp, "heatmaps.png")

end



l = @layout [a b ; c d]
#pheatobs = Plots.heatmap(obs2matrix(Xobs),title="observed")
pheatforw = Plots.heatmap(-obs2matrix(Xtrue),title="forward")
pheatguidedinit = Plots.heatmap(-obs2matrix(Xinit),title="guided init")
pheatguidedfinal = Plots.heatmap(-obs2matrix(Xfinal),title="guided final")

u = -1.0*obs2matrix(Xtrue)
for j in axes(u, 1), k in axes(u, 2)
    if j ∈ 2:J:J*nobs || j ∈ 3:J:J*nobs || j ∈ 4:J:J*nobs
        u[j,k] = u[j-1,k]
    elseif j ∉ 1:J:J*nobs
        u[j,k] = NaN
    end
end
pobs = heatmap(u, title="observed")

pp = Plots.plot(pobs, pheatguidedinit, pheatforw, pheatguidedfinal, layout=l, size=(1000,1200))
png(pp, "heatmaps.png")
