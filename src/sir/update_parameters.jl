


#     μᵒ = P.μ  * exp(propσ*randn())
#     Pᵒ = SIRguided(P.λ, μᵒ, P.ν, P.τ, P.𝒩)
#     for k in 1:nseg
#         Xᵒ[k], lrᵒ[k] = sample_segment!(Pᵒ, Xᵒ[k], Xobs[k],Xobs[k+1],Z[k],Q[k],J)
#     end
#     if log(rand()) < sum(lrᵒ .- lr)  + (log(Pᵒ.μ) - log(P.μ))+ logpdf(prior[2],μᵒ) - logpdf(prior[2],P.μ)
#         if mod(it, skip_print)==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    μ ", μᵒ)  end
#         X, Xᵒ, lr, lrᵒ, P, Pᵒ = Xᵒ, X, lrᵒ, lr, Pᵒ, P
#         accpar[2] += 1
#     else
#         if mod(it, skip_print )==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    μ ", P.μ)  end
#     end

#     νᵒ = P.ν  * exp(propσ*randn())
#     Pᵒ = SIRguided(P.λ, P.μ, νᵒ, P.τ, P.𝒩)
#     for k in 1:nseg
#         Xᵒ[k], lrᵒ[k] = sample_segment!(Pᵒ, Xᵒ[k], Xobs[k],Xobs[k+1],Z[k],Q[k],J)
#     end
#     if log(rand()) < sum(lrᵒ .- lr)  + (log(Pᵒ.ν) - log(P.ν)) + logpdf(prior[3],νᵒ) - logpdf(prior[3],P.ν)
#         if mod(it, skip_print)==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    ν ", νᵒ)  end
#         X, Xᵒ, lr, lrᵒ, P, Pᵒ = Xᵒ, X, lrᵒ, lr, Pᵒ, P
#         accpar[3] += 1
#     else
#         if mod(it, skip_print )==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    ν ", P.ν)  end
#     end


#     return X, Xᵒ, lr, lrᵒ, P, Pᵒ

# end
