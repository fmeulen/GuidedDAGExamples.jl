proc(x) =vcat(obs2matrix(x)...)


xgr = repeat(1:n,inner=length(tobs))
ygr = repeat(tobs, n)

ddobs = DataFrame(x=xgr,y=ygr, z= proc(Xobs))

################
Xf = [Xfinal[1]]
Xi = [Xinit[1]]
Xm = [Xmid[1]]
for i ∈ 2:length(Xfinal)
    push!(Xf,Xfinal[i][2:end])
    push!(Xi,Xinit[i][2:end])
    push!(Xm,Xmid[i][2:end])
end
Xfinal_fl= collect(Iterators.flatten(Xf))
Xinit_fl= collect(Iterators.flatten(Xi))
Xmid_fl= collect(Iterators.flatten(Xm))

xgr = repeat(1:n,inner=length(tfull))
ygr = repeat(tfull, n)

Xtrueobs = deepcopy(Xtrue)
for i ∈ eachindex(Xtrueobs)
    if !(i ∈ collect(tobs) )
        Xtrueobs[i] .= fill(State(7),n)
    end
end

ddsim = DataFrame(x=xgr,y=ygr, zinit= proc(Xinit_fl), zmid = proc(Xmid_fl),
 zfinal= proc(Xfinal_fl), ztrue = proc(Xtrue),ztrueobs = proc(Xtrueobs))


@rput ddobs
@rput ddsim
@rput figdir

R"""
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())
colors <- c("khaki2", "indianred3", "black")
colors <- c("khaki2", "red2", "black")


pobs <- ddobs %>% mutate(type=as.factor(z)) %>%
    mutate(type=fct_recode(type, S="1", I= "2", R="3"))%>%
     ggplot(aes(x=x,y=y,fill=type)) + geom_raster() +
    scale_fill_manual(values=colors)+xlab("") + ylab("")+
    scale_x_continuous(breaks=seq(10,100,by=10)) + scale_y_continuous(breaks=unique(ddobs$y))+
    ggtitle("observed data interpolated") + theme(legend.position="none")

#pdf(paste0(figdir,"/obs.pdf"),width=4,height=3)
#    show(pobs)
#dev.off()

pinit <- ddsim %>% mutate(type=as.factor(zinit)) %>%
    mutate(type=fct_recode(type, S="1", I= "2", R="3"))%>%
     ggplot(aes(x=x,y=y,fill=type)) + geom_raster() +
    scale_fill_manual(values=colors)+xlab("") + ylab("")+#xlab("individual") + ylab("time")+
    scale_x_continuous(breaks=seq(10,100,by=10)) + scale_y_continuous(breaks=unique(ddobs$y))+
    ggtitle("initial iterate")+ theme(legend.position="none")

pmid <- ddsim %>% mutate(type=as.factor(zmid)) %>%
    mutate(type=fct_recode(type, S="1", I= "2", R="3"))%>%
     ggplot(aes(x=x,y=y,fill=type)) + geom_raster() +
    scale_fill_manual(values=colors)+xlab("") + ylab("")+#xlab("individual") + ylab("time")+
    scale_x_continuous(breaks=seq(10,100,by=10)) + scale_y_continuous(breaks=unique(ddobs$y))+
    ggtitle("middle iterate")+ theme(legend.position="none")


pfinal <- ddsim %>% mutate(type=as.factor(zfinal)) %>%
    mutate(type=fct_recode(type, S="1", I= "2", R="3"))%>%
     ggplot(aes(x=x,y=y,fill=type)) + geom_raster() +
    scale_fill_manual(values=colors)+xlab("") + ylab("")+#xlab("individual") + ylab("time")+
    scale_x_continuous(breaks=seq(10,100,by=10)) + scale_y_continuous(breaks=unique(ddobs$y))+
    ggtitle("final iterate")+ theme(legend.position="none")

ptrue <- ddsim %>% mutate(type=as.factor(ztrue)) %>%
    mutate(type=fct_recode(type, S="1", I= "2", R="3"))%>%
     ggplot(aes(x=x,y=y,fill=type)) + geom_raster() +
    scale_fill_manual(values=colors)+xlab("") + ylab("")+#xlab("individual") + ylab("time")+
    scale_x_continuous(breaks=seq(10,100,by=10)) + scale_y_continuous(breaks=unique(ddobs$y))+
    ggtitle("true forward simulated")+ theme(legend.position="none")

colors <- c("khaki2", "indianred3", "black", "white")
colors <- c("khaki2", "indianred3", "black", "azure1")
colors <- c("khaki2", "red2", "black", "azure1")
ptrueobs <- ddsim %>% mutate(type=as.factor(ztrueobs)) %>%
    mutate(type=fct_recode(type, S="1", I= "2", R="3", Latent="7"))%>%
     ggplot(aes(x=x,y=y,fill=type)) + geom_raster() +
    scale_fill_manual(values=colors)+xlab("") + ylab("")+#xlab("individual") + ylab("time")+
    scale_x_continuous(breaks=seq(20,100,by=20)) + scale_y_continuous(breaks=unique(ddobs$y))+
    ggtitle("observed data")



library(gridExtra)
pdf(paste0(figdir,"/all.pdf"),width=6,height=7, pointsize=10)
#grid.arrange(ptrueobs, pobs, pinit, pmid, pfinal,ptrue)
grid.arrange(ptrueobs, pinit, pmid, pfinal,ptrue, ncol=2, layout_matrix=cbind(c(1,2,4),c(1,3,5)))
dev.off()


#png(paste0(figdir,"/all.png"),bg = "transparent", res=101)#,width=6,height=7, pointsize=10)
#grid.arrange(ptrueobs, pobs, pinit, pmid, pfinal,ptrue)
#dev.off()

"""

# trace plots: to be done
ec(x,i) = map(u->u[i],x)
dpar = DataFrame(iterate=repeat(1:ITER,3),
        parameter=vcat(map(x->x.λ, θθ),map(x->x.μ, θθ),map(x->x.ν, θθ)),
        type = repeat(["lambda", "mu", "nu"], inner = ITER))


dtrue = DataFrame(type=["lambda", "mu", "nu"], parameter=[Ptrue.λ, Ptrue.μ, Ptrue.ν])
@rput dpar
@rput dtrue
R"""
dpar$type <- factor(dpar$type, labels=c("lambda","mu","nu"))
p <- ggplot() + geom_path(data=dpar, mapping= aes(x=iterate,colour=type,y=parameter),size=1.05) +
 geom_hline(data=dtrue, mapping=aes(yintercept=parameter)) +
 facet_grid(type~.,labeller=label_parsed,scales='free') + theme(legend.position="none") + ylab("")+
 theme(strip.text.y = element_text(size=12,angle=0))
pdf(paste0(figdir,"/traceplots.pdf"), width = 8, height=5)
show(p)
dev.off()

png(paste0(figdir,"/traceplots.png"))
show(p)
dev.off()

"""

# neighbours in aux process plot
ppaux = heatmap([N[i[2]][i[1]][k] for i in vec(CartesianIndices((J, nobs-1))), k in 1:n])
png(ppaux, joinpath(figdir,"infected_neighbours_aux.png"))
