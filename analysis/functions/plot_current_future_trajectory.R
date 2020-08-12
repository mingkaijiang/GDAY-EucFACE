plot_current_future_trajectory <- function() {
    
    ### obs
    obsDF1 <- read.csv("outputs/EUC_simulated_DRY_AMB_2012_2019.csv", skip=1)
    obsDF2 <- read.csv("outputs/EUC_simulated_WET_AMB_2012_2019.csv", skip=1)

    ### future
    futDF1 <- read.csv("outputs/EUC_simulated_DRY_AMB_NOP_2020_2069.csv", skip=1)
    futDF2 <- read.csv("outputs/EUC_simulated_DRY_AMB_MDP_2020_2069.csv", skip=1)
    futDF3 <- read.csv("outputs/EUC_simulated_DRY_AMB_HIP_2020_2069.csv", skip=1)
    
    futDF4 <- read.csv("outputs/EUC_simulated_WET_AMB_NOP_2020_2069.csv", skip=1)
    futDF5 <- read.csv("outputs/EUC_simulated_WET_AMB_MDP_2020_2069.csv", skip=1)
    futDF6 <- read.csv("outputs/EUC_simulated_WET_AMB_HIP_2020_2069.csv", skip=1)
    
    ### assign factor
    obsDF1$Trt <- "DRY_AMB_2012_2019"
    obsDF2$Trt <- "WET_AMB_2012_2019"
    
    futDF1$Trt <- "DRY_AMB_NOP_2020_2069"
    futDF2$Trt <- "DRY_AMB_MDP_2020_2069"
    futDF3$Trt <- "DRY_AMB_HIP_2020_2069"
    futDF4$Trt <- "WET_AMB_NOP_2020_2069"
    futDF5$Trt <- "WET_AMB_MDP_2020_2069"
    futDF6$Trt <- "WET_AMB_HIP_2020_2069"
    
    ### merge
    myDF <- rbind(obsDF1, obsDF2, futDF1, futDF2, futDF3, futDF4, futDF5, futDF6)
    
    ### calculate annual mean and/or sum
    sumDF1 <- summaryBy(shoot+lai+branch+stem+root+croot+shootn+branchn+stemn+rootn+crootn+
                            shootp+branchp+stemp+rootp+crootp+soilc+soiln+soilp+inorgn+inorgp+
                            inorgavlp+inorglabp+inorgsorbp+inorgssorbp+inorgoccp+inorgparp+fertilizerp+
                            litterc+littercag+littercbg+litternag+litternbg+litterpag+litterpbg+
                            activesoil+slowsoil+passivesoil+activesoiln+slowsoiln+passivesoiln+
                            activesoilp+slowsoilp+passivesoilp~year+Trt, FUN=mean,
                        data=myDF, na.rm=T, keep.names=T)
    
    
    sumDF2 <- summaryBy(et+transpiration+soil_evap+canopy_evap+runoff+nep+gpp+npp+hetero_resp+
                            auto_resp+cpleaf+cpbranch+cpstem+cproot+
                            nuptake+ngross+nmineralisation+nloss+puptake+pgross+
                            pmineralisation+ploss+p_slow_biochemical+leafretransn+leafretransp~year+Trt,
                        FUN=sum, data=myDF, na.rm=T, keep.names=T)
    
    ### merge
    plotDF <- merge(sumDF1, sumDF2, by=c("year", "Trt"))
    
    
    
    ## number of columns
    n <- dim(plotDF)[2]
    
    
    pdf("outputs/analysis/current_future_trajectory.pdf")
    for (i in 3:n) {
    p <- ggplot(plotDF) +
        geom_point(aes(x = year, y = plotDF[,i]*100, fill = Trt, pch = Trt), size=4)+
        geom_line(aes(x = year, y = plotDF[,i]*100, col=Trt))+
        theme_linedraw() +
        theme(panel.grid.minor=element_blank(),
              axis.text.x=element_text(size=12),
              axis.title.x=element_blank(),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_text(size=16),
              panel.grid.major=element_blank(),
              legend.position="bottom",
              legend.box = 'horizontal',
              legend.box.just = 'left',
              plot.title = element_text(size=16, face="bold.italic", 
                                        hjust = 0.5))+
        ylab(paste0(colnames(plotDF)[i], " [g m-2 (yr-1 for flux)]"))+
        scale_color_manual(name="",
                           limits=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                    "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                    "DRY_AMB_HIP_2020_2069",
                                    "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                    "WET_AMB_HIP_2020_2069"),
                           labels=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                    "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                    "DRY_AMB_HIP_2020_2069",
                                    "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                    "WET_AMB_HIP_2020_2069"),
                           values=c("yellow", "cyan", "orange", "red2", "brown",
                                    "deepskyblue1", "blue", "darkblue"),
                           guide=guide_legend(nrow=4))+
        scale_fill_manual(name="",
                          limits=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                   "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                   "DRY_AMB_HIP_2020_2069",
                                   "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                   "WET_AMB_HIP_2020_2069"),
                          labels=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                   "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                   "DRY_AMB_HIP_2020_2069",
                                   "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                   "WET_AMB_HIP_2020_2069"),
                          values=c("yellow", "cyan", "orange", "red2", "brown",
                                   "deepskyblue1", "blue", "darkblue"),
                          guide=guide_legend(nrow=4))+
        scale_linetype_manual(name="",
                              limits=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                       "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                       "DRY_AMB_HIP_2020_2069",
                                       "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                       "WET_AMB_HIP_2020_2069"),
                              labels=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                       "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                       "DRY_AMB_HIP_2020_2069",
                                       "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                       "WET_AMB_HIP_2020_2069"),
                              values=c("dotted", "dotted", "solid", "solid",
                                       "solid", "solid", "solid", "solid"),
                              guide=guide_legend(nrow=4))+
        scale_shape_manual(name="",
                           limits=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                    "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                    "DRY_AMB_HIP_2020_2069",
                                    "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                    "WET_AMB_HIP_2020_2069"),
                           labels=c("DRY_AMB_2012_2019", "WET_AMB_2012_2019",
                                    "DRY_AMB_NOP_2020_2069", "DRY_AMB_MDP_2020_2069", 
                                    "DRY_AMB_HIP_2020_2069",
                                    "WET_AMB_NOP_2020_2069", "WET_AMB_MDP_2020_2069", 
                                    "WET_AMB_HIP_2020_2069"),
                           values=c(24,24,21,21, 21, 21, 21, 21),
                           guide=guide_legend(nrow=4))+
        ggtitle(colnames(plotDF)[i])+
        xlab("Year")+
        scale_x_continuous(limits=c(2010, 2070),
                           breaks=c(2010, 2020, 2030, 2040, 2050, 2060, 2070))
    
    
    plot(p)
    }
    
    dev.off()
    
 
}