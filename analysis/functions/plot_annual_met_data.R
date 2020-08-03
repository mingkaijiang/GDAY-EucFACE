plot_annual_met_data <- function() {
    
    
    ### read in DF
    spinDF <- read.csv("met_data/EUC_met_spinup_daily_50yrs.csv",skip=4)
    histDF <- read.csv("met_data/EUC_met_historic_daily_1750_2011.csv",skip=4)
    obsDF <- read.csv("met_data/EUC_met_DRY_AMB_daily_2012_2019.csv",skip=4)
    futDF <- read.csv("met_data/EUC_met_DRY_AMB_NOP_daily_2020_2069.csv",skip=4)
    
    ### colnames
    colnames(spinDF)[1] <- colnames(histDF)[1] <- colnames(obsDF)[1] <- colnames(futDF)[1] <- "year"
    spinDF$Trt <- "Spinup"
    histDF$Trt <- "Hist"
    obsDF$Trt <- "Obs"
    futDF$Trt <- "Fut"
    
    ### merge
    myDF <- rbind(spinDF, histDF, obsDF, futDF)
    
    ### annual
    plotDF <- summaryBy(.~year+Trt, data=myDF, FUN=mean, keep.names=T, na.rm=T)
    
    n <- dim(plotDF)[2]
    
    
    pdf("outputs/analysis/annual_met_data_validation.pdf", width = 8, height=6)
    for (i in 4:n) {
        p <- ggplot(plotDF) +
            geom_point(aes(x = year, y = plotDF[,i], fill = Trt, pch = Trt), size=4)+
            geom_line(aes(x = year, y = plotDF[,i], col=Trt))+
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
            ylab(colnames(plotDF)[i])+
            scale_color_manual(name="",
                               limits=c("Spinup", "Hist", "Obs", "Fut"),
                               labels=c("Spinup", "Hist", "Obs", "Fut"),
                               values=c("grey", "black", "green", "purple"),
                               guide=guide_legend(nrow=1))+
            scale_fill_manual(name="",
                              limits=c("Spinup", "Hist", "Obs", "Fut"),
                              labels=c("Spinup", "Hist", "Obs", "Fut"),
                              values=c("grey", "black", "green", "purple"),
                              guide=guide_legend(nrow=1))+
            scale_linetype_manual(name="",
                                  limits=c("Spinup", "Hist", "Obs", "Fut"),
                                  labels=c("Spinup", "Hist", "Obs", "Fut"),
                                  values=c("dotted", "dotted", "solid", "solid"),
                                  guide=guide_legend(nrow=1))+
            scale_shape_manual(name="",
                               limits=c("Spinup", "Hist", "Obs", "Fut"),
                               labels=c("Spinup", "Hist", "Obs", "Fut"),
                               values=c(24,24,21,21),
                               guide=guide_legend(nrow=1))+
            ggtitle(colnames(plotDF)[i])+
            xlab("Year")+
            scale_x_continuous(limits=c(1700, 2070),
                               breaks=c(1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050))
        
        
        plot(p)
    }
    
    dev.off()
    
}