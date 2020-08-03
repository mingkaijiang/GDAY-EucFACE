plot_daily_met_data <- function() {
    
    
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
    
    ### merge spinup and historic
    plotDF1 <- rbind(spinDF, histDF)
    plotDF1 <- subset(plotDF1, year <= 1755 & year >= 1745)
    
    ### merge historic and observed and future
    plotDF2 <- rbind(histDF, obsDF, futDF)
    plotDF2 <- subset(plotDF2, year <= 2025 & year >= 2006)
    
    ### add Date
    plotDF1$Date <- as.Date((plotDF1$doy-1), origin = paste0(plotDF1$year, "-01-01"))
    plotDF2$Date <- as.Date((plotDF2$doy-1), origin = paste0(plotDF2$year, "-01-01"))
    
    
    ### variable sample size
    n <- dim(plotDF1)[2]
    
    
    pdf("outputs/analysis/daily_met_spinup_historic_data_validation.pdf", width = 8, height=6)
    for (i in 4:(n-1)) {
        p <- ggplot(plotDF1) +
            geom_point(aes(x = Date, y = plotDF1[,i], fill = Trt, pch = Trt), size=4)+
            geom_line(aes(x = Date, y = plotDF1[,i], col=Trt))+
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
            ylab(colnames(plotDF1)[i])+
            scale_color_manual(name="",
                               limits=c("Spinup", "Hist"),
                               labels=c("Spinup", "Hist"),
                               values=c("grey", "black"),
                               guide=guide_legend(nrow=1))+
            scale_fill_manual(name="",
                              limits=c("Spinup", "Hist"),
                              labels=c("Spinup", "Hist"),
                              values=c("grey", "black"),
                              guide=guide_legend(nrow=1))+
            scale_linetype_manual(name="",
                                  limits=c("Spinup", "Hist"),
                                  labels=c("Spinup", "Hist"),
                                  values=c("dotted", "dotted"),
                                  guide=guide_legend(nrow=1))+
            scale_shape_manual(name="",
                               limits=c("Spinup", "Hist"),
                               labels=c("Spinup", "Hist"),
                               values=c(24,24),
                               guide=guide_legend(nrow=1))+
            ggtitle(colnames(plotDF1)[i])+
            xlab("Year")
        
        
        plot(p)
    }
    
    dev.off()
    
    
    
    pdf("outputs/analysis/daily_met_historic_obs_future_data_validation.pdf", width = 8, height=6)
    for (i in 4:(n-1)) {
        p <- ggplot(plotDF2) +
            geom_point(aes(x = Date, y = plotDF2[,i], fill = Trt, pch = Trt), size=4)+
            geom_line(aes(x = Date, y = plotDF2[,i], col=Trt))+
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
            ylab(colnames(plotDF2)[i])+
            scale_color_manual(name="",
                               limits=c("Hist", "Obs", "Fut"),
                               labels=c("Hist", "Obs", "Fut"),
                               values=c("black", "green", "purple"),
                               guide=guide_legend(nrow=1))+
            scale_fill_manual(name="",
                              limits=c("Hist", "Obs", "Fut"),
                              labels=c("Hist", "Obs", "Fut"),
                              values=c("black", "green", "purple"),
                              guide=guide_legend(nrow=1))+
            scale_linetype_manual(name="",
                                  limits=c("Hist", "Obs", "Fut"),
                                  labels=c("Hist", "Obs", "Fut"),
                                  values=c( "dotted", "solid", "solid"),
                                  guide=guide_legend(nrow=1))+
            scale_shape_manual(name="",
                               limits=c("Hist", "Obs", "Fut"),
                               labels=c("Hist", "Obs", "Fut"),
                               values=c(24,21,21),
                               guide=guide_legend(nrow=1))+
            ggtitle(colnames(plotDF2)[i])+
            xlab("Year")
        
        
        plot(p)
    }
    
    dev.off()
    
}