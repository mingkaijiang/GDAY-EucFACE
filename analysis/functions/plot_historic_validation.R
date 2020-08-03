plot_historic_validation <- function() {
    ### purpose: 
    ### plot 1750 - to 2069 to just smoothness of the re-start breaks,
    ### and the magnitude of some key response variables
    
    ### historic df
    histDF <- read.csv("outputs/EUC_amb_equilib.csv",skip=1)
    
    ### obs
    obsDF <- read.csv("outputs/EUC_simulated_DRY_AMB_2012_2019.csv", skip=1)
    
    ### future
    futDF <- read.csv("outputs/EUC_simulated_DRY_AMB_NOP_2020_2069.csv", skip=1)
    
    
    ### merge
    myDF <- rbind(histDF, obsDF, futDF)
    
    ### calculate annual mean and/or sum
    sumDF1 <- summaryBy(shoot+lai+branch+stem+root+croot+shootn+branchn+stemn+rootn+crootn+
                            shootp+branchp+stemp+rootp+crootp+soilc+soiln+soilp+inorgn+inorgp+
                            inorgavlp+inorglabp+inorgsorbp+inorgssorbp+inorgoccp+inorgparp+
                            litterc+littercag+littercbg+litternag+litternbg+litterpag+litterpbg+
                            activesoil+slowsoil+passivesoil+activesoiln+slowsoiln+passivesoiln+
                            activesoilp+slowsoilp+passivesoilp~year, FUN=mean,
                        data=myDF, na.rm=T, keep.names=T)
    
    
    sumDF2 <- summaryBy(et+transpiration+soil_evap+canopy_evap+runoff+nep+gpp+npp+hetero_resp+
                            auto_resp+nuptake+ngross+nmineralisation+nloss+puptake+pgross+
                            pmineralisation+ploss+leafretransn+leafretransp~year,
                        FUN=sum, data=myDF, na.rm=T, keep.names=T)
    
    ### merge
    plotDF <- merge(sumDF1, sumDF2, by="year")
    
    ## number of columns
    n <- dim(plotDF)[2]
    
    
    ### plotting
    pdf("outputs/analysis/historic_validation.pdf")
    for (i in 2:n) {
        plot(plotDF[,i]~plotDF$year, xlab="year", ylab=colnames(plotDF)[i])
        title(colnames(plotDF)[i])
    }
    dev.off()
}