time0 = proc.time()
###----- 1. Set parameters
project = "Heart progerin" # mir29   Heart progerin   Heart lamin  ISO Challenge  PCTX Treatment  Zmpste-Rankl  HGPS Amanda DBU Alberto
technique = "ECG" # Telemetry ECG Thermochallenge.ECG ECG.extra
processall = F # Default is F
frequency = 2000  # Acquisition frequency in Hz 2000 5000
segment = F # F for ECG
lengthsegment = 5 # Duration in min of each segment 1 5
normalize = T # Default is T
noisespan = 3 # Span of noise correction  3
basespan = 100 # Span for baseline calculation  100

###----- 2. Import data
## Load libraries
packages = c("data.table", "matrixStats")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) install.packages(setdiff(packages, rownames(installed.packages())))  
for (i in packages) library(i, character.only = T)


## Set directories
if (.Platform$OS.type == "unix") setwd("/Volumes/Victor/") else setwd("S:/LAB_VA/LAB/Victor/")
# setwd("/Users/Victor/Downloads/")
# setwd("S:/LAB_VA/LAB/Alvaro Macias/Results/Esperanza de vida/Challenge de Isoproterenol/Electros/")
# setwd("S:/LAB_VA/LAB/Alvaro Macias/Results/Esperanza de vida/Paclitaxel 0.5mg_kg 3_semana/ECGs/")
# setwd("/Volumes/LAB/Alberto M/Alberto/Ratones endotelio y SMC/ECG DBU/")
if(project == "PCTX Treatment") setwd("/Volumes/LAB_VA/LAB/Alvaro Macias/Results/Tratamientos/Paclitaxel 0.5mg_kg 3_semana/ECGs sucesivos tras PCTX agudo - REAL/")
route = paste0(project, " project/Raw data/", technique, "/")
origin = "raw txt/"
destination = ifelse(normalize, "normalized/", "not normalized/")
dir.create(paste0(route, destination), showWarnings = F)

## List of files to process
registers = sapply(dir(paste0(route, origin)), function (x) list(gsub(".txt", "", dir(paste0(route, origin, x))[grep(".txt", dir(paste0(route, origin, x)))])))
ecglist = c()
for (i in names(registers)) if(length(registers[[i]]) != 0) ecglist = rbind.data.frame(ecglist, cbind(i,registers[[i]], paste(i,gsub(" ", "", toupper(registers[[i]])))), stringsAsFactors = F)
if (!processall) {
  ecglista = c()
  if (file.exists(paste0(route, destination, "processed.log.txt"))) ecglista = data.frame(fread(paste0(route, destination, "processed.log.txt"), encoding = "UTF-8", header = T, stringsAsFactors = F, dec = ".", na.strings = c("N/A", "NA", "nan"), sep = "\t"))[,1:3]
  ecglist = merge(ecglist, setdiff(as.vector(ecglist[,3]), as.vector(ecglista[,3])), by.x = 3, by.y = 1)
  ecglist = ecglist[,c(2,3,1)]
}
# ecglist = ecglist[grepl("2016-10",ecglist$i),]

## Read files
if (file.exists(paste0(route, "correct.polarity.txt"))) {
  correctpol = data.frame(fread(paste0(route, "correct.polarity.txt"), encoding = "UTF-8", header = T, stringsAsFactors = F, dec = ".", na.strings = c("N/A", "NA", "nan"), sep = "\t"))
} else rm(correctpol)
Date = unique(ecglist[,1])[1]
Id = ecglist[ecglist[,1] == Date, 2][1]
for (Date in unique(ecglist[,1])) {
  dir.create(paste0(route, destination, Date), showWarnings = F)
  for (Id in ecglist[ecglist[,1] == Date, 2]) {
    tryCatch({
      rawdata = data.frame(fread(paste0(route, origin, Date, "/", Id, ".txt"), encoding = "UTF-8", header = F, stringsAsFactors = F, dec = ".", na.strings = c("N/A", "NA", "nan"), sep = ifelse(technique == "Telemetry", ",", "\t"), skip = 20, data.table = F))
      if (technique == "Telemetry") {
        names(rawdata) = c("Time", "CH2")
      } else rawdata = rawdata[,1:4] # Remove calculated leads

      
      ###----- 3. Correct frequency and polarity
      ## Set acquisition frequency to 2000 Hz
      if (frequency == 5000) {
        changefreq = function (x,y) rawdata[seq(x, by = 5, length.out = nrow(rawdata)/5),y]
        rawdata = rbind.data.frame(changefreq(1), cbind.data.frame("Time" = rowMeans(cbind(changefreq(3,1), changefreq(4,1)), na.rm = T), "CH2" = rowMeans(cbind(changefreq(3,2), changefreq(4,2)), na.rm = T)))
        rawdata = rawdata[order(rawdata$Time),]
      } else if (frequency == 2000) {
        rawdata[,1] = seq(0, by = 0.0005, length.out = nrow(rawdata)) # Time expressed in seconds and adjusted
      }
      
      ## Correct wave polarity
      if (exists("correctpol")) cpol = as.numeric(correctpol[correctpol$Date == Date & correctpol$Id == Id, -1:-2]) else cpol = (1)[-1]
      # if ((class(cpol) == "numeric" & length(cpol) > 0) | (class(cpol) == "data.frame" & nrow(cpol) > 0)) {
      cpol = cpol[!is.na(cpol)]
      if (length(cpol) > 0) {
        for(i in 1 + 1:length(cpol)) rawdata[,i] = rawdata[,i]*sign(cpol[i-1])
        rawdata = rawdata[,c(1,as.numeric(1 + abs(cpol)))]
      }
      ## Rename columns
      if (technique == "Telemetry") {
        names(rawdata) = c("Time", "CH2")
      } else names(rawdata) = c("Time", "CH1", "CH2", "CH3")
      
      
      ###----- 4. Segment files and normalize Waves
      ## Segment files
      segments = seq(0, max(rawdata[,1]), lengthsegment*60)
      if ((max(rawdata[,1]) - segments[length(segments)]*60)/(lengthsegment*60) < 0.5) segments = segments[-length(segments)]
      if (segment == F) segments = 0
      for (segm in segments) {
        tryCatch({
          if (segment) {
            segmdata = rawdata[rawdata[,1] >= segm & rawdata[,1] < segm + lengthsegment*60,]
            segmdata$Time = seq(0, by = 0.0005, length.out = nrow(segmdata)) # Reset time in each segment
          } else segmdata = rawdata
          
          ## Normalize Waves
          if (normalize) {
            if (technique == "Telemetry") channels = "CH2" else channels = c("CH1", "CH2", "CH3")
            for (channel in channels) {
              span = function (y,z) do.call("cbind", c(lapply(seq(y,1,-z), function (x) c(segmdata[-seq(1,x,z),channel],seq(1,x,z)/seq(1,x,z)*segmdata[nrow(segmdata),channel])),
                                                       list(segmdata[,channel]),
                                                       lapply(seq(1,y,z), function (x) c(seq(1,x,z)/seq(1,x,z)*segmdata[1,channel],segmdata[-seq(nrow(segmdata),(nrow(segmdata)-x+1),-z),channel]))))
              segmdata[,channel] = rowMeans2(span(noisespan, 1)) # Reduce noise
              segmdata[,paste0("Baseline", channel)] = rowMedians(span(basespan, 2)) # Estimate baseline
              segmdata[,channel] =  segmdata[,channel] - segmdata[,paste0("Baseline", channel)] # Baseline at 0
            }
          }

          ###----- 5. Export files and plots
          pdf(paste0(route, destination, Date, "/", gsub(" ", "", toupper(Id)), ifelse(segment, paste0(".", segm/60), ""), " ECG.pdf"), width = 8.27, height = 11.69, pointsize = 1)
          plotlines = seq(0,max(segmdata$Time), 2)
          par(mfrow = c(5, 1), bty = "l", cex.axis = 2, cex.lab = 2, pch = 19, cex = 0.5, lwd = 0.5, mar = c(4,5,3,3))
          invisible(sapply (plotlines, function (j){
            plot(NULL, type = "l", xlab = "Time (s)", yaxt = "n", ylab = "Lead", xaxs = "i", yaxs = "i", xlim = c(j, j + 2), ylim = c(-1,1), xaxp = c(j, j + 2, 10))
            axis(2, at = c(-0.5,0,0.5), labels = c("III", "II", "I"), tick = 0, las = 2)
            axis(2, at = c(0.5,0.75,1), labels = c("", "0.5 mV", ""), tick = 1, las = 3, line = 2, lwd = 1, col = 1)
            abline(h = c(-0.5, 0, 0.5), col = 1)
            for (channel in channels) lines((get(channel) + 1 - 0.5*which(c("CH1", "CH2", "CH3") == channel)) ~ Time, segmdata[segmdata$Time >= j & segmdata$Time <= j + 2,], col = ifelse(channel == "CH1", 2, ifelse(channel == "CH2", 4,3)))
          }))
          dev.off()

          fwrite(segmdata, file = paste0(route, destination, Date, "/", gsub(" ", "", toupper(Id)), ifelse(segment, paste0(".", segm/60), ""), ".txt"), row.names = F, col.names = F, sep = "\t", append = F, eol = "\r\n")
          setTxtProgressBar(txtProgressBar(style = 3), ((which(ecglist[,1] == Date & ecglist[,2] == Id) - 1)*length(segments) + which(segments == segm))/nrow(ecglist)*length(segments))
        }, error = function (e) cat("ERROR :", Date, " ", Id, " ", segm, " ", conditionMessage(e), "\n"))
      }
      if (!file.exists(paste0(route, destination, "processed.log.txt"))) fwrite(cbind.data.frame("Date", "Id", "Register", "Segmented", "Segment.length", "Normalized", "Processing.time"), file = paste0(route, destination, "processed.log.txt"), row.names = F, col.names = F, sep = "\t", append = F, eol = "\r\n")
      fwrite(cbind.data.frame(Date, gsub(" ", "", toupper(Id)), paste(Date, gsub(" ", "", toupper(Id))), segment, ifelse(segment, lengthsegment, NA), normalize, Sys.time()), file = paste0(route, destination, "processed.log.txt"), row.names = F, col.names = F, sep = "\t", append = T, eol = "\r\n")
      setTxtProgressBar(txtProgressBar(style = 3), which(ecglist[,1] == Date & ecglist[,2] == Id)/nrow(ecglist))
    }, error = function (e) cat("ERROR:", Date, " ", Id, " ", conditionMessage(e), "\n"))
  }
}

time1 = (proc.time() - time0)[[3]]
paste0(round(round(round(time1/60)/60)/24), "d ", round(round(time1/60)/60)%%24, "h ", round(time1/60)%%60, "min ", round(time1%%60), "s")
