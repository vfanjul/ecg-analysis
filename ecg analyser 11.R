time0 = proc.time()

###----- 1. Set parameters
project = "Heart progerin" # mir29   Heart progerin   Heart lamin  ISO Challenge  PCTX Treatment  Zmpste-Rankl  HGPS Amanda DBU Alberto
technique = "ECG" # Telemetry ECG Thermochallenge.ECG ECG.extra
channel = "CH2" # CH1 CH2 CH3
processall = F # Default is F

peakspan = 5
thresholdspan = 4000
qthreshold = 0.02
q0threshold = -0.01 # -0.02 is too much (41) and 0 too little (13), -0.01 is 9
p0threshold = 0.02
qslope = 0.01
pslope = 0.001
tslope = 0.01
rrmin = 30 # 30 ms/b is 4000 bpm
rrmax = 800 # 800 ms/b is 75 bpm
P0Pemin = 7
PQ0min = 20
Q0Qmax = 4
Q0Temin = 30
STnmin = 10
QRmin = 2 # (with 3, 4 excluded incorrectly, 13 excluded and 1 not excluded)

waves = c("P0", "P", "Pe", "Q0", "Q", "R", "S", "J", "T", "T90", "Tn", "Te")
wavesplot = c(expression("P"[0]), "P", expression("P"["e"]), expression("Q"[0]), "Q", "R", "S", "J", "T", expression("T"[90]), expression("T"["n"]), expression("T"["e"]))

packages = c("data.table", "matrixStats", "tidyr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) install.packages(setdiff(packages, rownames(installed.packages())))  
for (i in packages) library(i, character.only = T)

if (.Platform$OS.type == "unix") setwd("/Volumes/Victor/") else setwd("S:/LAB_VA/LAB/Victor/")
# setwd("/Users/Victor/Downloads/")
# setwd("S:/LAB_VA/LAB/Alvaro Macias/Results/Esperanza de vida/Challenge de Isoproterenol/Electros/")
# setwd("S:/LAB_VA/LAB/Alvaro Macias/Results/Esperanza de vida/Paclitaxel 0.5mg_kg 3_semana/ECGs/")
# setwd("/Volumes/LAB/Alberto M/Alberto/Ratones endotelio y SMC/ECG DBU/")
if(project == "PCTX Treatment") setwd("/Volumes/LAB_VA/LAB/Alvaro Macias/Results/Tratamientos/Paclitaxel 0.5mg_kg 3_semana/ECGs sucesivos tras PCTX agudo - REAL/")

baseroute = paste0(project, " project/Raw data/", technique, "/")
dataroute = paste0(baseroute,"normalized/")
if (technique == "Telemetry") channel = "CH2"
route = paste0(baseroute, "analyzed ", channel, "/")
dir.create(route, showWarnings = F)

registers = sapply(dir(dataroute), function (x) list(gsub(".txt", "", dir(paste0(dataroute, x))[grep(".txt", dir(paste0(dataroute, x)))])))
ecglist = c()
for (i in names(registers)) if(length(registers[[i]]) != 0) ecglist = rbind.data.frame(ecglist, cbind(i,registers[[i]], paste(i,registers[[i]])), stringsAsFactors = F)
if (!processall) {
  ecglista = c()
  if (file.exists(paste0(route, "processed.log.txt"))) ecglista = data.frame(fread(paste0(route, "processed.log.txt"), encoding = "UTF-8", header = T, stringsAsFactors = F, dec = ".", na.strings = c("N/A", "NA", "nan"), sep = "\t"))[,1:3]
  ecglist = merge(ecglist, setdiff(as.vector(ecglist[,3]), as.vector(ecglista[,3])), by.x = 3, by.y = 1)
  ecglist = ecglist[,c(2,3,1)]
}
# ecglist = ecglist[grepl("2018",ecglist$i),]


Date = unique(ecglist[,1])[1]
Id = ecglist[ecglist[,1] == Date, 2][1]
for (Date in unique(ecglist[,1])) {
  dir.create(paste0(route, Date), showWarnings = F)
  for (Id in ecglist[ecglist[,1] == Date,2]) {
    tryCatch({
      ###----- 2. Import data
      rawdata = fread(paste0(dataroute, Date, "/", Id, ".txt"), encoding = "UTF-8", header = F, stringsAsFactors = F, dec = ".", na.strings = c("N/A", "NA", "nan"))
      if (technique == "Telemetry") {
        names(rawdata) = c("Time", "CH2", "BaselineCH2")
      } else names(rawdata) = c("Time", "CH1", "CH2", "CH3", "BaselineCH1", "BaselineCH2", "BaselineCH3")
      rawdata = data.frame(rawdata)
      
      ## Peaks and point slope
      rawdata[is.na(rawdata[,channel]),channel] = 0
      rawdata[is.na(rawdata[,paste0("Baseline", channel)]),paste0("Baseline", channel)] = 10
      rawdata$Peak = 0
      span = function (y,z) do.call("cbind", c(lapply(seq(y,1,-z), function (x) c(rawdata[-seq(1,x,z),channel],seq(1,x,z)/seq(1,x,z)*rawdata[nrow(rawdata),channel])),
                                               list(rawdata[,channel]),
                                               lapply(seq(1,y,z), function (x) c(seq(1,x,z)/seq(1,x,z)*rawdata[1,channel],rawdata[-seq(nrow(rawdata),(nrow(rawdata)-x+1),-z),channel]))))
      rawdata[rawdata[,channel] == rowMaxs(span(peakspan, 1)),"Peak"] = 1 # Positive peaks
      rawdata[rawdata[,channel] == rowMins(span(peakspan, 1)),"Peak"] = -1 # Negative peaks
      rawdata$Slope = (c(rawdata[-1,channel],0) - c(0,rawdata[-nrow(rawdata),channel]))/2 # Point slope
      
      ###----- 3. Detect peaks & waves
      ## R and S thresholds
      interval = seq(1,(nrow(rawdata) - thresholdspan),thresholdspan)
      outsr = rawdata[(nrow(rawdata) - thresholdspan):nrow(rawdata),][rawdata[(nrow(rawdata) - thresholdspan):nrow(rawdata),"Peak"] == 1 & rawdata[(nrow(rawdata) - thresholdspan):nrow(rawdata),channel] > 0,channel]
      outss = rawdata[(nrow(rawdata) - thresholdspan):nrow(rawdata),][rawdata[(nrow(rawdata) - thresholdspan):nrow(rawdata),"Peak"] == -1 & rawdata[(nrow(rawdata) - thresholdspan):nrow(rawdata),channel] < 0,channel]
      if (length(outsr) > 5 & min(outsr) != max(outsr) & length(outss) > 5 & min(outss) != max(outss)) {
        clustr = cbind(outsr,kmeans(outsr,2)$cluster)
        rawdata$Rthreshold = min(max(clustr[clustr[,2] == 1,1]), max(clustr[clustr[,2] == 2,1]))
        clusts = cbind(outss,kmeans(outss,2)$cluster)
        rawdata$Sthreshold = max(min(clusts[clusts[,2] == 1,1]), min(clusts[clusts[,2] == 2,1]))
      } else {
        rawdata$Rthreshold = 1
        rawdata$Sthreshold = -1
      }
      invisible(sapply(interval, function (i) {
        outsr = rawdata[i:(i + interval[2]),][rawdata[i:(i + interval[2]),"Peak"] == 1 & rawdata[i:(i + interval[2]),channel] > 0,channel]
        outss = rawdata[i:(i + interval[2]),][rawdata[i:(i + interval[2]),"Peak"] == -1 & rawdata[i:(i + interval[2]),channel] < 0,channel]
        if (length(outsr) > 5 & min(outsr) != max(outsr) & length(outss) > 5 & min(outss) != max(outss)) {
          clustr = cbind(outsr,kmeans(outsr,2)$cluster)
          # rawdata$Rthreshold[i:(i + interval[2])] = min(max(clustr[clustr[,2] == 1,1]), max(clustr[clustr[,2] == 2,1]))
          rawdata$Rthreshold[i:(i + interval[2])] <<- max(min(clustr[clustr[,2] == 1,1]), min(clustr[clustr[,2] == 2,1]))
          clusts = cbind(outss,kmeans(outss,2)$cluster)
          # rawdata$Sthreshold[i:(i + interval[2])] = max(min(clusts[clusts[,2] == 1,1]), min(clusts[clusts[,2] == 2,1]))
          rawdata$Sthreshold[i:(i + interval[2])] <<- min(max(clusts[clusts[,2] == 1,1]), max(clusts[clusts[,2] == 2,1]))
        } else {
          rawdata$Rthreshold[i:(i + interval[2])] <<- 1
          rawdata$Sthreshold[i:(i + interval[2])] <<- -1
        }
      }))
      rm(list = c(ls(pattern = "outs"), ls(pattern = "clust")))

      
      ## QRS Complex
      rawdata$Wave = ""
      # rawdata[(rawdata$Peak == 1 & rawdata[,channel] > rawdata$Rthreshold & rawdata$Rthreshold > abs(rawdata$Sthreshold)),"Wave"] = "R"
      rawdata[(rawdata$Peak == 1 & rawdata[,channel] >= rawdata$Rthreshold & rawdata$Rthreshold > abs(rawdata$Sthreshold)),"Wave"] = "R"
      rawdata[(rawdata$Peak == 1 | rawdata$Wave == "R"),"Wave"][c("", rawdata[(rawdata$Peak == 1 | rawdata$Wave == "R"),"Wave"][-nrow(rawdata[(rawdata$Peak == 1 | rawdata$Wave == "R"),])]) == "R"] = ""
      rawdata[(rawdata$Peak == -1 | rawdata$Wave == "R"),"Wave"][c("", rawdata[(rawdata$Peak == -1 | rawdata$Wave == "R"),"Wave"][-nrow(rawdata[(rawdata$Peak == -1 | rawdata$Wave == "R"),])]) == "R"] = "S"
      # rawdata[(rawdata$Peak == -1 & rawdata[,channel] < rawdata$Sthreshold & rawdata$Rthreshold < abs(rawdata$Sthreshold)),"Wave"] = "S"
      rawdata[(rawdata$Peak == -1 & rawdata[,channel] <= rawdata$Sthreshold & rawdata$Rthreshold < abs(rawdata$Sthreshold)),"Wave"] = "S"
      rawdata[(rawdata$Peak == -1 | rawdata$Wave == "S"),"Wave"][c(rawdata[(rawdata$Peak == -1 | rawdata$Wave == "S"),"Wave"][-1], "") == "S"] = ""
      rawdata[(rawdata$Peak == 1 | rawdata$Wave == "S"),"Wave"][c(rawdata[(rawdata$Peak == 1 | rawdata$Wave == "S"),"Wave"][-1], "") == "S"] = "R"
      rawdata[rawdata[,paste0("Baseline", channel)] > 9.9 | rawdata[,paste0("Baseline", channel)] < -9.9,"Wave"] = "" # Remove artifacts when signal is lost
      rawdata = subset(rawdata, select = -c(Rthreshold, Sthreshold))
      

      rawdata[(rawdata$Peak == 1 | rawdata$Wave == "S"),"Wave"][c("", rawdata[(rawdata$Peak == 1 | rawdata$Wave == "S"),"Wave"][-nrow(rawdata[(rawdata$Peak == 1 | rawdata$Wave == "S"),])]) == "S"] = "T"
      # rawdata[(rawdata$Slope < qslope | rawdata$Wave == "R"),"Wave"][c(rawdata[(rawdata$Slope < qslope | rawdata$Wave == "R"),"Wave"][-1], "") == "R"] = "Q"
      rawdata[(rawdata$Slope < qslope & rawdata[,channel] < qthreshold | rawdata$Wave == "R" ),"Wave"][c(rawdata[(rawdata$Slope < qslope & rawdata[,channel] < qthreshold | rawdata$Wave == "R"),"Wave"][-1], "") == "R"] = "Q"
      # rawdata[((abs(rawdata[,channel]) < qthreshold & abs(rawdata$Slope) < qslope) | rawdata$Wave == "Q"),"Wave"][c(rawdata[((abs(rawdata[,channel]) < qthreshold & abs(rawdata$Slope) < qslope) | rawdata$Wave == "Q"),"Wave"][-1:-2], "", "") == "Q"] = "Q0"
      # rawdata[(abs(rawdata[,channel]) < qthreshold & (abs(rawdata$Slope) < qslope/2 | rawdata$Peak == 1) & rawdata$Wave == "") | rawdata$Wave == "Q","Wave"][c(rawdata[(abs(rawdata[,channel]) < qthreshold & (abs(rawdata$Slope) < qslope/2 | rawdata$Peak == 1) & rawdata$Wave == "") | rawdata$Wave == "Q","Wave"][-1:-2], "", "") == "Q"] = "Q0"
      # rawdata[(rawdata[,channel] > 0 & (abs(rawdata$Slope) < qslope | rawdata$Peak == 1) & rawdata$Wave == "") | rawdata$Wave == "Q","Wave"][c(rawdata[(rawdata[,channel] > 0 & (abs(rawdata$Slope) < qslope | rawdata$Peak == 1) & rawdata$Wave == "") | rawdata$Wave == "Q","Wave"][-1:-2], "", "") == "Q"] = "Q0"
      rawdata[(rawdata[,channel] > q0threshold & (abs(rawdata$Slope) < qslope | rawdata$Peak == 1) & rawdata$Wave == "") | rawdata$Wave == "Q","Wave"][c(rawdata[(rawdata[,channel] > q0threshold & (abs(rawdata$Slope) < qslope | rawdata$Peak == 1) & rawdata$Wave == "") | rawdata$Wave == "Q","Wave"][-1:-2], "","") == "Q"] = "Q0"
      

      ## Number of beats
      rawdata[rawdata$Wave == "R","BeatR"] = 1:(nrow(rawdata[rawdata$Wave == "R",]))
      rawdata[1,"BeatR"] = 0
      rawdata = fill(rawdata, BeatR)
      rawdata$BeatR[abs(rawdata[,paste0("Baseline", channel)]) > 9.9] = rawdata$BeatR[abs(rawdata[,paste0("Baseline", channel)]) > 9.9] + 0.1 # Regions where signal is offline
      
      ## P wave
      wavedata = rawdata[rawdata$Peak == 1 & rawdata$Wave == "",]
      pdata = data.frame("Pwave" = 1, sapply(0:(max(rawdata[,"BeatR"])-1), function (x) max(wavedata[wavedata$BeatR == x, channel])), "BeatR" = 0:(max(rawdata$BeatR)-1))
      rawdata = merge(rawdata, pdata, by.x = c(channel,"BeatR"), by.y = c(2,3), all.x = T, sort = F)
      rawdata = rawdata[order(rawdata$Time),]
      rawdata[!is.na(rawdata$Pwave) & rawdata$Wave == "","Wave"] = "P"
      rawdata[duplicated(rawdata[,c("BeatR", "Wave")]),"Wave"] = ""
      # rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][c(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][-1:-2], "", "") == "P"] = "P0"
      # rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][c("", "", rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][-nrow(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata$Wave == "") | rawdata$Wave == "P",]):-(nrow(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata$Wave == "") | rawdata$Wave == "P",])-1)]) == "P"] = "Pe"
      # rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & abs(rawdata[,channel]) < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][c(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & abs(rawdata[,channel]) < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][-1:-2], "", "") == "P"] = "P0"
      rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][c(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][-1:-2], "", "") == "P"] = "P0"
      # rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & abs(rawdata[,channel]) < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][c("", "", rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & abs(rawdata[,channel]) < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][-nrow(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & abs(rawdata[,channel]) < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P",]):-(nrow(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & abs(rawdata[,channel]) < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P",])-1)]) == "P"] = "Pe"
      # rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][c("", "", rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][-nrow(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P",]):-(nrow(rawdata[((abs(rawdata$Slope) < pslope | rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P",])-1)]) == "P"] = "Pe"
      rawdata[((rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][c("", rawdata[((rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P","Wave"][-nrow(rawdata[((rawdata$Peak == -1) & rawdata[,channel] < p0threshold & rawdata$Wave == "") | rawdata$Wave == "P",])]) == "P"] = "Pe"
      rawdata[duplicated(rawdata[,c("BeatR", "Wave")]),"Wave"] = ""
      rawdata = subset(rawdata, select = -Pwave)
      
      ## Correct begining of beat
      rawdata[rawdata$Wave == "P0" | rawdata$Wave == "P" | rawdata$Wave == "Pe" | rawdata$Wave == "Q0" | rawdata$Wave == "Q","BeatP"] = rawdata[rawdata$Wave == "P0" | rawdata$Wave == "P" | rawdata$Wave == "Pe" | rawdata$Wave == "Q0" | rawdata$Wave == "Q","BeatR"] + 1
      rawdata[rawdata$Wave == "R","BeatP"] = rawdata[rawdata$Wave == "R","BeatR"]
      rawdata[1,"BeatP"] = 0
      rawdata = fill(rawdata, BeatP)
      rawdata[rawdata$Wave == "T","BeatT"] = rawdata[rawdata$Wave == "T","BeatR"] + 1
      rawdata[1,"BeatT"] = 0
      rawdata = fill(rawdata, BeatT)
      
      ## T wave
      wavedata = rawdata[rawdata$Peak != 0 & rawdata$Wave == "",]
      tndata = data.frame("Tnwave" = 1, sapply(0:(max(rawdata[,"BeatR"])-1), function (x) min(wavedata[wavedata$BeatR == x & wavedata$BeatP == x & wavedata$BeatT == x+1, channel])), "BeatR" = 0:(max(rawdata$BeatR)-1))
      rawdata = merge(rawdata, tndata, by.x = c(channel,"BeatR"), by.y = c(2,3), all.x = T, sort = F)
      rawdata = rawdata[order(rawdata$Time),]
      rawdata[!is.na(rawdata$Tnwave) & rawdata$Wave == "","Wave"] = "Tn"
      rawdata[duplicated(rawdata[,c("BeatR", "Wave")]),"Wave"] = ""
      rawdata[(rawdata[,channel] >= 0 | rawdata$Wave == "Tn"),"Wave"][c("", rawdata[(rawdata[,channel] >= 0 | rawdata$Wave == "Tn"),"Wave"][-nrow(rawdata[(rawdata[,channel] >= 0 | rawdata$Wave == "Tn"),])]) == "Tn"] = "Te"
      rawdata = subset(rawdata, select = -Tnwave)
      
      ## Correct T wave
      rawdata[rawdata$Wave == "Tn","BeatTn"] = rawdata[rawdata$Wave == "Tn","BeatR"] + 1
      rawdata[1,"BeatTn"] = 0
      rawdata = fill(rawdata, BeatTn)
      rawdata[rawdata$Wave == "T", "Wave"] = ""
      
      wavedata = rawdata[rawdata$Peak == 1 & rawdata$Wave == "",]
      tdata = data.frame("Twave" = 1, sapply(0:(max(rawdata[,"BeatR"])-1), function (x) max(wavedata[wavedata$BeatR == x & wavedata$BeatP == x & wavedata$BeatT == x+1 & wavedata$BeatTn == x, channel])), "BeatR" = 0:(max(rawdata$BeatR)-1))
      rawdata = merge(rawdata, tdata, by.x = c(channel,"BeatR"), by.y = c(2,3), all.x = T, sort = F)
      rawdata = rawdata[order(rawdata$Time),]
      rawdata[!is.na(rawdata$Twave) & rawdata$Wave == "","Wave"] = "T"
      rawdata[duplicated(rawdata[,c("BeatR", "Wave")]),"Wave"] = ""
      rawdata = subset(rawdata, select = -Twave)
      
      # rawdata[(abs(rawdata$Slope) <= tslope | rawdata$Wave == "S"),"Wave"][c("", rawdata[(abs(rawdata$Slope) <= tslope | rawdata$Wave == "S"),"Wave"][-nrow(rawdata[(abs(rawdata$Slope) <= tslope | rawdata$Wave == "S"),])]) == "S"] = "J"
      # rawdata[(abs(rawdata$Slope) > tslope | rawdata$Wave == "T"),"Wave"][c(rawdata[(abs(rawdata$Slope) > tslope | rawdata$Wave == "T"),"Wave"][-1], "") == "T"] = "J"
      # rawdata[(rawdata$Slope > tslope | rawdata$Wave == "T"),"Wave"][c(rawdata[(rawdata$Slope > tslope | rawdata$Wave == "T"),"Wave"][-1], "") == "T"] = "J"
      # rawdata[(rawdata$Slope <= tslope | rawdata$Wave == "S"),"Wave"][c("", rawdata[(rawdata$Slope <= tslope | rawdata$Wave == "S"),"Wave"][-nrow(rawdata[(rawdata$Slope <= tslope | rawdata$Wave == "S"),])]) == "S"] = "J"
      # rawdata[(rawdata$Slope <= tslope & rawdata$Peak == 0) | rawdata$Wave == "S","Wave"][c("", rawdata[(rawdata$Slope <= tslope & rawdata$Peak == 0) | rawdata$Wave == "S","Wave"][-nrow(rawdata[(rawdata$Slope <= tslope & rawdata$Peak == 0) | rawdata$Wave == "S",])]) == "S"] = "J"
      rawdata[(rawdata$Slope <= tslope & rawdata$Peak == 0) | rawdata$Wave == "S","Wave"][c("", "", rawdata[(rawdata$Slope <= tslope & rawdata$Peak == 0) | rawdata$Wave == "S","Wave"][-nrow(rawdata[(rawdata$Slope <= tslope & rawdata$Peak == 0) | rawdata$Wave == "S",]):-(nrow(rawdata[(rawdata$Slope <= tslope & rawdata$Peak == 0) | rawdata$Wave == "S",])-1)]) == "S"] = "J"
      
      t90data = data.frame("BeatR" = unique(rawdata$BeatR))
      for (i in c("T","Tn")) t90data = merge(t90data, rawdata[rawdata$Wave == i,c("BeatR", channel)], by = "BeatR", all.x = T, sort = F)
      names(t90data) = c("BeatR", "T","Tn")
      t90data$T90threshold = 0.1*t90data$T + 0.9*t90data$Tn
      t90data = t90data[!duplicated(t90data[,c("BeatR")]),]
      rawdata = merge(rawdata, t90data[,c(1,4)], by = "BeatR", all.x = T, sort = F)
      rawdata = rawdata[order(rawdata$Time),]
      rawdata$T90threshold[is.na(rawdata$T90threshold)] = 0
      rawdata[((rawdata[,channel] - rawdata$T90threshold) <= 0 | rawdata$Wave == "T"),"Wave"][c("", rawdata[((rawdata[,channel] - rawdata$T90threshold) <= 0 | rawdata$Wave == "T"),"Wave"][-nrow(rawdata[((rawdata[,channel] - rawdata$T90threshold) <= 0 | rawdata$Wave == "T"),])]) == "T"] = "T90"
      rawdata = subset(rawdata, select = -T90threshold)
      
      
      ## Eliminate duplicated waves and empty rows
      rawdata[duplicated(rawdata[,c("BeatR", "Wave")]),"Wave"] = ""
      rawdata = rawdata[complete.cases(rawdata$Time),]
      
      ## Time relative to R in each BeatR
      RTime = rawdata[rawdata$Wave == "R", c("BeatP", "Time")]
      names(RTime) = c("BeatP", "RTime")
      rawdata = merge(rawdata, RTime, by = "BeatP", all.x = T, sort = F)
      rawdata = rawdata[order(rawdata$Time),]
      rawdata$TimeR = rawdata$Time - rawdata$RTime
      rawdata = subset(rawdata, select = -RTime)
      rm(list = c("wavedata", "pdata", "tndata", "tdata", "t90data", "RTime"))
      
      
      # nwaves = c(); sapply(waves, function (x) nwaves = c(nwaves,nrow(rawdata[rawdata$Wave == x,])))
      
      ###----- 4. Exclude outliers and calculate intervals
      ## Time of wave points in each beat
      ecgtime = data.frame("BeatR" = unique(rawdata$BeatR))
      for (i in waves) ecgtime = merge(ecgtime, rawdata[rawdata$Wave == i,c("BeatR", "Time")], by = "BeatR", all.x = T, sort = F)
      names(ecgtime) = c("BeatR", waves)
      for (i in waves[1:5]) ecgtime[,i] = c(NA,ecgtime[-nrow(ecgtime),i])

      ## Voltage of wave points in each beat and outliers
      ecgvoltage = data.frame("BeatR" = unique(rawdata$BeatR))
      for (i in waves) ecgvoltage = merge(ecgvoltage, rawdata[rawdata$Wave == i,c("BeatR", channel)], by = "BeatR", all.x = T, sort = F)
      names(ecgvoltage) = c("BeatR", waves)
      for (i in waves[1:5]) ecgvoltage[,i] = c(NA,ecgvoltage[-nrow(ecgvoltage),i])
      for (i in waves) {
        ecgvoltage[,paste0(i,"out")] = 1
        ecgvoltage[!is.na(ecgvoltage[,i]) & ecgvoltage[,i] >= boxplot(ecgvoltage[,i], plot = F)$stats[1] & ecgvoltage[,i] <= boxplot(ecgvoltage[,i], plot = F)$stats[5],paste0(i,"out")] = 0
      }
      ecgvoltage$Rout = 0
      
      ## Time of wave points in each beat (related to R peak) and outliers
      ecgRTime = ecgtime - ecgtime$R
      ecgRTime$BeatR = ecgtime$BeatR
      for (i in waves) {
        ecgRTime[,paste0(i,"out")] = 1
        ecgRTime[!is.na(ecgRTime[,i]) & ecgRTime[,i] >= boxplot(ecgRTime[,i], plot = F)$stats[1] & ecgRTime[,i] <= boxplot(ecgRTime[,i], plot = F)$stats[5],paste0(i,"out")] = 0
      }
      ecgRTime$Rout = 0
      
      ## Outliers exclusion
      ecgdata = merge(cbind.data.frame(ecgtime, ecgRTime[,as.vector(sapply(waves, function (x) paste0(x, "out")))]), ecgvoltage, by = "BeatR", all = T)
      for (i in waves) {
        ecgdata[!is.na(ecgdata[,paste0(i,"out.x")]) & !is.na(ecgdata[,paste0(i,"out.y")]) & (ecgdata[,paste0(i,"out.x")] == 1 | ecgdata[,paste0(i,"out.y")] == 1) ,paste0(i,".x")] = NA
        ecgdata[!is.na(ecgdata[,paste0(i,"out.x")]) & !is.na(ecgdata[,paste0(i,"out.y")]) & (ecgdata[,paste0(i,"out.x")] == 1 | ecgdata[,paste0(i,"out.y")] == 1) ,paste0(i,".y")] = NA
      }
      ecgdata = ecgdata[order(ecgdata$BeatR),]
      ecgdata = ecgdata[-1,]
      rm(list = c("ecgtime", "ecgvoltage", "ecgRTime"))

      ## Interval calculation
      ecgint = data.frame("Beat" = ecgdata$BeatR)
      ecgint$Time = ecgdata$R.x
      ecgint$RR = (c(ecgdata$R.x[-1],0) - ecgdata$R.x)*1000
      ecgint[(!is.na(ecgint$RR) & (ecgint$RR < rrmin | ecgint$RR > rrmax)),"RR"] = NA
      ecgint$HR = 60*1000/ecgint$RR
      ecgint$PQ = (ecgdata$Q0.x - ecgdata$P.x)*1000
      
      ecgint$QRS = (ecgdata$S.x - ecgdata$Q0.x)*1000
      ecgint$QJ = (ecgdata$J.x - ecgdata$Q0.x)*1000
      ecgint$QT = (ecgdata$T.x - ecgdata$Q0.x)*1000
      ecgint$QT90 = (ecgdata$T90.x - ecgdata$Q0.x)*1000
      ecgint$QTe = (ecgdata$Te.x - ecgdata$Q0.x)*1000
      ecgint$QT90c = ecgint$QT90/sqrt(ecgint$RR/100)
      ecgint$QTec = ecgint$QTe/sqrt(ecgint$RR/100)
      
      ecgint$P0Pe = (ecgdata$Pe.x - ecgdata$P0.x)*1000
      ecgint$JT = (ecgdata$T.x - ecgdata$J.x)*1000
      ecgint$JT = ifelse(ecgint$JT<0, 0, ecgint$JT)
      ecgint$TT90 = (ecgdata$T90.x - ecgdata$T.x)*1000
      # ecgint$TTe = (ecgdata$Te.x - ecgdata$T.x)*1000
      namesecgint = sapply(names(ecgint)[-1:-4], function(x) paste(x, "(ms)"))
      ecgint$Tsteep = abs(ecgdata$T.y - ecgdata$T90.y)/(ecgdata$T90.x - ecgdata$T.x)
      ecgint$Psteep = abs(ecgdata$P.y - (ecgdata$Pe.y + ecgdata$P0.y)/2)/(ecgdata$Pe.x - ecgdata$P0.x)
      namesecgint = c(namesecgint, "T steepness (mV/s)", "P steepness (mV/s)")
      ecgint = ecgint[-nrow(ecgint),]
      ecgint = ecgint[ecgint$Beat != 0,]

      ## RR variation
      RRvar = ecgint$RR[-1] - ecgint$RR[-nrow(ecgint)]
      HRV = data.frame(
        "RMSSD" = sqrt(mean(RRvar^2, na.rm = T)),
        "SDRR" = sd(ecgint$RR, na.rm = T),
        "SDSD" = sd(RRvar, na.rm = T),
        "pNN6 (%)" = length(RRvar[abs(RRvar) > 6 & !is.na(RRvar)])/(length(RRvar[!is.na(RRvar)])-1)*100
      )
      
      ## Median wave
      mediandata = rawdata[!is.na(rawdata$TimeR) & rawdata$BeatR == ceiling(rawdata$BeatR), c(channel,"TimeR","Wave", "BeatR", "BeatP")]
      mediandata$TimeR = round(round(mediandata$TimeR*1000, 1)*2, 0)/2
      mediandata = mediandata[mediandata$TimeR >=  min(sapply(waves, function (x) median((ecgdata[,paste0(x,".x")] - ecgdata$R.x)*1000, na.rm = T))) - 5 & mediandata$TimeR <=  max(sapply(waves, function (x) median((ecgdata[,paste0(x,".x")] - ecgdata$R.x)*1000, na.rm = T))) + 5,]
      mediandata = mediandata[complete.cases(mediandata),]
      mediandata = mediandata[order(mediandata$TimeR),]
      
      medianvoltage = data.frame(unique(mediandata$TimeR), do.call("rbind", lapply(unique(mediandata$TimeR), function(x) boxplot.stats(mediandata[mediandata$TimeR == x,channel])$stats)))
      names(medianvoltage) = c("TimeR", "Q0V", "Q1V", "MV", "Q3V", "Q4V")
      for (i in names(medianvoltage)[-1]) medianvoltage[,i] = apply(cbind(c(medianvoltage[-1,i],0), medianvoltage[,i], c(0,medianvoltage[-nrow(medianvoltage),i])), 1, mean)
      # mediandata = merge(mediandata, medianvoltage, by = "TimeR", all.x = T)
      # mediandata = mediandata[!duplicated(mediandata[,c("TimeR", "BeatR")]),]
      medianwaves = data.frame("Wave" = "", "TimeR" = min(medianvoltage$TimeR), stringsAsFactors = F)
      for (i in waves) medianwaves = rbind.data.frame(medianwaves, cbind.data.frame("Wave" = i, "TimeR" = round(round(median((ecgdata[,paste0(i, ".x")] - ecgdata[,"R.x"])*1000, na.rm = T), 1)*2, 0)/2), stringsAsFactors = F)
      medianvoltage = merge(medianvoltage, medianwaves, by = "TimeR", all = T)
      medianvoltage[is.na(medianvoltage$Wave),"Wave"] = ""
      
      ## Median intervals
      medecgint = rbind.data.frame(c(max(rawdata$Time), median(sapply(2:(length(waves)+1), function (x) length(ecgdata[!is.na(ecgdata[,x]), x])), na.rm = T), max(ecgint$Beat), 0, sapply(3:ncol(ecgint), function (x) median(ecgint[,x], na.rm = T))))
      names(medecgint) = c("Duration (min)", "MedianBeats", "TotalBeats", "RelBeats (%)", "RR (ms)", "HR (bpm)", namesecgint)
      medecgint$`RelBeats (%)` = round(medecgint$MedianBeats/medecgint$TotalBeats*100,2)
      medecgint = cbind.data.frame(medecgint, HRV)
      if ((medianvoltage[medianvoltage$Wave == "Pe","TimeR"] - medianvoltage[medianvoltage$Wave == "P0","TimeR"]) < P0Pemin | medecgint$`PQ (ms)` < PQ0min | medecgint$`P0Pe (ms)` > 30) {
        medecgint[grepl("^P", names(medecgint))] = NA
        medecgint$`Absent P` = 1
      } else medecgint$`Absent P` = 0
      
      # for (i in c("P", "Q", "R", "S", "T")) medecgint[paste0(i, " (mV)")] = median(ecgdata[,paste0(i,".y")], na.rm = T) + 1
      if (medecgint$`QTe (ms)` < Q0Temin | median((ecgdata$Tn.x - ecgdata$S.x)*1000, na.rm = T) < STnmin) medecgint[grepl("T", names(medecgint))][-1:-2] = NA
      if (median((ecgdata$Q.x - ecgdata$Q0.x)*1000, na.rm = T) > Q0Qmax & medecgint$`QRS (ms)` > 15) medecgint[grepl("Q", names(medecgint))] = NA
      
      ## ECG waves
      ecgwaves = merge(rawdata[rawdata$Wave != "", c("Time", channel, "BeatP", "BeatR", "Wave")], data.frame("Time" = sort(unique(unlist(ecgdata[,2:(length(waves)+1)]))), "Outlier" = 0), all.x = T)
      ecgwaves$Outlier[is.na(ecgwaves$Outlier)] = 1
      ecgwaves = cbind.data.frame("Id" = Id, "Date" = Date, ecgwaves)

      
      
      ###----- 5. Export plots and results
      Id = gsub(" ", "", toupper(Id))
      ## Wave point identification plots
      pdf(paste0(route, Date, "/", Id, " ", Date, " wave points variation.pdf"), 12, 9, pointsize = 12)
      par(mfrow = c(3,4), bty = "l", cex.axis = 1.5, cex.lab = 1.5, pch = 19, cex = 0.5, mar = c(5,5,2,2))
      for (i in waves) {
        plot((ecgdata[,paste0(i, ".x")] - ecgdata[,"R.x"])*1000, ecgdata[,paste0(i, ".y")], ylab = paste0("Voltage ", i, " (mV)"), xlab = paste0("Time ", i, " - R (ms)"), col = rainbow(max(ecgdata$BeatR), start = 0, end = 2/6))
        abline(h = boxplot.stats(ecgdata[,paste0(i, ".y")])$stats, col = 8, lty = c(3,2,1,2,3))
        abline(v = boxplot.stats((ecgdata[,paste0(i, ".x")] - ecgdata[,"R.x"])*1000)$stats, col = 8, lty = c(3,2,1,2,3))
      }
      dev.off()
      
      pdf(paste0(route, Date, "/", Id, " ", Date, " ECG.pdf"), width = 8.27, height = 11.69, pointsize = 1)
      plotlines = seq(0,max(rawdata$Time), 2)
      par(mfrow = c(5, 1), bty = "l", cex.axis = 2, cex.lab = 2, pch = 19, cex = 0.5, lwd = 0.5, mar = c(4,5,3,3))
      invisible(sapply (plotlines, function (j){
        plot(NULL, type = "l", xlab = "Time (s)", ylab = "Voltage (mV)", xaxs = "i", yaxs = "i", xlim = c(j, j + 2), ylim = c(-1,1), xaxp = c(j, j + 2, 10))
        abline(h = 0, col = 1)
        # for (i in waves) {
          # abline(v = rawdata[rawdata$Time >= j & rawdata$Time <= j + 2 & rawdata$Wave == i,"Time"], col = rainbow(length(waves))[which(waves == i)])
          # abline(v = rawdata[rawdata$Wave == i,"Time"], col = rainbow(length(waves))[which(waves == i)])
          # points(ecgdata[,paste0(i,".y")] ~ ecgdata[,paste0(i,".x")], col = 2, na.omit = T)
          # points(ecgdata[,paste0(i,".y")] ~ ecgdata[,paste0(i,".x")], col = rainbow(length(waves))[which(waves == i)], na.omit = T)
          # points(get(paste0(i,".y")) ~ get(paste0(i,".x")), ecgdata[ecgdata[,paste0(i,".x")] >= j & ecgdata[,paste0(i,".x")] <= j + 2,], col = rainbow(length(waves))[which(waves == i)], na.omit = T)
        # }
          sapply(waves, function(i) {
            abline(v = rawdata[rawdata$Wave == i,"Time"], col = rainbow(length(waves))[which(waves == i)])
            points(ecgdata[,paste0(i,".y")] ~ ecgdata[,paste0(i,".x")], col = rainbow(length(waves))[which(waves == i)], na.omit = T, cex = 2)
            })
          
        lines(get(channel) ~ Time, rawdata[rawdata$Time >= j & rawdata$Time <= j + 2,])
      }))
      dev.off()
      
      
      ## Interval calculation, RR variation and summary plots
      plots_per_col = ifelse(length(names(ecgint)[-1:-4]) == 3, 1, round(sqrt(length(names(ecgint)[-1:-4])),0))
      plots_per_row = ifelse(length(names(ecgint)[-1:-4]) == 3, 3, ceiling(length(names(ecgint)[-1:-4])/plots_per_col))
      pdf(paste0(route, Date, "/", Id, " ", Date, " parameter variation.pdf"), plots_per_row*2, plots_per_col*2, pointsize = 12)
      par(mfrow = c(plots_per_col,plots_per_row), bty = "l", cex.axis = 1.5, cex.lab = 1.5, pch = 19, cex = 0.5, mar = c(5,5,2,2))
      for (i in 5:ncol(ecgint)) {
        plot(ecgint[,i], ylab = namesecgint[i-4], xlab = "Beat", col = rainbow(max(ecgdata$BeatR), start = 0, end = 2/6))
        abline(h = boxplot.stats(ecgint[,i])$stats, col = 8, lty = c(3,2,1,2,3))
      }
      dev.off()
      
      pdf(paste0(route, Date, "/", Id, " ", Date, " RR variation.pdf"), 8, 8, pointsize = 18)
      par(mfrow = c(2,2), bty = "l", cex.axis = 1.5, cex.lab = 1.5, pch = 19, cex = 0.5, mar = c(5,5,2,2))
      plot(ecgint$RR, ylab = "RR duration (ms)", xlab = "Beat", col = rainbow(max(ecgdata$BeatR), start = 0, end = 2/6))
      abline(h = boxplot.stats(ecgint$RR)$stats, col = 8, lty = c(3,2,1,2,3))
      
      plot(ecgint$RR[-1] ~ ecgint$RR[-nrow(ecgint)], ylab = expression("RR"["n+1"]* " (ms)"), xlab = expression("RR"["n"]* " (ms)"), col = rainbow(max(ecgint$Beat), start = 0, end = 2/6))
      abline(lm(ecgint$RR[-1] ~ ecgint$RR[-nrow(ecgint)]))
      abline(h = boxplot.stats(ecgint$RR)$stats, col = 8, lty = c(3,2,1,2,3))
      abline(v = boxplot.stats(ecgint$RR)$stats, col = 8, lty = c(3,2,1,2,3))

      plot(RRvar, ylab = expression("RR"["n+1"]* " - RR"["n"]* " (ms)"), xlab = "Beat", col = rainbow(max(ecgint$Beat), start = 0, end = 2/6))
      abline(h = boxplot.stats(RRvar)$stats, col = 8, lty = c(3,2,1,2,3))
      
      plot(RRvar[-1] ~ RRvar[-length(RRvar)], ylab = expression("RRvar"["n+1"]* " (ms)"), xlab = expression("RRvar"["n"]* " (ms)"), col = rainbow(max(ecgint$Beat), start = 0, end = 2/6))
      abline(lm(RRvar[-1] ~ RRvar[-length(RRvar)]))
      abline(h = boxplot.stats(RRvar)$stats, col = 8, lty = c(3,2,1,2,3))
      abline(v = boxplot.stats(RRvar)$stats, col = 8, lty = c(3,2,1,2,3))
      dev.off()
    
      pdf(paste0(route, Date, "/", Id, " ", Date, " summary.pdf"), 7, 7, pointsize = 18)
      par(bty = "l", cex.axis = 1.5, cex.lab = 1.5, pch = 19, cex = 0.5, mar = c(5,5,2,2))
      plot(NULL, ylab = paste0("Voltage (mV)"), xlab = paste0("Time - R (ms)"), xlim = summary(medianvoltage$TimeR)[c(1,6)], ylim = c(min(medianvoltage$MV, na.rm = T) - 0.3, max(medianvoltage$MV, na.rm = T) + 0.3))
      abline(h = 0, col = 8)
      polygon(c(medianvoltage$TimeR,rev(medianvoltage$TimeR)),c(medianvoltage$Q4V,rev(medianvoltage$Q0V)), col = adjustcolor(1, alpha.f = 0.2), border = NA)
      lines(medianvoltage$TimeR, medianvoltage$MV, lwd = 2)
      points(medianvoltage[medianvoltage$Wave != "","TimeR"], medianvoltage[medianvoltage$Wave != "","MV"], col = 2)
      text(medianvoltage[medianvoltage$Wave != "","TimeR"], medianvoltage[medianvoltage$Wave != "","MV"], labels = sapply(1:length(medianvoltage[medianvoltage$Wave != "","Wave"]), function(x) wavesplot[which(medianvoltage[medianvoltage$Wave != "","Wave"][x] == waves)]), pos = ifelse(sapply(1:length(medianvoltage[medianvoltage$Wave != "","Wave"]), function(x) which(medianvoltage[medianvoltage$Wave != "","Wave"][x] == waves)) %% 2 == 0, 3, 1), offset = 1, cex = 1.5)
      axis(side = 3, line = -3, at = c(medianvoltage[medianvoltage$Wave == "P","TimeR"],medianvoltage[medianvoltage$Wave == "Q0","TimeR"]) + c(1,-1), labels = F, tck = 0)
      mtext(side = 3, line = -3, at = (medianvoltage[medianvoltage$Wave == "P","TimeR"] + medianvoltage[medianvoltage$Wave == "Q0","TimeR"])/2, bquote(paste("PQ = ", .(round(medecgint$`PQ (ms)`,1)))), cex = 0.7)
      axis(side = 3, line = -3, at = c(medianvoltage[medianvoltage$Wave == "Q0","TimeR"],medianvoltage[medianvoltage$Wave == "T","TimeR"]) + c(1,-1), labels = F, tck = 0)
      mtext(side = 3, line = -3, at = (medianvoltage[medianvoltage$Wave == "Q0","TimeR"] + medianvoltage[medianvoltage$Wave == "T","TimeR"])/2, bquote(paste("QT = ", .(round(medecgint$`QT (ms)`,1)))), cex = 0.7)
      axis(side = 3, line = -1, at = c(medianvoltage[medianvoltage$Wave == "Q0","TimeR"],medianvoltage[medianvoltage$Wave == "T90","TimeR"]) + c(1,-1), labels = F, tck = 0)
      mtext(side = 3, line = -1, at = (medianvoltage[medianvoltage$Wave == "Q0","TimeR"] + medianvoltage[medianvoltage$Wave == "T90","TimeR"])/2, bquote(paste(Q, T[90], c, " = ", .(round(medecgint$`QT90c (ms)`,1)))), cex = 0.7)
      axis(side = 1, line = -4, at = c(medianvoltage[medianvoltage$Wave == "P0","TimeR"],medianvoltage[medianvoltage$Wave == "Pe","TimeR"]) + c(1,-1), labels = F, tck = 0)
      mtext(side = 1, line = -3.5, at = (medianvoltage[medianvoltage$Wave == "P0","TimeR"] + medianvoltage[medianvoltage$Wave == "Pe","TimeR"])/2, bquote(paste(P[0], P[e], " = ", .(round(medecgint$`P0Pe (ms)`,1)))), cex = 0.7)
      axis(side = 1, line = -4, at = c(medianvoltage[medianvoltage$Wave == "Q0","TimeR"],medianvoltage[medianvoltage$Wave == "S","TimeR"]) + c(1,-1), labels = F, tck = 0)
      mtext(side = 1, line = -3.5, at = (medianvoltage[medianvoltage$Wave == "Q0","TimeR"] + medianvoltage[medianvoltage$Wave == "S","TimeR"])/2, bquote(paste("QRS = ", .(round(medecgint$`QRS (ms)`,1)))), cex = 0.7)
      axis(side = 1, line = -2, at = c(medianvoltage[medianvoltage$Wave == "Q0","TimeR"],medianvoltage[medianvoltage$Wave == "Te","TimeR"]) + c(1,-1), labels = F, tck = 0)
      mtext(side = 1, line = -1.5, at = (medianvoltage[medianvoltage$Wave == "Q0","TimeR"] + medianvoltage[medianvoltage$Wave == "Te","TimeR"])/2, bquote(paste(Q, T[e], " = ", .(round(medecgint$`QTe (ms)`,1)))), cex = 0.7)
      mtext(side = 3, line = - 1, at = (medianvoltage[medianvoltage$Wave == "Tn","TimeR"] + medianvoltage[medianvoltage$Wave == "Te","TimeR"])/2, paste0("RR = ", round(medecgint$`RR (ms)`,1)), cex = 0.7)
      mtext(side = 3, line = - 3, at = (medianvoltage[medianvoltage$Wave == "Tn","TimeR"] + medianvoltage[medianvoltage$Wave == "Te","TimeR"])/2, paste0("HR = ", round(medecgint$`HR (bpm)`,1), " bpm"), cex = 0.7)
      dev.off()
    
      ## Relevant results files
      fwrite(ecgwaves, file = paste0(route, Date, "/", Id, " ", Date, " ECG waves.txt"), col.names = T, sep = "\t", append = F, eol = "\r\n")
      fwrite(ecgwaves[ecgwaves$Wave == "R",], file = paste0(route, Date, "/", Id, " ", Date, " ECG R waves.txt"), col.names = T, sep = "\t", append = F, eol = "\r\n")
      fwrite(cbind.data.frame(Id, "Date" = Date, medianvoltage), file = paste0(route, Date, "/", Id, " ", Date, " ECG median wave.txt"), col.names = T, sep = "\t", append = F, eol = "\r\n")
      if (!file.exists(paste0(route, "ECG intervals.txt"))) fwrite(cbind.data.frame(Id, "Date" = Date, ecgint)[1,][-1,], file = paste0(route, "ECG intervals.txt"), row.names = F, col.names = T, sep = "\t", append = F, eol = "\r\n")
      fwrite(cbind.data.frame(Id, "Date" = Date, ecgint), file = paste0(route, "ECG intervals.txt"), row.names = F, col.names = F, sep = "\t", append = T, eol = "\r\n")
      if (!file.exists(paste0(route, "rawdata.txt"))) fwrite(cbind.data.frame(Id, "Date" = Date, medecgint)[-1,], file = paste0(route, "rawdata.txt"), row.names = F, col.names = T, sep = "\t", append = F, eol = "\r\n")
      fwrite(cbind.data.frame(Id, "Date" = Date, medecgint), file = paste0(route, "rawdata.txt"), row.names = F, col.names = F, sep = "\t", append = T, eol = "\r\n")

      if (!file.exists(paste0(route, "processed.log.txt"))) fwrite(cbind.data.frame("Date", "Id", "Register", "Processing.time"), file = paste0(route, "processed.log.txt"), row.names = F, col.names = F, sep = "\t", append = F, eol = "\r\n")
      fwrite(cbind.data.frame(Date, Id, paste(Date, Id), Sys.time()), file = paste0(route, "processed.log.txt"), row.names = F, col.names = F, sep = "\t", append = T, eol = "\r\n")
      setTxtProgressBar(txtProgressBar(style = 3), which(ecglist[,1] == Date & toupper(ecglist[,2]) == Id)/nrow(ecglist))
    }, error = function (e) cat("ERROR:", Date, " ", Id, " ", conditionMessage(e), "\n"))
  }
}

time1 = (proc.time() - time0)[[3]]
paste0(round(round(round(time1/60)/60)/24), "d ", round(round(time1/60)/60)%%24, "h ", round(time1/60)%%60, "min ", round(time1%%60), "s")


