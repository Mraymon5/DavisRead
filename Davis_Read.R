#This is a script that should be run as-is, which will then add the function Davis_Read() to the global environment. This is a program that sifts through a folder of .ms8.txt files, and extracts some data from them, saving it to CSVs. MAR 8/1/23

#Settings====
Davis_Read <- function(Folder = NA,  #sets the folder containing the data
                       KeyWord = ".txt", #sets a search term for identifying files
                       File = NA, #A file path to a specific single file to be run
                       Cutoff = 50, #minimum acceptable ILI
                       OutputFolder = F, #sets a folder to save .csv of output
                       ILIFolder = F, #sets a folder to save a .csv of ILIs
                       Title = NULL, #sets a title for data output
                       PlotRawInterval = F, #plot raw interval histogram
                       PlotInterval = F, #plot filtered interval hist.
                       PlotCumulative = F, #plot licks over time
                       Timer = F,
                       StJohnFile = NULL #a .txt file provided by Steven St John, requires that .ms8 be supplied to File
                       ){   #reports runtime

if(0){ #For debugging
  Folder = NA;  #sets the folder containing the data
  KeyWord = ".txt"; #sets a search term for identifying files
  File = "~/Documents/MAR_Data/MR07/BAT/20241223/20241223MR07.txt"; #A file path to a specific single file to be run
  Cutoff = 50; #minimum acceptable ILI
  OutputFolder = F; #sets a folder to save .csv of output
  ILIFolder = F; #sets a folder to save a .csv of ILIs
  Title = NULL; #sets a title for data output
  PlotRawInterval = F; #plot raw interval histogram
  PlotInterval = F; #plot filtered interval hist.
  PlotCumulative = F; #plot licks over time
  Timer = F;
  StJohnFile = NULL #a .txt file provided by Steven St John, requires that .ms8 be supplied to File
}
  
require(tidyverse)

if (Timer != F) {
  ptm <- proc.time()
} # This starts a timer

ILI <- ILIFolder
O   <- OutputFolder

if (is.na(Folder) & is.na(File)){
  stop("You must specify a File or a Folder")
}

if (!is.na(Folder)){
  FileNames <- list.files(path = Folder, pattern = KeyWord, full.names = T) #creates a list of file names with full path
  Directory <- data.frame(Files = list.files(path = Folder, pattern = KeyWord, full.names = F), Date = NA, Animal = NA) #creates a list of file names without the path
}

if (!is.na(File)){
  FileNames <- File #write the path of the specified file over the list of file names, with full path
  Directory <- data.frame(Files = gsub(".*/","",File), Date = NA, Animal = NA) #creates a list of file names without the path
  if (!is.na(Folder)){warning("File overwrites Folder")}
}

isMS8 = grepl(x = Directory$Files, pattern = ".MS8")
if (isMS8){
  Directory$Files <- sub(pattern = ".MS8.txt", replacement = "", x = Directory$Files) #trims the extension(s) from the directory list
  for (i in 1:length(Directory$Files)) { #This loop parses the date and animal ID from the file name
    Directory$Date[i] <- paste(substr(x = Directory$Files[i], start = 1, stop = 2), substr(x = Directory$Files[i], start = 3, stop = 4),sep = "/")
    Directory$Animal[i] <- substr(x = Directory$Files[i], start = 5, stop = 20)
  }
} else {
  Directory$Files <- sub(pattern = ".txt", replacement = "", x = Directory$Files) #trims the extension(s) from the directory list
  for (i in 1:length(Directory$Files)) { # i = 1 #This loop parses the date and animal ID from the file name
    Directory$Date[i] <- paste(substr(x = Directory$Files[i], start = 5, stop = 6), substr(x = Directory$Files[i], start = 7, stop = 8),substr(x = Directory$Files[i], start = 1, stop = 4),sep = "/")
    Directory$Animal[i] <- substr(x = Directory$Files[i], start = 9, stop = 200)
  }
}


for (file_N in 1:length(FileNames)) { #Start whole folder loop
  
  if (!is.na(File)) {
    D <- File #if file has been supplied, store File to D, which is used as a placeholder for files
  } else {
    D <- FileNames[file_N] #if file has not been supplied, pull the appropriate file name from the list
  }
  
  
#Data Import====
VersionAt = read_lines(file = D) %>% grepl(pattern = "Version") %>% which
Version <- read.csv(file = D, skip = VersionAt-1, header = F, sep = ",", as.is = T, nrows = 1) #check the first line of the MS8 output to see if it's the old version or the new(er) version

if (!is.null(StJohnFile)){
  Temp_Lat <- read.table(file = StJohnFile, nrows = 1, skip = 3, sep = " ", strip.white = T) %>% t
  if (Temp_Lat[6] == "s") Temp_Lat = as.numeric(Temp_Lat[5]) * 1000 else Temp_Lat = as.numeric(Temp_Lat[5])
  Licks_Dav <- read.table(file = StJohnFile, skip = 9) %>% t()
  Meta <- read.table(file = File, skip = 5, nrows = 2, sep = ",", strip.white = T, row.names = 1)
  Meta_Data <- read.table(file = File, skip = 8, nrows = Meta[2,1], sep = ",", strip.white = T, col.names = c("PRES","TUBE","CONCENTRATION","SOLUTION",  "IPI"  , "LENGTH", "LICKS"))
  Meta_Data$Latency = Temp_Lat
} else {
  MetaAt = read_lines(file = D) %>% grepl(pattern = "Max Wait") %>% which
  Meta <- read.csv(file = D, sep = ",", as.is = T, header = F, row.names = 1, skip = MetaAt-1, nrows = 3) #in theory, capture the number of trials and the trial time limit
  Meta = Meta[c(1,3),] %>% data.frame
  DataAt = read_lines(file = D) %>% grepl(pattern = "PRESENTATION,TUBE") %>% which
  BetterMeta = read.csv(file = D, nrows = (DataAt-2), skip = 1, row.names = 1, header = F) %>% t
  Meta_Data <- read.csv(file = D, skip = DataAt-1, nrows = Meta[2,1])
  Licks_Dav <- (read.csv(file = D, skip = (DataAt + Meta[2,1]), sep = ",", fill = T, header = F, row.names = 1, flush = T, allowEscapes = T, col.names = c(1:max(Meta_Data$LICKS))))
}


Licks_R <- data.frame(t(x = Licks_Dav), row.names = c(1:ncol(Licks_Dav))) #transpose the licks and change the row names back to regular

for (i in 1:ncol(Licks_R)) { #turn the ILIs into sequential timecodes, where each represents the time from trial start for that lick
  for (is in 1:nrow(Licks_R)) {
    if (is.na(Licks_R[is,i]) == F) {
      if (is == 1) {
        Licks_R[is,i] <- Meta_Data$Latency[i] + Licks_R[is,i]
      } else {
        Licks_R[is,i] <- Licks_R[is,i] + Licks_R[is - 1,i]
      }
    }
    
  }
}

Licks_R <- rbind(Meta_Data$Latency, Licks_R) #append latencies to the lick time intervals
Licks_R[1,][Licks_R[1,] == (Meta[1,1] * 1000)] <- NA #if the latency is the trial time limit, store as NA

#Microstructure Stuff====
for (trial_N in 1:ncol(Licks_R)) {#trial_N = 1 #start per trial loop
  
  TrialLab <- paste(Meta_Data$CONCENTRATION[trial_N], Meta_Data$SOLUTION[trial_N], sep = " ") #title of the trial
  
  if (Meta_Data$LICKS[trial_N] == 0){
    Lick_Data <- data.frame(RawTime = 0)
  } else{
  Lick_Data <- na.omit(data.frame(RawTime = Licks_R[,trial_N]))
  }
Lick_Data[, 2:7] <- 0
  colnames(Lick_Data)[2:7] <- c("FilteredTime", "BurstCount", "RawInterval", "FilteredInterval", "DoubleContact", "LickTotal")
  Lick_Data <- rbind(0, Lick_Data)
  for (X in 2:nrow(Lick_Data)) { 
    if ((Lick_Data$RawTime[X]-Lick_Data[X-1,2])>Cutoff) {
      Lick_Data[X,2] <- Lick_Data[X,1]
    } else {
      Lick_Data[X,2] <- Lick_Data[X-1,2]
    }} #This loop Processes the Filtered Time column
  
    LL <- 0
  for (X in 2:nrow(Lick_Data)) {
    if ((Lick_Data[X, 2] - Lick_Data[X - 1, 2]) > Cutoff) {
      LL <- LL + 1
    }
    Lick_Data[X, 7] <- LL
  } #This loop processes the Lick Total column of Lick_Data
  
  # Burst Analysis, Left ==== 
  DC <- 0
  for(X in 2:nrow(Lick_Data)) {
    if ((Lick_Data[X, 2] - Lick_Data[X - 1, 2]) > 0 & (Lick_Data[X, 2] - Lick_Data[X - 1, 2]) < 1000) {
      Lick_Data[X, 3] <- (Lick_Data[X - 1, 3] + 1)
      Lick_Data[X - 1, 6] <- DC
    } else {
      if ((Lick_Data[X, 2] - Lick_Data[X - 1, 2]) == 0) {
        Lick_Data[X, 3] <- Lick_Data[X - 1, 3]
        DC <- DC + 1; Lick_Data[X - 1, 6] <- DC
      } else {
        Lick_Data[X, 3] <- 1
        Lick_Data[X - 1, 6] <- DC
        DC <- 0
      }}} #This Loop processes the Burst Count column of Lick_Data
  
  Burst_Data <- data.frame(Licks = 0, Duration = 0, Start = 0, Stop = 0)
  N <- 0
  for(X in 2:nrow(Lick_Data)) {
    if ((Lick_Data[X, 3] - Lick_Data[X - 1, 3]) < (-1)) {
      N <- N + 1
      Burst_Data[N, 1] <- Lick_Data[X - 1, 3]
    }
    if (X == nrow(Lick_Data) & Lick_Data[X, 3] > 2){
      N <- N + 1
      Burst_Data[N, 1] <- Lick_Data[X, 3]
    }} #This loop processes Licks in Burst_Data
  
  LB <- N
  N <- 0
  for(X in 2:nrow(Lick_Data)){
    if ((Lick_Data[X, 3] - Lick_Data[X - 1, 3]) < (-1)) {
      N <- N + 1
      Burst_Data[N, 2] <- Lick_Data[X - 1, 2] - Lick_Data[X - (Lick_Data[X - 1, 3] + Lick_Data[X - 1, 6]), 2]
      Burst_Data[N, 4] <- Lick_Data[X - 1, 1]
      Burst_Data[N, 3] <- Lick_Data[(X - (Lick_Data[X - 1, 3] + Lick_Data[X - 1, 6])) ,1]
    }
    if (X == nrow(Lick_Data) & Lick_Data[X, 3] > 2){
      N <- N + 1
      Burst_Data[N, 2] <- Lick_Data[X, 2] - Lick_Data[X - (Lick_Data[X - 1, 3] + Lick_Data[X - 1, 6]), 2]
      Burst_Data[N, 4] <- Lick_Data[X, 1]
      Burst_Data[N, 3] <- Lick_Data[(X - (Lick_Data[X - 1, 3] + Lick_Data[X - 1, 6])) ,1]
    }} #This loop processes Burst_Data[,2:4]
  
  if (Meta_Data$LICKS[trial_N] == 0) {
    PauseL <- data.frame(Start = 0, Stop = (Meta[1,1] * 1000))
  } else {
  PauseL <- data.frame(Start = c(0, Burst_Data$Stop), Stop = c(Burst_Data$Start, Licks_R[1,trial_N]+(Meta_Data$LENGTH[trial_N]*1000)))
  }
  PauseL$Pause <- PauseL$Stop - PauseL$Start
  
  # Interval Analyses, Left ==== 
  Lick_Data[1, 4] <- 0
  for(X in 2:nrow(Lick_Data)) {
    Lick_Data[X, 4] <- (Lick_Data[X, 1] - Lick_Data[X - 1, 1])
  } #Processes Lick_Data$RawInterval
  Lick_Data[1, 5] <- 0
  for(X in 2:nrow(Lick_Data)) {
    Lick_Data[X, 5] <- (Lick_Data[X, 2] - Lick_Data[X - 1, 2])
  }

  #Plotting====
    if (PlotRawInterval == "A"||PlotRawInterval == trial_N) {
    hist(Lick_Data$RawInterval[Lick_Data$RawInterval <= 1000], main = paste(trial_N, TrialLab, sep = " "), xlab = "Interval (mS)", xlim = c(0, 1000), breaks = seq(0, 1000, by = 10))
  }

if (PlotInterval == "A"||PlotInterval == trial_N) {
    hist(Lick_Data$RawInterval[Lick_Data$FilteredInterval <= 1000 & Lick_Data$FilteredInterval > 0], main = paste(trial_N, TrialLab, sep = " "), xlab = "Interval (mS)", xlim = c(0, 1000), breaks = seq(0, 1000, by = 10), col = "green4")
  }
  
  if (PlotCumulative == "A" || PlotCumulative == trial_N) {
    plot(x = ((Lick_Data$RawTime) / 1000), y = c(1:nrow(Lick_Data)), type = "l", xlab = "Seconds", ylab = "Licks", main = paste(trial_N, TrialLab, sep = " "))
  }
  
  #Data Export ====
  Lickout <- data.frame(
    "Trial" = trial_N,
    "Concentration" = Meta_Data$CONCENTRATION[trial_N],
    "Solution" = Meta_Data$SOLUTION[trial_N],
    "Latency" = Meta_Data$Latency[trial_N],
    "LickTotal" = LL,
    "Min.Lick" = max(Lick_Data$LickTotal[Lick_Data$FilteredTime < 60000]),
    "BurstTotal" = LB,
    "BurstLicks" = mean(Burst_Data$Licks),
    "BurstDur." = mean(Burst_Data$Duration),
    "BurstLatency" = Burst_Data[1, 3],
    "PauseTotal" = nrow(PauseL),
    "PauseDur." = mean(PauseL$Pause),
    "MPI" = mean(Lick_Data$FilteredInterval[Lick_Data$FilteredInterval > Cutoff & Lick_Data$FilteredInterval <= 160]),
    "Efficiency" = (NROW(Lick_Data$FilteredInterval[Lick_Data$FilteredInterval > 0 & Lick_Data$FilteredInterval <= 160])) / LL,
    "Date" = Directory$Date[file_N],
    "Animal" = Directory$Animal[file_N],
    "File" = D,
    stringsAsFactors = F)
  
  
  if (trial_N > 1) { 
    LickExport[(nrow(LickExport) + 1), ] <- Lickout[1, ]
  } else {
    LickExport <- Lickout
  }
  
  if (OutputFolder != F) {
    if (is.null(Title)){
      write.csv(LickExport, file = paste(O,Directory$Files[file_N],".csv", sep = ""), row.names = F)
    } else {
      write.csv(LickExport, file = paste(O,Title,file_N,".csv", sep = ""), row.names = F)
    }
  }
  
  ILItabX <- data.frame(Lick_Data$FilteredInterval, trial_N)
  
  colnames(ILItabX) <- c("ILI", "ID")
  
  if (ILIFolder != F) {
    if (trial_N > 1)
    {
      ILItab <- read.csv(paste(ILI,Title,"ILI",file_N,".csv", sep = ""),header = T)
      ILItab <- data.frame(c(ILItab$ILI, ILItabX$ILI), c(ILItab$ID, ILItabX$ID))
      colnames(ILItab) <- c("ILI", "ID")
      write.csv(ILItab,file = paste(ILI,Title,"ILI",file_N,".csv", sep = ""), row.names = F)
      #View(ILItab)
    } else {
      ILItab <- ILItabX
      colnames(ILItab) <- c("ILI","ID")
      write.csv(ILItab,file = paste(ILI,Title,"ILI",file_N,".csv", sep = ""),row.names = F)
    }}
  
  if (is.na(File)) {
    assign(x = paste(Title, "Directory", sep = ""), value = Directory, envir = .GlobalEnv)
  } else {  
    if (is.null(Title)){
      assign(x = Directory$Files[file_N], value = LickExport, envir = .GlobalEnv)
    } else {
      assign(x = Title, value = LickExport, envir = .GlobalEnv)
    }
  }
}
if (!is.na(File)) {
  break
} else {
  print(file_N)
}

}
if (Timer!= F) {
  proc.time() - ptm
}}
