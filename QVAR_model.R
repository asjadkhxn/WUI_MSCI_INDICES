### REFERENCES
#Tomohiro Ando, Matthew Greenwood-Nimmo, Yongcheol Shin (2022) 
#Quantile Connectedness: Modeling Tail Behavior in the Topology of Financial Networks. Management Science 68(4):2401-2431.
#Available at SSRN: https://doi.org/10.1287/mnsc.2021.3984

library("openxlsx")
library("parallel")

#—————————————————————————0.05—————————————————————————————————————
options(mc.cores=detectCores())
#setwd("D:\\Desktop")
source("C:\\Users\\asjad\\OneDrive\\Desktop\\aswini sop\\1\\functionsQVAR.R")

DATA = read.xlsx('C:\\Users\\asjad\\OneDrive\\Desktop\\Main_Data.xlsx',sheet = 'Main')
DATE <- as.Date(DATA[, 1], format = '%Y-%m-%d')
Y = DATA[,-1]
k = ncol(Y)
NAMES = colnames(Y)

### STATIC CONNECTEDNESS APPROACH
tau = 0.05
nlag = 1 # VAR(1)
nfore = 10 # 10-step ahead forecast
qvar_full = QVAR(Y, p=nlag, tau=tau)
CV_full = GFEVD(qvar_full$B, qvar_full$Q, n.ahead=nfore)$GFEVD
rownames(CV_full)=colnames(CV_full)=NAMES
print(DCA(CV_full)$TABLE)

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
space = 200 + nlag # 200 days rolling window estimation
t0 = t-space

total = matrix(NA, ncol=1, nrow=t0)
gfevd = ct = npso = array(NA, c(k, k, t0))
net = from = to = matrix(NA, ncol=k, nrow=t0)
colnames(gfevd)=rownames(gfevd)=colnames(ct)=rownames(ct)=colnames(Y)
for (i in 1:(t-space)){
  qvar = QVAR(Y[i:(space+i-1),], p=nlag, tau=tau)
  gfevd[,,i] = GFEVD(qvar$B, qvar$Q, n.ahead=nfore)$GFEVD
  vd = DCA(gfevd[,,i])
  total[i,] = vd$TCI
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
}
### STORE THE ESTIMATES
#dir.create("D:///Quantile connectedness_codes/tau_v0.05")
#setwd("D://Quantile connectedness_codes/tau_v0.05") 
write.csv(total,'tau_v0.05.csv') 
### END

#—————————————————————————0.95—————————————————————————————————————

#Run this code to get tau_v0.95
options(mc.cores=detectCores())
#setwd("D:\\Desktop")
source("C:\\Users\\asjad\\OneDrive\\Desktop\\aswini sop\\1\\functionsQVAR.R")

DATA = read.xlsx('C:\\Users\\asjad\\OneDrive\\Desktop\\Main_Data.xlsx',sheet = 'Main')
DATE = as.Date(as.character(DATA[,1]),format='%d-%Y-%m')
Y = DATA[,-1]
k = ncol(Y)
NAMES = colnames(Y)

### STATIC CONNECTEDNESS APPROACH
tau = 0.95
nlag = 1 # VAR(1)
nfore = 10 # 10-step ahead forecast
qvar_full = QVAR(Y, p=nlag, tau=tau)
CV_full = GFEVD(qvar_full$B, qvar_full$Q, n.ahead=nfore)$GFEVD
rownames(CV_full)=colnames(CV_full)=NAMES
print(DCA(CV_full)$TABLE)


### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
space = 200 + nlag # 200 days rolling window estimation
t0 = t-space

total = matrix(NA, ncol=1, nrow=t0)
gfevd = ct = npso = array(NA, c(k, k, t0))
net = from = to = matrix(NA, ncol=k, nrow=t0)
colnames(gfevd)=rownames(gfevd)=colnames(ct)=rownames(ct)=colnames(Y)
for (i in 1:(t-space)){
  qvar = QVAR(Y[i:(space+i-1),], p=nlag, tau=tau)
  gfevd[,,i] = GFEVD(qvar$B, qvar$Q, n.ahead=nfore)$GFEVD
  vd = DCA(gfevd[,,i])
  total[i,] = vd$TCI
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
}

write.csv(total,'tau_v0.95.csv') 

#--------------------------fig. 1----------------------
library(openxlsx)
library(parallel)
library(ggplot2)

DATA <- read.xlsx('C:\\Users\\asjad\\OneDrive\\Desktop\\Main_Data.xlsx',sheet = 'Main')
DATE <- as.Date(DATA[201:nrow(DATA), 1], format = '%Y-%m-%d')

DATA_0.05 <- read.csv('tau_v0.05.csv', header = TRUE)
DATA_0.95 <- read.csv('tau_v0.95.csv', header = TRUE)


differences <- DATA_0.95 - DATA_0.05

if (is.vector(differences)) {
  differences <- data.frame(Value = logdifferences)
}

names(differences)[names(differences) == "V1"] <- "Value"
names(differences)[names(differences) == "X"] <- "Date"
differences_df <- data.frame(Date = DATE[-1], Value = differences[-1])



ggplot(differences_df, aes(x = Date, y = Value)) +
  geom_line() +
  labs(title = "Energy market vulnerability index",
       x = "Date",
       y = "Vulnerability of energy markets (%)") +
  theme_minimal()
