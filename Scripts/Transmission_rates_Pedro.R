library(readxl)
library(data.table)

ppv <- 0.974
npv <- 0.998

dt <- as.data.table(read_excel("./Data/FarmData.xlsx", sheet="Model", skip =3, col_names =FALSE, range ="A4:AB1000")) 

hdr <- c("DamID", "Farm1", "D_Pre_E", "D_Pre_PCR", "D_Post_E", "D_Post_PCR", 
         "Sample Type", "Age_Class", "Age_Months", "D_Pre_Date", "D_Post_Date", 
         "Preg_Outcome", "Calf_ID", "Farm2", 
         "C_Birth_E", "C_Birth_PCR", "C_Post_E", "C_Post_PCR", "Sample_Type_C", "Calf_Sex", 
         "Calf_Pre_Date", "Calf_Post_Date",
         "Calf_Birth_Date", "Colostrum_ID", "Farm3", "Colostrum_E", "Colostrum_PCR", "Notes")

stopifnot(ncol(dt) >= length(hdr))
setnames(dt, hdr)

#### delta times ####
date_cols <- c("D_Pre_Date", "D_Post_Date", "Calf_Pre_Date", "Calf_Post_Date", "Calf_Birth_Date")

# Convert and clean bad entries
dt[, (date_cols) := lapply(.SD, function(x) {
  x <- as.character(x)
  x[!grepl("^\\d{1,2}[-/\\.]\\d{1,2}[-/\\.]\\d{2,4}|\\d{4}[-/\\.]\\d{1,2}[-/\\.]\\d{1,2}$", x)] <- NA
  as.Date(x, tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y"))
}), .SDcols = date_cols]

dt[, D_Delta_Dates := as.numeric(D_Post_Date - D_Pre_Date)]
dt[, Calf_Delta_Dates := as.numeric(Calf_Pre_Date - Calf_Birth_Date)]
dt[, Heif_Delta_Dates := as.numeric(Calf_Post_Date - Calf_Pre_Date)]

png("calf_sampling_histogram.png", width = 800, height = 600)
hist(dt[, Calf_Delta_Dates], main = "Calf Sampling Interval", xlab = "Days")
dev.off()

dt[, unique(D_Delta_Dates)]
dt[, unique(Calf_Delta_Dates)]
dt[, unique(Heif_Delta_Dates)]

#### Prevalence Cow pre and PTcow ####
Pos <- dt[D_Pre_E=="Pos" , .(Pos=.N), by=Farm1]
ntested <- dt[!is.na(D_Pre_E) & D_Pre_E!="Doub" , .(n=.N), by=Farm1]
Prev_cowPre <- merge(ntested, Pos, by="Farm1", all.x=TRUE)

Prev_cowPre[is.na(Pos), Pos :=0]
Prev_cowPre[, Prev:=Pos/n]
print(Prev_cowPre)
print(Prev_cowPre[1:4, mean(n)])
print(Prev_cowPre[1:4, mean(Prev)])

cows <- dt[!is.na(D_Pre_E) & !is.na(D_Post_E),]
print(cows[, table(D_Pre_E, D_Post_E, Farm1)])
cont_table <- cows[, table(D_Pre_E, D_Post_E, Farm1)]
# D_Post_E
# D_Pre_E Doub Neg Pos
# Doub    2   1   3
# Neg    10 100  41
# Pos     2   5  50

# raw ht
41/141
# 0.2907801
# Note cells a,b,c,d are for columns/rows with results in Pos and then Neg order
#     Pos   Neg
# Pos   a     b
# Neg   c     d

# Numerators Cows
a <- c(50, 36, 5, 5, 1, 3)
b <- c(5, 0, 4, 1, 0, 0)
c <- c(41, 22, 14, 1, 0, 4)
d <- c(100, 45, 26, 8, 17, 4)

# Denominator Cows
neg_pre <- c(141, 67, 40, 9, 17, 8)
pos_pre <- c(55, 36, 9, 6, 1, 3)

# Calculate numerator and denominator
numerat <- a*(1-ppv)*ppv + b*(1-ppv)*npv + c*npv*ppv + d*npv*(1-npv)
denom <- neg_pre*npv + pos_pre*(1-ppv)

# PT cows
PTcows <- numerat / denom
names(PTcows) <- c("TOTAL", "PH", "G", "NM", "W", "NC")
print(PTcows)

# cows P(HT) 
PTcows <- numerat/denom
# [1] 0.2915952 

#### P(HT) after first week ####
calves <- dt[!is.na(C_Birth_E) & !is.na(C_Post_E),]
print(calves[, table(C_Birth_E, C_Post_E)])
# C_Post_E
# C_Birth_E Doub Neg Pos
# Doub    0   3   0
# Neg     1  28   1
# Pos     0   3   7

# Numerators
a <- c(7, NA, 2, 4, 1, NA)   # Pos–Pos
b <- c(3, NA, 3, 0, 0, NA)   # Pos–Neg
c <- c(1, NA, 0, 1, 0, NA)   # Neg–Pos
d <- c(28, NA, 18, 0, 10, NA) # Neg–Neg

# Denominator
neg_birth <- c(29, NA, 18, 1, 10, NA)
pos_birth <- c(10, NA, 5, 4, 1, NA)

# Compute numerator and denominator
numerat <- a*(1-ppv)*ppv + b*(1-ppv)*npv + c*npv*ppv + d*npv*(1-npv)
denom <- neg_birth*npv + pos_birth*(1-ppv)

# P(HT) for calves
PTcalves <- numerat / denom
names(PTcalves) <- c("TOTAL", "PH", "G", "NM", "W", "NC")
print(PTcalves)

# calves P(HT)
PTcalvesThreeMonths <- numerat/denom
print(PTcalvesThreeMonths)
# [1] 0.04393713

# Compare with P(HT) for calves at birth (~ 4 days)
PTcalvesCol <- 0.172
PTcalves4days <- 1-(1-PTcalvesThreeMonths)^(1/(365/4))
PTcolostrum <- PTcalvesCol - PTcalves4days
print(PTcolostrum)
