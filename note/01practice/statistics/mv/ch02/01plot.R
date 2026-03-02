# 1. plot
par(mfrow=c(2,2), mai=c(0.6, 0.6, 0.4, 0.4), cex=0.7, cex.main=1)
plot(`PH_M_cm-Rep1`, `PH_M_cm-Rep2`, main = "(a) scatter plot")
plot(as.factor(`COUNTRY 
of origin`), xlab="country", main="(b) bar plot")
plot(`Hd_dto_days-CFLN06`~as.factor(`GrwHabit_E_sw-CFLN06`), xlab="Grow Habit", main="(c) box plot")
plot(as.factor(`COUNTRY 
of origin`)~as.factor(`GrwHabit_E_sw-CFLN06`), main="(d) ridge plot")

# 2. other
par(mfrow=c(2,2), mai=c(0.6,0.6,0.4,0.4), cex=0.6)
data <- data.frame(PH_M_cm_Rep1 = `PH_M_cm-Rep1`, PH_M_cm_Rep2 = `PH_M_cm-Rep2`)
data <- data[!grepl("\\*", data$PH_M_cm_Rep1) & !grepl("\\*", data$PH_M_cm_Rep2), ]
data$PH_M_cm_Rep1 <- sapply(strsplit(as.character(data$PH_M_cm_Rep1), "-"), `[`, 1)
data$PH_M_cm_Rep2 <- sapply(strsplit(as.character(data$PH_M_cm_Rep2), "-"), `[`, 1)
lm.fit <- lm(PH_M_cm_Rep1 ~ PH_M_cm_Rep2, data = data)
plot(lm.fit)









