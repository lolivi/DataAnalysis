#-----------------------------------------------------
#---------------------ESERCIZIO 1---------------------
#-----------------------------------------------------

'''
Data contained in "grb_data.txt" are related to the afterglow in the X-ray of a gamma-ray burst, as seen by the Swift satellite.
Provide a model to data.
'''

#lmtest per test statistici
library("lmtest")
library("MASS")

grb = read.table(file = "C:/Users/leona/Documenti/Università Fisica/Master MPM Aviation/Data Analysis/Progetto Esame - Analisi Dati/grb_data.txt", header=TRUE, fill=TRUE)
#grb = read.table(file.choose(), header=TRUE, fill=TRUE)
time = grb[,1]
gamma = grb[,2]
title = "Gamma Ray Bursts from Swift Satellite"
xaxis = "Time [s]"
yaxis = "Flux [10^-11 erg/cm2/s]"

#no correlation withouth log
plot(time,gamma,main = title, xlab = xaxis, ylab = yaxis)
cor(grb)
cor.test(time, gamma, alternative="two.sided", conf.level=0.95)

#trying with log
plot(time,gamma,main = title, xlab = xaxis, ylab = yaxis, log = "xy")
#plot(timelog,gammalog,main = title, xlab = xaxis, ylab = yaxis)
timelog = log10(time)
gammalog = log10(gamma)
cor(timelog,gammalog) #almost 1
cor.test(timelog, gammalog, alternative="two.sided", conf.level=0.95)

#---------------
#MODELLO LINEARE
#---------------

lm_full = lm(gammalog ~ timelog)
summary(lm_full) #ci fa un test su coeff ang e intercetta e testa se sono 0 -> Non lo sono perché p-value è 10^-16
#ci sono poi anche r^2 a adj r^2, circa 1 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante
abline(lm_full, col = "red") #printa nello stesso plot una linea rossa

#------------------------------
#STUDIO RESIDUI MODELLO LINEARE
#------------------------------

res<- residuals(lm_full)

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #non sono normali

#2) Test Correlazione

#autocorrelazione
acf(res, lag.max = 100, main = "Autocorrelation - Linear Model")
#le righe blu definiscono gli intervalli di confidenza con distro normale ma distro NON è normale!
#dovrei restringerli ma in ogni caso non sono affidabili

modello <- formula(lm_full)
dwtest(modello, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(modello, order = 40)
#sono correlati

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(lm_full$fitted.values, res, xlab = "Fitted Values of Model", ylab = "Residuals", main = "Residual Plot - Linear Model") 
abline(0,0, col="red") 
bptest(modello) #non potrei usarla perché devo avere residui non correlati e normali

#------------------
#MODELLO PARABOLICO
#------------------

timelog2 <- timelog^2
par_full = lm(gammalog ~ poly(timelog,2))
#par_full = lm(gammalog ~ poly(timelog,5))

summary(par_full) #migliora r^2 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante
ord = order(timelog)
plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],par_full$fit[ord],col="green")
lines(timelog[ord],lm_full$fit[ord],col="red")

#---------------------------------
#STUDIO RESIDUI MODELLO PARABOLICO
#---------------------------------

res<- residuals(par_full)

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #sono normali ma col 16% di p-value

#2) Test Correlazione

#autocorrelazione
acf(res, main = "Autocorrelation - Quadratic Model", lag.max = 100)
#le righe blu definiscono gli intervalli di confidenza con distro normale

modello <- formula(par_full)
dwtest(modello, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(modello, order = 34)
#sono correlati!

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(par_full$fitted.values, res, main = "Residual Plot - Quadratic Model", xlab = "Fitted Values of Model", ylab = "Residuals") 
abline(0,0, col="red") 
bptest(modello) #non potrei usarla perché devo avere residui non correlati e normali

#-------------------------------
#STUDIO OUTLIERS MODELLO LINEARE
#-------------------------------

#62 e 63
plot(lm_full,4) #cook's distance
abline(0.06,0,col="red")
plot(lm_full,5) #leverage

outliers = c(50,51,52,53,54,55,56,57,58,59,60,61,62,63)
time_out =  timelog[outliers] #isola gli indici degli outlier
gamma_out = gammalog[outliers] 

index <- 1 : length(timelog)
clean_idx <- index[!index %in% outliers]
time_clean = timelog[clean_idx]
gamma_clean = gammalog[clean_idx]

#plot outliers
plot(timelog,gammalog,main = title, xlab = xaxis, ylab = yaxis)
lines(time_out, gamma_out, col='red', type = 'p') #li disegno con mod "p" che sta per punti
lines(time_clean, gamma_clean, col='blue', type = 'p')

#--------------------------------
#MODELLO LINEARE SENZA OUTLIER
#--------------------------------

lm_clean = lm(gamma_clean ~ time_clean)
summary(lm_clean) #ci fa un test su coeff ang e intercetta e testa se sono 0 -> Non lo sono perché p-value è 10^-16
#ci sono poi anche r^2 a adj r^2, circa 1 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante

ord_clean = order(time_clean)
#plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
#lines(timelog[ord],lm_full$fit[ord],col="red")
#lines(time_clean[ord_clean],lm_clean$fit[ord_clean],col="blue")

xaxis = "Log Time [s]"
yaxis = "Log Flux [10^-11 erg/cm2/s]"

plot(time_clean[ord],gamma_clean[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],lm_full$fit[ord],col="red")
lines(time_clean[ord_clean],lm_clean$fit[ord_clean],col="blue")

plot(lm_clean,4) #new outliers?
plot(lm_clean,5)

#------------------------------------------
#STUDIO RESIDUI MODELLO LIN SENZA OUTLIER
#------------------------------------------

res<- residuals(lm_clean)

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #sono normali

#2) Test Correlazione

#autocorrelazione
acf(res, lag.max = 100, main = "Autocorrelation - Linear Model (t < 5000 s)")
#le righe blu definiscono gli intervalli di confidenza con distro normale ma distro NON è normale!
#dovrei restringerli ma in ogni caso non sono affidabili

modello <- formula(lm_clean)
dwtest(modello, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(modello, order = 15)
#non sono correlati

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(lm_clean$fitted.values, res, , main = "Residual Plot - Linear Model (t < 5000 s)", xlab = "Fitted Values of Model", ylab = "Residuals") 
abline(0,0, col="red") 
bptest(modello) #residui ok

#----------------------------------
#STUDIO OUTLIERS MODELLO PARABOLICO
#----------------------------------

#61 e 60
plot(par_full,4) #cook's distance
#abline(0.1,0,col = "red")
plot(par_full,5) #leverage

#sequenza togliendo 3 alla volta
#outliers = c(58,60,61)
#outliers = c(50,58,60,61,62,63)
#outliers = c(1,50,52,56,58,60,61,62,63) #dw = 1.16 funziona bene
#outliers = c(1,2,50,52,53,56,57,58,60,61,62,63)
#outliers = c(1,2,3,50,52,53,56,57,58,60,61,62,63)
#outliers = c(1,2,3,50,52,53,56,57,58,59,60,61,62,63)
#outliers = c(1,2,3,4,50,51,52,53,54,56,57,58,59,60,61,62,63)
#outliers = c(1,2,3,4,50,51,52,53,54,55,56,57,58,59,60,61,62,63)

#sequenza togliendo all'inizio tutti quelli sopra 0.06
#outliers = c(50,52,56,58,60,61,62,63) #funziona! 1.05
#outliers = c(1,50,52,53,56,58,60,61,62,63)
#outliers = c(1,2,3,50,52,53,56,58,60,61,62,63)
#outliers = c(1,2,3,4,50,52,51,53,56,57,58,60,61,62,63)

#sequenza togliendo all'inizio tutti quelli sopra 0.1
#outliers = c(50,58,60,61) #funziona! 0.96 di dw
#outliers = c(50,52,58,60,61,62,63) #1.14
#outliers = c(1,50,52,53,56,58,60,61,62,63)
#outliers = c(1,2,50,52,53,56,58,60,61,62,63)
#outliers = c(1,2,3,50,52,53,56,58,60,61,62,63)
outliers = c(1,50,52,53,57,58,60,61,62,63) #1.3 di dw

time_out =  timelog[outliers] #isola gli indici degli outlier
gamma_out = gammalog[outliers] 

index <- 1 : length(timelog)
clean_idx <- index[!index %in% outliers]
time_clean = timelog[clean_idx]
gamma_clean = gammalog[clean_idx]

#plot outliers
plot(timelog,gammalog,main = title, xlab = xaxis, ylab = yaxis)
lines(time_out, gamma_out, col='red', type = 'p') #li disegno con mod "p" che sta per punti
lines(time_clean, gamma_clean, col='blue', type = 'p')

#--------------------------------
#MODELLO PARABOLICO SENZA OUTLIER
#--------------------------------

time_clean2 <- time_clean^2
par_clean = lm(gamma_clean ~ poly(time_clean,2))
#par_clean = lm(gamma_clean ~ poly(time_clean,5))
summary(par_clean) #ci fa un test su coeff ang e intercetta e testa se sono 0 -> Non lo sono perché p-value è 10^-16
#ci sono poi anche r^2 a adj r^2, circa 1 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante

ord_clean = order(time_clean)
plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],par_full$fit[ord],col="red")
lines(time_clean[ord_clean],par_clean$fit[ord_clean],col="blue")

plot(par_clean,4) #new outliers?
plot(par_clean,5)

#------------------------------------------
#STUDIO RESIDUI MODELLO PAR SENZA OUTLIER
#------------------------------------------

res<- residuals(par_clean)

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #sono normali

#2) Test Correlazione

#autocorrelazione
acf(res, lag.max = 100, main = "Autocorrelation - Quadratic Model (No Outliers)")
#le righe blu definiscono gli intervalli di confidenza con distro normale ma distro NON è normale!
#dovrei restringerli ma in ogni caso non sono affidabili

modello <- formula(par_clean)
dwtest(modello, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(modello, order = 7)
#sono correlati

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(par_clean$fitted.values, res, xlab = "Fitted Values of Model", ylab = "Residuals", main = "Residual Plot - Quadratic Model (No Outliers)") 
abline(0,0, col="red") 
bptest(modello) #non potrei usarla perché devo avere residui non correlati e normali

#---------------------------
#MODELLO PARABOLICO CON MASS
#---------------------------

timelog2 <- timelog^2
par_full_M = rlm(gammalog ~ poly(timelog,2), method = "M")
summary(par_full_M) #migliora r^2 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante

ord = order(timelog)
plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],par_full$fit[ord],col="green")
lines(timelog[ord],par_full_M$fit[ord],col="red")

#--------------------------------------
#STUDIO RESIDUI MODELLO PARABOLICO MASS
#--------------------------------------

res<- par_full_M$residuals

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #sono normali ma col 16% di p-value

#2) Test Correlazione

#autocorrelazione
acf(res, lag.max = 100, main = "Autocorrelation - Quadratic Model (M-Estimator)")
#le righe blu definiscono gli intervalli di confidenza con distro normale

modello <- formula(par_full_M)
dwtest(modello, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(modello, order = 45)
#sono correlati!

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(par_full_M$fitted.values, res,  xlab = "Fitted Values of Model", ylab = "Residuals", main = "Residual Plot - Quadratic Model (M-Estimator)") 
abline(0,0, col="red") 
bptest(modello) #non potrei usarla perché devo avere residui non correlati e norma

#---------------------------
#MODELLO PARABOLICO CON LTS
#---------------------------

timelog2 <- timelog^2
par_full_lqs = lqs(gammalog ~ poly(timelog,2), method = "lts")
summary(par_full_lqs) #migliora r^2 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante

ord = order(timelog)
plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],par_full$fit[ord],col="green")
lines(timelog[ord],par_full_lqs$fit[ord],col="red")

#--------------------------------------
#STUDIO RESIDUI MODELLO PARABOLICO LQS
#--------------------------------------

res<- par_full_lqs$residuals

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #sono normali ma col 16% di p-value

#2) Test Correlazione

#autocorrelazione
acf(res)
#le righe blu definiscono gli intervalli di confidenza con distro normale

modello <- formula(par_full_lqs)
dwtest(modello, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(modello, order = 10)
#sono correlati!

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(par_full_lqs$fitted.values, res) 
abline(0,0, col="red") 
bptest(modello) #non potrei usarla perché devo avere residui non correlati e norma

#----------------------------
#MODELLO CON LOESS
#----------------------------

help(loess)
loess_fit = loess(gammalog ~ timelog, degree = 2)
summary(loess_fit) #migliora r^2 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante

ord = order(timelog)
ord = 50:63
plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],lm_full$fit[ord],col="blue")
lines(timelog[ord],par_full$fit[ord],col="green")
lines(timelog[ord],loess_fit$fit[ord],col="red")

ord = 1:49
plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],lm_clean$fit[ord],col="blue")
lines(timelog[ord],par_full$fit[ord],col="green")
lines(timelog[ord],loess_fit$fit[ord],col="red")

#---------------------------------------
#STUDIO RESIDUI MODELLO LOESS
#---------------------------------------

res<- loess_fit$residuals
res = residuals(loess_fit)

dw_den = 0.
dw_num = 0.

for (i in 1:63) {
  dw_den = dw_den + res[i]*res[i]
}

for (i in 2:63) {
  dw_num = dw_num + (res[i]-res[i-1])^2
}

print(dw_num)
print(dw_den)
dw = dw_num/dw_den
print(dw)

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #sono normali ma col 16% di p-value

#2) Test Correlazione

#autocorrelazione
acf(res,lag.max = 100,main = "Autocorrelation - Loess Model")
#le righe blu definiscono gli intervalli di confidenza con distro normale

modello <- formula(loess_fit)
dwtest(loess_fit, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(loess_fit, order = 45)
#sono correlati!

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(loess_fit$fit, res, xlab = "Fitted Values of Model", ylab = "Residuals", main = "Residual Plot - Loess Model") 
abline(0,0, col="red") 
bptest(modello) #non potrei usarla perché devo avere residui non correlati e norma

#-------------------------------
#STUDIO OUTLIERS MODELLO LOESS
#-------------------------------

#62 e 63
plot(loess_fit,4) #cook's distance
abline(0.06,0,col="red")
plot(loess_fit,5) #leverage

outliers = c(58,59,60,61,62,63)
time_out =  timelog[outliers] #isola gli indici degli outlier
gamma_out = gammalog[outliers] 

index <- 1 : length(timelog)
clean_idx <- index[!index %in% outliers]
time_clean = timelog[clean_idx]
gamma_clean = gammalog[clean_idx]

#plot outliers
plot(timelog,gammalog,main = title, xlab = xaxis, ylab = yaxis)
lines(time_out, gamma_out, col='red', type = 'p') #li disegno con mod "p" che sta per punti
lines(time_clean, gamma_clean, col='blue', type = 'p')


#--------------------------------
#MODELLO LOESS SENZA OUTLIER
#--------------------------------

loess_clean = loess(gammalog ~ timelog, degree = 2)
summary(loess_clean) #ci fa un test su coeff ang e intercetta e testa se sono 0 -> Non lo sono perché p-value è 10^-16
#ci sono poi anche r^2 a adj r^2, circa 1 
#F-stat che ha come ipotesi nulla che i dati sono descritti da costante

ord_clean = order(time_clean)
plot(timelog[ord],gammalog[ord],main = title, xlab = xaxis, ylab = yaxis)
lines(timelog[ord],loess_fit$fit[ord],col="red")
lines(time_clean[ord_clean],loess_clean$fit[ord_clean],col="blue")

plot(loess_clean,4) #new outliers?
plot(loess_clean,5)

#------------------------------------------
#STUDIO RESIDUI MODELLO PAR SENZA OUTLIER
#------------------------------------------

res<- residuals(loess_clean)

#1) Test Normalità Residui

#controlliamo se i residui sono normali
qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals quantiles")
qqline(res, distribution=qnorm, col="red") #aggiungiamo la retta di riferimento
shapiro.test(res) #sono normali

#2) Test Correlazione

#autocorrelazione
acf(res,lag.max = 100)
#le righe blu definiscono gli intervalli di confidenza con distro normale ma distro NON è normale!
#dovrei restringerli ma in ogni caso non sono affidabili

modello <- formula(loess_clean)
dwtest(modello, alternative="two.sided") #ipotesi nulla -> autocorrelazione è 0
bgtest(modello, order = 10)
#sono correlati

#residual plot res vs valori stimati per capire se è costante varianza o no
plot(loess_clean$fitted.values, res) 
abline(0,0, col="red") 
bptest(modello) #non potrei usarla perché devo avere residui non correlati e normali

