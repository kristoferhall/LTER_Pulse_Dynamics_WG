# Wavelet and cross-wavelet example script using the WaveletComp R package
#
# See WaveletComp 1.1: A guided tour through the R package
# http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf
#
# Basic examples and discussion come from this document
#
#
#
# First, there is an example of wavelet analysis with synthetic data.
# Second, there is an example of cross-wavelet analysis with synthetic data.
# Third, conduct wavelet and cross-wavelet analyses using Ameriflux 
#   daily precipitation and NEE data from SEV LTER US-Seg tower (desert grassland site)
# Fourth, conduct wavelet and cross-wavelet analyses of FCE LTER monthly water
#   level and dissolved organic carbon (DOC).

library(WaveletComp)
library(dplyr)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(lubridate)
library(readr)


# 1. Wavelet example with synthetic data --------------------------

# series with constant period of 50 and added noise
x <- periodic.series(start.period = 50, length = 1000)
x <- x + 0.2*rnorm(1000) 

my.data <- data.frame(x = x)

my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, 
                        n.sim = 10)
# see guided tour p. 12 for info on the function arguments. 
# dt is number of observations per time period

# analyze.wavelet returns a number of attributes
attributes(my.w)

# wavelet image
wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))

# returns high power around a period of 50, as expected. White bands
# show where power is significantly different than white noise.


# look at the average power of the series
maximum.level = 1.001*max(my.w$Power.avg)
wt.avg(my.w, maximum.level = maximum.level)


# it is possible to reconstruct the original signal from the wavelet transform
reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2),
            legend.coords = "bottomleft", ylim = c(-1.8, 1.8))







# 2. Cross-wavelet example with synthetic data -----------------------

x1 <- periodic.series(start.period = 1*24, length = 24*96)
x2 <- periodic.series(start.period = 2*24, length = 24*96)
x3 <- periodic.series(start.period = 4*24, length = 24*96)
x4 <- periodic.series(start.period = 8*24, length = 24*96)
x5 <- periodic.series(start.period = 16*24, length = 24*96)

x <- x1 + x2 + 3*x3 + x4 + x5 + 0.5*rnorm(24*96)
y <- x1 + x2 - 3*x3 + x4 + 3*x5 + 0.5*rnorm(24*96)


# plot raw data to show differences in periodicity
time <- 1:(24*96)

example_x <- data.frame(time = time, variable = x, name = "x")
example_y <- data.frame(time = time, variable = y, name = "y")

example_data <- rbind(example_x, example_y)

ggplot(example_data) + 
        geom_line(aes(x = time, y = variable)) +
        facet_wrap(~name, nrow = 2)



my.data <- data.frame(x = x, y = y)
my.wc <- analyze.coherency(my.data, my.pair = c("x","y"),
                           loess.span = 0,
                           dt = 1/24, dj = 1/100,
                           lowerPeriod = 1/2,
                           make.pval = TRUE, 
                           n.sim = 10)
# - loess.span = 0 because no need to detrend series
# - dt = 1/24 means 24 observations per time unit
# - lowerPeriod = 1/2 means lowest period is 12 hours

wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", 
         periodlab = "period (days)")
# - arrows to right - indicate the two series are in phase at the respective
#   period
# - arrows to left - indicate the two series are in anti-phase
# - arrows are only plotted within white contour lines indicating significance
#      at the 10% level (with respect to null hypothesis of white noise)

# individual wavelet power spectra
wt.image(my.wc, my.series = "x")
wt.image(my.wc, my.series = "y")

# time-averaged cross-wavelet power
wc.avg(my.wc, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (days)")




# 3. SEV flux data -------------------------------------------------
# Daily Ameriflux US-Seg data set - SEV LTER desert grassland
# 
# This is an example of running the analysis on Amerifulx data
# from the Sevilleta LTER desert grassland site summarized at the
# daily scale. Variables are NEE and precipitation.
#
# *****  Data should not have any gaps

# US-Seg
# user: change to the file path where the file is located on your computer
seg_all <- read_csv("~/Documents/LTER_Pulse_Time_Series_WG/SEV_daily_flux/US-Seg_daily_aflx.csv")





# retaining only variables of interest, formatting date, and changing name 
# to lower case 'date' so that the variable works with the WaveletComp package.
seg_all <- seg_all %>%
        select(Date, P_int, NEE_int) %>% 
        mutate(date = mdy(Date)) %>% 
        select(-Date) %>% 
        rename(precip = P_int, NEE = NEE_int) %>% 
        select(date, precip, NEE)


summary(seg_all)
str(seg_all)
head(seg_all)


# plot of data
s <- ggplot(seg_all, aes(x = date))

s <- s + geom_line(aes(y = precip, color = "Precipitation (mm)"), alpha = 0.8)
s <- s + geom_line(aes(y = NEE * 10 , color = "NEE (g C/m^2)"), alpha = 0.5)

s <- s + scale_y_continuous(sec.axis = sec_axis(~./10, name = "NEE (g C / m^2)")) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank())

s <- s + scale_color_manual(values = c("green", "blue"))
s <- s + labs(title = "Seg", y = "Precipitation (mm)", x = "Date", color = "Variable")
s 
# don't pay attention to the scale, just trying to see general pattern of data

# CO2 input to the ecosystem is negative NEE



# separating variables to run wavelet analyses separately
seg_precip <- seg_all %>% 
        select(date, precip)

seg_nee <- seg_all %>% 
        select(date, NEE)


# precip wavelet
my.wt <- analyze.wavelet(seg_precip, "precip",make.pval = TRUE, n.sim = 10) # n.simulations will need to be higher for the final figure, but takes very long

wt.image(my.wt, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

wt.image(my.wt, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8),
         plot.contour = FALSE)   # without contour lines


# NEE wavelet
my.wt <- analyze.wavelet(seg_nee, "NEE", make.pval = TRUE, n.sim = 10) # n.simulations will need to be higher for the final figure, but takes very long

wt.image(my.wt, main = "US-Seg Daily NEE",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

wt.image(my.wt, main = "US-Seg Daily NEE",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8),
         plot.contour = FALSE) # without contour lines




# cross-wavelet for precip and NEE
seg_all_2 <- data.frame(seg_all)

my.wc <- analyze.coherency(seg_all_2, c("precip","NEE"), n.sim = 2) # n.simulations will need to be higher for the final figure, but takes very long

wc.image(my.wc, main = "cross-wavelet power spectrum, precip and NEE",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)")

wc.image(my.wc, main = "cross-wavelet power spectrum, precip and NEE",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)",
         plot.contour = FALSE) # without contour lines


# plot of average coherence:
wc.avg(my.wc, which.avg = "wc", 
       siglvl = 0.05, sigcol = 'red', 
       legend.coords = "topleft", 
       periodlab = "period (days)")








# 4. Monthly FCE LTER surface water level and dissolved organic carbon ----------------------
level <- read_csv("~/Documents/LTER_Pulse_Time_Series_WG/FCE_data/fce.srs2.level.csv") %>% 
        rename(level = y)

doc <- read_csv("~/Documents/LTER_Pulse_Time_Series_WG/FCE_data/fce.srs2.doc.csv") %>% 
        rename(doc = y)


fce_all <- level %>% 
        left_join(doc) %>% 
        mutate(date = mdy(date))

s <- ggplot(fce_all, aes(x = date))

s <- s + geom_line(aes(y = doc, color = "DOC"), alpha = 0.8)
s <- s + geom_line(aes(y = level * 10 , color = "Water Level"), alpha = 0.5)

s <- s + scale_y_continuous(sec.axis = sec_axis(~./10, name = "DOC")) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank())

s <- s + scale_color_manual(values = c("green", "blue"))
s <- s + labs(title = "FCE", y = "Water Level", x = "Date", color = "Variable")
s 
# don't pay attention to the scale, just trying to see general pattern of data



# separating variables to run wavelet analyses
fce_level <- fce_all %>% 
        select(date, level)

fce_doc <- fce_all %>% 
        select(date, doc)


# water level wavelet
my.wt <- analyze.wavelet(fce_level, "level",make.pval = TRUE, n.sim = 10) # n.simulations will need to be higher for the final figure, but takes very long

wt.image(my.wt, main = "FCE Water Level",
         periodlab = "period (monthly)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))



# DOC wavelet
my.wt <- analyze.wavelet(fce_doc, "doc",make.pval = TRUE, n.sim = 10) # n.simulations will need to be higher for the final figure, but takes very long

wt.image(my.wt, main = "FCE DOC",
         periodlab = "period (monthly)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))



# cross-wavelet for precip and NEE
fce_all_2 <- data.frame(fce_all)

my.wc <- analyze.coherency(fce_all_2, c("level","doc"), n.sim = 2) # n.simulations will need to be higher for the final figure, but takes very long

wc.image(my.wc, main = "cross-wavelet power spectrum, water level and doc",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (months)")


# plot of average coherence:
wc.avg(my.wc, which.avg = "wc", 
       siglvl = 0.05, sigcol = 'red', 
       legend.coords = "topleft", 
       periodlab = "period (months)")









