# 2022 LTER ASM pulse dynamics workshops WaveletComp R package template
#
# (Author: KM Hall)
#
# The purpose of this R script is to help folks get up and running quickly
# using WaveletComp on your own time series data to conduct basic wavelet analyses. 
#
#
# This template will require modifications and adjustments by the USER to
# conduct wavelet and cross-wavelet analyses, but it is intended to 
# require very minimal tweaking. The places in the code where the USER needs to
# make modifications are commented as "**USER: info about required modification...**".
# The rest of the code should run without modification. 
#
# This accompanies the LTER_Pulses_WaveletComp_R_Package_Tutorial.html document,
# an R Markdown document that walks through some WaveletComp analysis examples.
#
# This template script only conducts very basic analyses, and additional
# code will need to be written to conduct in depth analyses of model results. 
#
# See the WaveletComp guided tour for more information on the package:
# http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf






### The Setup  ---------------------------------------------


# Load R packages
library(WaveletComp)
library(dplyr)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(lubridate)
library(readr)



# Import data

# **USER: You need to change this to the file path for your data on your computer**
#
# This example uses US-Seg-daily_aflx.csv data available on the workshop Google Drive. It
# is daily flux tower data from the AmeriFlux US-Seg site on the Sevilleta NWR in New Mexico
# from 2007-2020.
df <- read_csv("~/Documents/LTER_Pulse_Time_Series_WG/SEV_daily_flux/US-Seg_daily_aflx.csv") %>% 
  mutate(date = mdy(Date)) %>%         # WaveletComp expects date var named "date"
  select(-Date) %>% 
  select(date, everything()) %>% 
  arrange(date)                  # make sure data are properly sorted in chronological order

# NOTE: WaveletComp expects the date variable to by named 'date'. It is best
#       to have dates in YYYY-MM-DD format ("%Y-%m-%d"). If your data is formatted
#       differently, you may need to adjust the read_csv pipeline.
#
# NOTE: WaveletComp expects there to be no missing data. You need to provide
#       the package with a complete time series data set without missing data.
#       This may require a certain amount of data preprocessing to gap fill and
#       make any other necessary adjustments.





# Take a peek at the data
names(df)

str(df)

summary(df)

head(df)


# **USER: Select variables from data for analysis** 
#   w_var is the single variable for wavelet analysis
#   wc_var1 is the first variable for cross-wavelet analysis
#   wc_var2 is the second variable for cross-wavelet analysis
w_var <- "VAR_NAME_HERE"

wc_var1 <- "VAR_NAME_HERE"
wc_var2 <- "VAR_NAME_HERE"


# vars from example data
# w_var <- "TA_avg"
# wc_var1 <- "P_int"
# wc_var2 <- "NEE_int"






### Conduct wavelet analysis --------------------------------------


# This section runs a wavelet analysis on w_var. The USER is encouraged
# to look at the documentation for the various function options 

# reformat data selecting date and w_var
df_univar <- df %>% 
  select(date, w_var) %>% 
  as.data.frame()

# conduct wavelet analysis 
df.w <- analyze.wavelet(df_univar, my.series = 2,
                        loess.span = 0.75,
                        dt = 1, dj = 1/20,
                        make.pval = TRUE, method = "white.noise",
                        n.sim = 25)   # should set n.sim higher - default is 100

# wavelet image
wt.image(df.w,
         col.contour = "white",
         # color.palette = mypalette,
         lwd = 1.5,
         show.date = TRUE,
         main = "Univariate Wavelet",
         legend.params = list(width = 1.8, label.digits = 3, lab = "Power"))


# plot of average power across the entire time series
maximum.level = 1.001*max(df.w$Power.avg)    # See WaveletComp guided tour for explanation
wt.avg(df.w, maximum.level = maximum.level,
       legend.coords = "bottomright",
       main = "Average Power for w_var")

# reconstruct the original time series
reconstruct(df.w,
            show.date = TRUE,
            legend.coords = "bottomright",
            ylim = c(-40, 20),
            main = "Reconstruction of w_var by WaveletComp",
            col = c('tan', 'black'))





### Conduct cross-wavelet analysis --------------------------------------


# This section runs a cross-wavelet analysis on wc_var1 and wc_var2.
# The USER is encouraged to look at the documentation for the various 
# function options.

# reformat data selecting date, wc_var1 and wc_var2
df_bivar <- df %>% 
  select(date, wc_var1, wc_var2) %>% 
  as.data.frame()


# conduct cross-wavelet analysis
df.wc <- analyze.coherency(df_bivar, my.pair = c(2, 3),
                                       loess.span = 0.75,
                                       dt = 1, dj = 1/20,
                                       make.pval = TRUE,
                                       n.sim = 20)

# cross-wavelet image
wc.image(df.wc,
         main = "Cross-wavelet Image",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)",
         col.contour = "white",
         lwd = .8,
         show.date = TRUE,
         col.arrow = "black",
         which.image = "wp"
)






# NOTE: 
#      analyze.wavelet and analyze.coherency functions have several optional
#      arguments that can be adjusted. See package documentation for more
#      details. If your data are not at a daily time scale, you may need to make
#      adjustments to the defaults. 
#
#      Also see package documentation for more information on what data the 
#      analyze.wavelet and analyze.coherency models return. This information
#      is embedded in a list structure that can be manipulated to examine
#      regions of interest as determined from the wavelet or cross-wavelet images.
#      This will require further programming by the USER.



