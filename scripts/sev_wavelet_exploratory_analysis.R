



library(WaveletComp)
library(dplyr)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(lubridate)
library(readr)








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


s <- ggplot(seg_all, aes(x = date))

s <- s + geom_line(aes(y = precip, color = "Precipitation (mm)"), alpha = 0.8)
s <- s + geom_line(aes(y = NEE * 10 , color = "NEE (g C/m^2)"), alpha = 0.5)

s <- s + scale_y_continuous(sec.axis = sec_axis(~./10, name = "NEE (g C / m^2)")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

s <- s + scale_color_manual(values = c("green", "blue"))
s <- s + labs(title = "Seg", y = "Precipitation (mm)", x = "Date", color = "Variable")
s 








# precip ------------

library(timetk)

seg_all %>% 
  plot_time_series(date, precip,
                   .color_var = year(date),
                   
                   .interactive = TRUE,
                   .plotly_slider = TRUE,
                   
                   .title = "US-Seg Daily Precipitation",
                   .x_lab = "Time",
                   .y_lab = "Precipitation (mm)")


library(TSstudio)

precip <- seg_all %>% 
  select(date, precip) 
  
ts_heatmap(precip, wday = FALSE, title = "US-Seg Precipitation")


library(feasts)

precip <- as_tsibble(precip)

precip %>% 
  model(STL(precip ~ season(period = "1 year"))) %>% 
  components() %>% 
  autoplot()


precip_STL.m <- precip %>% 
  model(STL(precip ~ season(period = "1 year")))


# precipitation after removing season and trend
#
# going to see if this provides anything interesting with wavelet analysis
precip_remainder <- precip_STL.m %>% components() %>% select(remainder)




# precip wavelets --------------------------

# only the remainder of the time series -
#   was curious about whether this would show any interesting patterns
#   after removing seasonality and trend from the data - based on STL decomposition
precip_remainder.wt <- analyze.wavelet(precip_remainder, "remainder", make.pval = TRUE, n.sim = 100) 

wt.image(precip_remainder.wt, main = "US-Seg Daily Precipitation Remainder \n(Season and Trend REMOVED)",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip_remainder.wt$Power.avg)
wt.avg(precip_remainder.wt, maximum.level = maximum.level)






# precip - base model - 0 loess
precip.wt.base0 <- analyze.wavelet(precip, "precip", 
                             make.pval = TRUE, 
                             n.sim = 100,
                             loess.span = 0) 

wt.image(precip.wt.base0, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base0$Power.avg)
wt.avg(precip.wt.base0, maximum.level = maximum.level)






# precip - base model with loess time series smoothing -
#
#    loess.span = 0.75 is the default
precip.wt.base <- analyze.wavelet(precip, "precip", 
                                  make.pval = TRUE, 
                                  n.sim = 100,
                                  loess.span = 0.75)
                                  

wt.image(precip.wt.base, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base$Power.avg)
wt.avg(precip.wt.base, maximum.level = maximum.level)


# looking at the low power region from late 2010 into 2011 -
#    for several months, there was very low precip
precip %>% 
  filter(year(date) %in% c(2009, 2010, 2011, 2012)) %>% 
  ggplot(aes(x = date, y = precip)) +
  geom_point()

precip %>% 
  as_tibble() %>% 
  filter(year(date) %in% c(2009, 2010, 2011, 2012)) %>%
  mutate(year = year(date)) %>% 
  group_by(year) %>% 
  summarize(sum(precip))

precip %>% 
  as_tibble() %>% 
  filter(year(date) %in% c(2009, 2010, 2011, 2012)) %>%
  mutate(year_month = ym(paste(year(date), " ", month(date)))) %>% 
  group_by(year_month) %>% 
  summarize(precip = sum(precip)) %>% 
  ggplot(aes(x = year_month, y = precip)) +
  geom_line()


# looking at a high power region in 2013
precip %>% 
  filter(year(date) %in% c(2012, 2013, 2014)) %>% 
  ggplot(aes(x = date, y = precip)) +
  geom_point()

precip %>% 
  as_tibble() %>% 
  filter(year(date) %in% c(2012, 2013, 2014)) %>% 
  mutate(year = year(date)) %>% 
  group_by(year) %>% 
  summarize(sum(precip))

precip %>% 
  as_tibble() %>% 
  filter(year(date) %in% c(2012, 2013, 2014)) %>% 
  mutate(year_month = ym(paste(year(date), " ", month(date)))) %>% 
  group_by(year_month) %>% 
  summarize(precip = sum(precip)) %>% 
  ggplot(aes(x = year_month, y = precip)) +
  geom_line()





# method 'shuffle' - 
precip.wt.base.shuffle <- analyze.wavelet(precip, "precip", 
                                        make.pval = TRUE, 
                                        n.sim = 100,
                                        loess.span = 0.75,
                                        method = "shuffle")


wt.image(precip.wt.base.shuffle, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base.shuffle$Power.avg)
wt.avg(precip.wt.base.shuffle, maximum.level = maximum.level)






# method 'Fourier.rand' - 
precip.wt.base.Fourier.rand <- analyze.wavelet(precip, "precip", 
                                          make.pval = TRUE, 
                                          n.sim = 100,
                                          loess.span = 0.75,
                                          method = "Fourier.rand")


wt.image(precip.wt.base.Fourier.rand, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base.Fourier.rand$Power.avg)
wt.avg(precip.wt.base.Fourier.rand, maximum.level = maximum.level)






# method 'AR' - 
precip.wt.base.AR <- analyze.wavelet(precip, "precip", 
                                          make.pval = TRUE, 
                                          n.sim = 100,
                                          loess.span = 0.75,
                                          method = "AR")


wt.image(precip.wt.base.AR, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base.AR$Power.avg)
wt.avg(precip.wt.base.AR, maximum.level = maximum.level)




# method 'AR' 2- 
precip.wt.base.AR2 <- analyze.wavelet(precip, "precip", 
                                     make.pval = TRUE, 
                                     n.sim = 100,
                                     loess.span = 0.75,
                                     method = "AR",
                                     params = list(p = 2))


wt.image(precip.wt.base.AR2, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base.AR2$Power.avg)
wt.avg(precip.wt.base.AR2, maximum.level = maximum.level)








# method 'ARIMA' - 
precip.wt.base.ARIMA <- analyze.wavelet(precip, "precip", 
                                          make.pval = TRUE, 
                                          n.sim = 100,
                                          loess.span = 0.75,
                                          method = "ARIMA")


wt.image(precip.wt.base.ARIMA, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base.ARIMA$Power.avg)
wt.avg(precip.wt.base.ARIMA, maximum.level = maximum.level)








# different octave - 
precip.wt.base.dj8 <- analyze.wavelet(precip, "precip", 
                                        make.pval = TRUE, 
                                        n.sim = 100,
                                        dj = 1/100)


wt.image(precip.wt.base.dj8, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base.dj8$Power.avg)
wt.avg(precip.wt.base.dj8, maximum.level = maximum.level)






# different time resolution - 
precip.wt.base.dt365 <- analyze.wavelet(precip, "precip", 
                                      make.pval = TRUE, 
                                      n.sim = 100,
                                      dt = 1/365)


wt.image(precip.wt.base.dt365, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))

# look at the average power of the series
maximum.level = 1.001*max(precip.wt.base.dt365$Power.avg)
wt.avg(precip.wt.base.dt365, maximum.level = maximum.level)















# looking at base model 

str(precip.wt.base)
class(precip.wt.base)


head(precip.wt.base$series)

plot(seq(from = 1, to = 195, by = 1), precip.wt.base$Power.avg)

precip.wt.base$Ampl
power <- as.matrix(precip.wt.base$Power)

precip.wt.base$Period

# one row of the power column
power_row8 <- power[8, ]


seg_all_date <- seg_all$date

power_r8_df <- tibble(date = seg_all_date,
                      power8 = power_row8)

power_r8_df %>% 
  ggplot(aes(date, power8)) +
  geom_line()



# another row of the power column - power 362 (close to 1 year=365)
power_row151 <- power[151, ]


seg_all_date <- seg_all$date

power_r151_df <- tibble(date = seg_all_date,
                      power = power_row151)

power_r151_df %>% 
  ggplot(aes(date, power)) +
  geom_line()


# NOTE: Period values are slices of the periods to be looked at
#    and correspond to the rows of the Power matrix. Columns
#    of the Power matric are the data records (days, in this case)




# another row of the power column - period 16
power_row61 <- power[61, ]


seg_all_date <- seg_all$date

power_r61_df <- tibble(date = seg_all_date,
                        power = power_row61)

power_r61_df %>% 
  ggplot(aes(date, power)) +
  geom_line()


wt.image(precip.wt.base, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))




# another row of the power column - period 64
power_row101 <- power[101, ]


seg_all_date <- seg_all$date

power_r101_df <- tibble(date = seg_all_date,
                       power = power_row101)

power_r101_df %>% 
  ggplot(aes(date, power)) +
  geom_line()


wt.image(precip.wt.base, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))



# another row of the power column - period 7 - about a week
power_row37 <- power[37, ]


seg_all_date <- seg_all$date

power_r37_df <- tibble(date = seg_all_date,
                        power = power_row37)

power_r37_df %>% 
  ggplot(aes(date, power)) +
  geom_line()


wt.image(precip.wt.base, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))



# another row of the power column - period 30 - about a month
power_row79 <- power[79, ]


seg_all_date <- seg_all$date

power_r79_df <- tibble(date = seg_all_date,
                        power = power_row79)

power_r79_df %>% 
  ggplot(aes(date, power)) +
  geom_line()


wt.image(precip.wt.base, main = "US-Seg Daily Precipitation",
         periodlab = "period (daily)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))










# looking at data at monthly scale

seg_m <- seg_all %>% 
  mutate(year_month = ym(paste(year(date), " ", month(date)))) %>% 
  group_by(year_month) %>% 
  summarize(precip = sum(precip),
            NEE = sum(NEE)) %>% 
  rename(date = year_month)

precip_m.wt.base <- analyze.wavelet(seg_m, "precip", 
                                    make.pval = TRUE, 
                                    n.sim = 100)

wt.image(precip_m.wt.base, main = "US-Seg Monthly Precipitation",
         periodlab = "period (monthly)",
         label.time.axis = T, show.date = T, date.format = "%Y-%m-%d",
         color.key = "quantile",legend.params = list(label.digits = 3, lab = "wavelet power levels", mar = 8))








# looking at details of WaveletComp functions
# precip.wt.base

wt.phase.image(precip.wt.base)


wt.image(precip.wt.base, 
         siglvl = .05,
         # col.contour = "darkgrey",
         lwd = 1,
         n.levels = 30,
         show.date = TRUE,
         main = "US-Seg Flux Tower Precipitation (mm)")
# this is kind of helpful in that it shows some of the changes
# in periodicity power over time a little more coarsely -
#    default is 100

reconstruct(precip.wt.base)
















# sev met 40 daily ------------------------------

met <- read_csv("~/Documents/SEV/Projects/OUP_Climate_Chapter/data/processed_data/met_daily_gap_filled.csv") %>% 
  mutate(sta = as.factor(sta)) %>% 
  filter(sta == "40")

met.wt <- analyze.wavelet(met, "airt", 
                          make.pval = TRUE, 
                          n.sim = 20)
# it smooths the time series - might be interesting to try 
# feeding it an STL decomposition's residuals.



library(viridis)
# mypalette <-  "mako(100)" #"turbo(100)"      #"magma(100)"       #"cividis(100)"         # "viridis(100)",                 # "topo.colors(100)",
mypalette <-  "terrain.colors(100)" 



wt.image(met.wt,
         col.contour = "white",
         color.palette = mypalette,
         lwd = 1,
         show.date = TRUE,
         main = "SEV Met40 Daily Mean Air Temperature (C)")

# 120 & 220

met.wt_power <- as.matrix(met.wt$Power)


met.wt$Period


# 120
power_119 <- tibble(date = met$date,
                    power = met.wt_power[119, ])
                     
power_119 %>% 
  ggplot(aes(date, power)) +
  geom_line()

# 220
power_136 <- tibble(date = met$date,
                    power = met.wt_power[136, ])

power_136 %>% 
  ggplot(aes(date, power)) +
  geom_line()




tibble(date = met$date,
       power = met.wt_power[134, ]) %>%
  ggplot(aes(date, power)) +
  geom_line()


tibble(date = met$date,
       power = met.wt_power[151, ]) %>%
  ggplot(aes(date, power)) +
  geom_line()

met %>% 
  ggplot(aes(x = date, y = airt)) +
  geom_line()

met %>% 
  ggplot(aes(x = date, y = maxair)) +
  geom_line()
# skeptical of 40's data prior to 1991 - is it the data or did I fuck up? 
#    look back at pre-1991 data from EDI - how bad are these variables, and did the majority
#    get filled with hourly means in gap filling? That could account for the lack of variance 
#    of the time series in that time range. It seems the noise is really reduced in that time range.




# solar 
met.sol.wt <- analyze.wavelet(met, "sol_total", 
                          make.pval = TRUE, 
                          n.sim = 20)

# mypalette <-  "cividis(100)"

wt.image(met.sol.wt,
         col.contour = "white",
         color.palette = mypalette,
         lwd = 1,
         show.date = TRUE,
         main = "SEV Met40 Daily Solar Radiation")

# Data for 40 prior to 1991 is problematic, as is using gap filling. There is a big chunk of 
# data missing just before 1991, min/maxair missing - thereby leading to problems with calculating
# daily mean airt. Solar also seems to be lacking enough variance pre-1991.

# 40 around period 220 is oddly low at fairly regular intervals. How do other stations look?



# need to look at raw met 40 data

met40_raw <- read_csv("~/Documents/SEV/Projects/OUP_Climate_Chapter/data/raw_data/sev_hrly_met_40_42_49_50.csv") %>% 
  filter(StationID == 40) %>% 
  mutate(sta = as.factor(StationID)) %>% 
  select(-StationID) %>% 
  select(sta, everything()) %>% 
  rename(dt = Date_Time, date = Date, airt = Temp_C, minair = Min_Temp_C, maxair = Max_Temp_C,
         sol = Solar_Radiation)

met40_raw %>% 
  ggplot(aes(dt, airt)) +
  geom_line()


met40_raw %>% 
  ggplot(aes(dt, minair)) +
  geom_line()

met40_raw %>% 
  ggplot(aes(dt, maxair)) +
  geom_line()

met40_raw %>% 
  ggplot(aes(dt, sol)) +
  geom_line()






# met 50 ----------  

met50 <- read_csv("~/Documents/SEV/Projects/OUP_Climate_Chapter/data/processed_data/met_daily_gap_filled.csv") %>% 
  mutate(sta = as.factor(sta)) %>% 
  filter(sta == "50")

met50.wt <- analyze.wavelet(met50, "airt", 
                          make.pval = TRUE, 
                          n.sim = 20)
# it smooths the time series - might be interesting to try 
# feeding it an STL decomposition's residuals.



# mypalette <-  "mako(100)" 

wt.image(met50.wt,
         col.contour = "white",
         color.palette = mypalette,
         lwd = 1,
         show.date = TRUE,
         main = "SEV met50 Daily Mean Air Temperature (C)",
         
         legend.params = list(width = 1.8, label.digits = 4))


reconstruct(met50.wt,
            show.date = TRUE,
            legend.coords = "bottomright",
            ylim = c(-40, 20),
            main = "Reconstruction of Met 40 \nDaily Mean Airt Temperature",
            col = c('tan', 'black'))

# 120 & 220

met50.wt$Period

met50.wt_power <- as.matrix(met50.wt$Power)

# 120
power_119 <- tibble(date = met50$date,
                    power = met50.wt_power[119, ])

power_119 %>% 
  ggplot(aes(date, power)) +
  geom_line()

# 220
power_136 <- tibble(date = met50$date,
                    power = met50.wt_power[136, ])

power_136 %>% 
  ggplot(aes(date, power)) +
  geom_line()




tibble(date = met50$date,
       power = met50.wt_power[134, ]) %>%
  ggplot(aes(date, power)) +
  geom_line()


tibble(date = met50$date,
       power = met50.wt_power[151, ]) %>%
  ggplot(aes(date, power)) +
  geom_line()

met50 %>% 
  ggplot(aes(x = date, y = airt)) +
  geom_line()

met50 %>% 
  ggplot(aes(x = date, y = maxair)) +
  geom_line()




# looking at power levels in periods of ~220-256
library(purrr)

met50_200s <- tibble(date = rep(met50$date),
                     transpose(met50.wt_power[131:141,]))

period_range <- rep(met50.wt$Period[131:141], each = 10*nrow(met50))



# met50_200s <- tibble(date = met50$date,
#                      period_range = met50.wt$Period[131],
#                      power = met50.wt_power[131])

met50_200s_power <- as.matrix(met50.wt_power[131:141,])

dimnames(met50_200s_power) <- list(paste0("period_", round(met50.wt$Period[131:141], 2)))

met50_200s_power_t <- as_tibble(t(met50_200s_power))

met50_200s_df <- tibble(date = met50$date, met50_200s_power_t)


met50_200s_df_long <- met50_200s_df %>% 
  pivot_longer(contains("period"))

met50_200s_df_long %>% 
  ggplot(aes(x = date, y = value, color = name)) +
  geom_line() +
  facet_wrap(~ name)


mean(met50.wt$Power)
max(met50.wt$Power)
quantile(met50.wt$Power, .99)






# 16-64
met50.wt$Period

met50_16_64_power <- as.matrix(met50.wt_power[61:101,])


met50_16_64 <- tibble(date = rep(met50$date),
                     transpose(met50.wt_power[61:101,]))

period_range <- rep(met50.wt$Period[61:101], each = 40*nrow(met50))


dimnames(met50_16_64_power) <- list(paste0("period_", round(met50.wt$Period[61:101], 2)))

met50_16_64_power_t <- as_tibble(t(met50_16_64_power))

met50_16_64_df <- tibble(date = met50$date, met50_16_64_power_t)


met50_16_64_df_long <- met50_16_64_df %>% 
  pivot_longer(contains("period"))

met50_16_64_df_long %>% 
  ggplot(aes(x = date, y = value, color = name)) +
  geom_line() +
  facet_wrap(~ name) +
  theme(legend.position="none")
  

met50 %>% 
  filter(year(date) == 2010:2015) %>% 
  ggplot(aes(x = date, y = airt)) +
  geom_line()

wt.image(met50.wt,
         col.contour = "white",
         color.palette = mypalette,
         lwd = 1,
         show.date = TRUE,
         main = "SEV met50 Daily Mean Air Temperature (C)",
         
         legend.params = list(width = 1.8, label.digits = 4))








# cross-wavelet - maxair and ppt

met50_maxair_ppt <- met50 %>% 
  select(date, maxair, ppt) %>% 
  as.data.frame()

met50_maxair_pppt.wc <- analyze.coherency(met50_maxair_ppt, c("maxair", "ppt"), n.sim = 20)






# cross-wavelet power - which.image = "wp"
wc.image(met50_maxair_pppt.wc, main = "Cross-wavelet power spectrum, maxair and ppt",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)",
         col.contour = "white",
         color.palette = mypalette,
         lwd = .8,
         show.date = TRUE,
         col.arrow = "red",
         which.image = "wp"
)

# cross-wavelet coherence - which.image = "wc"
wc.image(met50_maxair_pppt.wc, main = "Wavelet Coherence, maxair and ppt",
         # legend.params = list(lab = "cross-wavelet coherence"),
         periodlab = "period (days)",
         col.contour = "white",
         color.palette = mypalette,
         lwd = .8,
         show.date = TRUE,
         which.image = "wc"
)

# plot of average coherence:
wc.avg(met50_maxair_pppt.wc, which.avg = "wc", 
       siglvl = 0.05, sigcol = 'red', 
       legend.coords = "topleft", 
       periodlab = "period (days)")




# subsetting to 2010-2015 - curious how it looks with shorter modeling runs

met50_maxair_ppt_10_14 <- met50 %>% 
  filter(year(date) %in% 2010:2014) %>% 
  select(date, maxair, ppt) %>% 
  as.data.frame()

met50_maxair_ppt_10_14.wc <- analyze.coherency(met50_maxair_ppt_10_14, c("maxair", "ppt"), n.sim = 20)


met50_ppt_maxair_10_14.wc <- analyze.coherency(met50_maxair_ppt_10_14, c("ppt", "maxair"), n.sim = 20)


mypalette <-  "terrain.colors(100)" 




# cross-wavelet power - which.image = "wp"
wc.image(met50_maxair_ppt_10_14.wc, main = "Cross-wavelet power spectrum, maxair and ppt\n2010-2015",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)",
         col.contour = "white",
         color.palette = mypalette,
         lwd = .8,
         show.date = TRUE,
         col.arrow = "red",
         which.image = "wp"
)

# cross-wavelet coherence - which.image = "wc"
wc.image(met50_maxair_ppt_10_14.wc, main = "Wavelet Coherence, maxair and ppt \n2010-2015",
         # legend.params = list(lab = "cross-wavelet coherence"),
         periodlab = "period (days)",
         col.contour = "white",
         color.palette = mypalette,
         lwd = .8,
         show.date = TRUE,
         which.image = "wc"
)

# plot of average coherence:
wc.avg(met50_maxair_ppt_10_14.wc, which.avg = "wc", 
       siglvl = 0.05, sigcol = 'red', 
       legend.coords = "topleft", 
       periodlab = "period (days)")



wc.image(met50_ppt_maxair_10_14.wc, main = "Cross-wavelet power spectrum, ppt and maxair\n2010-2014",
         legend.params = list(lab = "cross-wavelet power levels"),
         periodlab = "period (days)",
         col.contour = "white",
         color.palette = mypalette,
         lwd = .8,
         show.date = TRUE,
         col.arrow = "red",
         which.image = "wp"
)


