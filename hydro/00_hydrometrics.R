### Hydrologic metrics ###
## Potential covariates for C-Q models derived from Discrete Fast Fourier Transform (DFFT)

## Tamara Harms, 8/24

library(here)
library(tidyverse)
library(googledrive)
library(discharge)
library(data.table)
library(lubridate)
library(zoo)

######################
### Discharge data ###
######################
qURL <- "https://drive.google.com/drive/folders/1Vtyv6AC-QD8Ss-KqbguSbiiDedgBYm14"

Q_new <- drive_get(as_id(qURL))

Q_glist <- drive_ls(Q_new, pattern = "discharge_daily.rds")

#I know this is hinky with the switching of working directory. The drive_download command does not currently allow downloading into subdirectories.
setwd(here("data"))
walk(Q_glist$id, ~ drive_download(as_id(.x), overwrite = TRUE))

setwd(here())

dat <- readRDS(here("data", "discharge_daily.rds"))

## Select test sites
t.sites <- read.csv(here("hydro", "fires_test_sites_ed.csv"))
t.sites <- t.sites %>% separate(usgs_site, c("USGS", "usgs_site"), "-")
t.sites <- t.sites %>% separate(ignition_date, c("igDate", "junk"), "T") %>%
                       mutate(igDate = as.Date(igDate, format = "%Y-%m-%d"))

t.dat <- dat %>% filter(site_no %in% t.sites$usgs_site)
t.dat <- t.dat %>% separate(usgs_site, c("trash", "usgs_site"), sep = "-")

##############################
### Prep data for analysis ###
##############################
## Pare data to 25 y pre-fire
# Join ignition dates
t.dat <- full_join(t.dat, t.sites, by = "usgs_site")

# Filter to 25 y before fire
t.dat <- t.dat %>% mutate(t.diff = Date %--% igDate/years(1)) %>%
                   filter(t.diff <= 25)

# Format for discharge package: Date (YYYY-MM-DD) in column1, CFS (val) in column 2
t.dat <- t.dat %>% select(usgs_site, Date, Flow) %>%
  setNames(c("site", "date", "cfs"))

tdat.lst <- split(t.dat, t.dat$site)

tdat.lst <- lapply(tdat.lst, select, -"site")

# remove too short records for 09367540 & 09367580
tdat.lst[["09367540"]] <- NULL
tdat.lst[["09367580"]] <- NULL

#####################
### DFFT analysis ###
#####################
flow.list <- lapply(tdat.lst, asStreamflow)
seas.list <- lapply(flow.list, fourierAnalysis)

## DFFT fit plots
lapply(names(seas.list), function(x) {
  pdf(file = paste0("hydro/plots/", "DFFT-", x, ".pdf"))
  print(plot(seas.list[[x]]))
  title(paste0("\n", x))
  dev.off()
})

#######################
### Summary metrics ###
#######################
## Metrics: FPExt, HSAF, HSAM, LSAF, LSAM, NAA, rms.noise, rms.signal, snr, Zdays, Seasonal, HFsigma, LFsigma
# All metrics that are not originally composite (such as HFsigma, SNR, etc.) are averaged over the full Q record
### ***Currently calculated on calendar years- should convert to water year basis?*** Or running values?
DFFTavg <- function(X.flow, X.seas, SiteName){
  X.bl <- prepareBaseline(X.flow)
  X.FPExt <- getFPExt(X.bl$resid.sig, X.flow$data$year)
  X.HSAF <- getHSAF(X.bl$resid.sig, X.flow$data$year)
  X.HSAM <- getHSAM(X.bl$resid.sig, X.flow$data$year)%>%
    select(year, HSAM)
  X.LSAF <- getLSAF(X.bl$resid.sig, X.flow$data$year)
  X.LSAM <- getLSAM(X.bl$resid.sig, X.flow$data$year)%>%
    select(year, LSAM)
  X.NAA <- getNAA(X.bl$resid.sig, X.flow$data$year)
  X.HFsigma <- sigmaHighFlows(X.flow) #Use X.HFsigma$sigma.hfb
  X.LFsigma <- sigmaLowFlows(X.flow) #Use X.LFsigma$sigma.hfb
  X.Zdays <- X.seas$signal%>% #Average number of zero flow days per year
    group_by(year)%>%
    summarize(Zdays=sum(discharge <0.1))
  #Find average of each DFFT metric
  X.fftavg <- bind_cols(Site=SiteName, HSAM=mean(X.HSAM$HSAM, na.rm=T),LSAM=mean(X.LSAM$LSAM, na.rm=T),
                       HSAF=mean(X.HSAF$HSAF, na.rm=T), LSAF=mean(X.LSAF$LSAF, na.rm=T), FPExt=mean(X.FPExt$FPExt),
                       NAA=mean(X.NAA$NAA), HFsigma=mean(X.HFsigma$sigma.hfb), LFsigma=mean(X.LFsigma$sigma.lfb))
  X.fftavg <- bind_cols(X.fftavg, rms.signal=X.seas$rms$rms.signal,
                      rms.noise=X.seas$rms$rms.noise, snr=X.seas$rms$snr,
                      Zdays=mean(X.Zdays$Zdays, na.rm=T), Seasonal=X.seas$seasonal)
}

met.out <- mapply(x = flow.list, y = seas.list, z = names(seas.list), FUN = function(x, y, z) DFFTavg(x, y, z))
met.out.df <- met.out[-1,] %>%
              as_tibble(., rownames = "metric") %>%
              unnest(-metric) %>%
              pivot_longer(!metric, names_to = "site", values_to = "value")

write.csv(met.out, here("hydro", "DFFT_avg_25years.csv"), row.names = FALSE)

## Plot
metrics.pl <- met.out.df %>% ggplot(aes(x = site, y = value)) +
                                geom_point() +
                                facet_wrap(~metric, scales = "free_y") +
                                theme_bw() +
                                theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank())

ggsave(metrics.pl, path = here("hydro", "plots"), file = "DFFT_metrics_25y.pdf", width = 10, height = 8, units = "in")

#################################################################
### Summarize metrics over 4 years pre- and 4 years post-fire ###
#################################################################
### ****Goals here- ave metrics pre- and post-fire, then sum of + anomalies post-fire-- filter flow and seas to 4 y pre & 4 y post fire first?
# Metrics by year
####****Summarizing by calendar year... Switch to years in fire timeline or water years?***####
DFFTyr <- function(X.flow, X.seas, SiteName){
  X.bl <- prepareBaseline(X.flow)
  X.FPExt <- getFPExt(X.bl$resid.sig, X.flow$data$year)
  X.HSAF <- getHSAF(X.bl$resid.sig, X.flow$data$year)
  X.HSAM <- getHSAM(X.bl$resid.sig, X.flow$data$year)%>%
    select(year, HSAM)
  X.LSAF <- getLSAF(X.bl$resid.sig, X.flow$data$year)
  X.LSAM <- getLSAM(X.bl$resid.sig, X.flow$data$year)%>%
    select(year, LSAM)
  X.NAA <- getNAA(X.bl$resid.sig, X.flow$data$year)
  X.HFsigma <- sigmaHighFlows(X.flow) #Use X.HFsigma$sigma.hfb
  X.LFsigma <- sigmaLowFlows(X.flow) #Use X.LFsigma$sigma.hfb
  X.Zdays <- X.seas$signal%>% #Average number of zero flow days per year
    group_by(year)%>%
    summarize(Zdays = sum(discharge <0.1))
  X.fftavg <- bind_cols(Site = SiteName, HSAM = X.HSAM$HSAM, LSAM = X.LSAM$LSAM,
                        HSAF = X.HSAF$HSAF, LSAF = X.LSAF$LSAF, FPExt = X.FPExt$FPExt,
                        NAA = X.NAA$NAA, HFsigma = X.HFsigma$sigma.hfb, LFsigma = X.LFsigma$sigma.lfb,
                        X.Zdays)
}

met.out.y <- mapply(x = flow.list, y = seas.list, z = names(seas.list), FUN = function(x, y, z) DFFTyr(x, y, z))

#met.out.y.df <- met.out.y %>%
#  as_tibble(., rownames = "metric") %>%
#  unnest(-metric) 
#%>%
#  pivot_longer(!metric, names_to = "site", values_to = "value")

## Try single site...need to get year as a column
## This works
ss.flow <- flow.list[["07103700"]]
ss.seas <- seas.list[["07103700"]]

ss.metrics <- DFFTyr(ss.flow, ss.seas, "07103700")

## filter to pre-fire & filter to post-fire
# Note change years(x) if revising pre-/post-fire window
t.sites.pp <- t.sites %>% mutate(before.fire = igDate - years(4)) %>%
                          mutate(after.fire = igDate + years(4)) %>%
                          mutate(ybeffire = format(as.Date(before.fire, format = "%Y-%m-%d"), "%Y")) %>%
                          mutate(yaftfire = format(as.Date(after.fire, format = "%Y-%m-%d"), "%Y"))

names(ss.metrics)[names(ss.metrics) == 'Site'] <- 'usgs_site'
ss.metrics.f <- left_join(ss.metrics, t.sites.pp, by = "usgs_site") 

###########################
### High-flow anomalies ###
###########################
# sum of pre-fire and post-fire high-flow anomalies

## Extract positive anomalies
panom <- lapply(seas.list, function(x) subset(x$signal, x$signal$resid.sig >= 0) )

# Join site IDs to anomalies
panom <- Map(cbind, panom, usgs_site = names(panom))
panom <- lapply(panom, function(x) left_join(x, t.sites.pp, by = "usgs_site"))

# sum pre- and post-fire positive anomalies
pre.sum <- panom %>% map(~ .x %>%
                           filter(date >= before.fire & date <= igDate) %>%
                           summarize(bpos = sum(resid.sig)))

post.sum <- panom %>% map(~ .x %>%
                           filter(date >= igDate & date <= after.fire) %>%
                           summarize(apos = sum(resid.sig)))

dates <- panom %>% map(~ .x %>%
                         select(c("usgs_site", "event_id", "igDate", "before.fire", "after.fire")) %>%
                         unique())

# combine
pre.sum.df <- bind_rows(pre.sum, .id = "usgs_site")
post.sum.df <- bind_rows(post.sum, .id = "usgs_site")
dates.df <- bind_rows(dates, .id = "usgs_site")

pos.anom <- left_join(dates.df, pre.sum.df)
pos.anom <- left_join(pos.anom, post.sum.df)

###***other possibilities:
###*days since fire of largest positive anomaly OR some fraction of cdf of net anomalies?
###*number days with positive anomalies post-fire OR number days with + anomalies above some threshold (like 2*SD mean of + anomalies)
###*Average time between positive anomalies above some threshold

