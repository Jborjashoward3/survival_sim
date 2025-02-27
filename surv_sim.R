
def <- defData(varname = "status", formula = 0.9, dist = "binary")

set.seed(NULL)

sdef <- defSurv(varname = "survTime",scale=1)
sdef2 <- defSurv(varname = "eventtime", scale=0.2)
# sdef <- defSurv(sdef, varname = "censorTime", scale = 80, shape = 1)

sdef


dtSurv <- genData(1000000, def)
dtSurv <- genSurv(dtSurv, sdef)

dtSurv2 <- genData(1000000,def)
dtSurv2 <- genSurv(dtSurv2,sdef2)

dtSurvmain <- left_join(dtSurv,dtSurv2 %>% select(id,eventtime))

df <- dtSurvmain %>% 
  mutate(time=survTime*365,
         vaxtime=runif(1000000,min=20,max=90),
         # vaxtime=rnorm(1000000,mean=55,sd=14),
         # vaxtime=(rt(1000000,10)+8)*7,
         # vaxtime=eventtime*365,
         vaxbinom=rbinom(1000000,1,0.9),
         event=1) %>% 
  mutate(status=ifelse(vaxtime>0 & time>vaxtime &vaxbinom==1,1,0))         
  


appendid <-
  df %>% mutate(id=row_number()) 

beforevaxset <-
  appendid %>% filter(status==1) %>% select(id,tstop=vaxtime) %>% mutate(status=0,tstart=0,event=0) %>% select(id,tstart,tstop,event,status)

aftervaxset <-
  appendid %>% filter(status==1) %>% mutate(tstop=time,event=1) %>% select(id,tstart=vaxtime,tstop,status,event)

nevervaxset <-
  appendid %>% filter(status==0) %>% mutate(tstart=0) %>%  select(id,tstart,status,tstop=time,event) %>% mutate(tstop=ifelse(tstop==0,1,tstop))


tdepset<-
  rbind(beforevaxset,aftervaxset,nevervaxset)


summary(coxph(Surv(tstart, tstop,event) ~  status, tdepset))

yearcensor <-
  tdepset %>% mutate(tstop=ifelse(tstart+tstop<50,50,tstop))

summary(coxph(Surv(tstart, tstop,event) ~  status, yearcensor))



landmark <- 45

setfirstweek <-
  df %>% filter(time>landmark) %>% 
  # filter((status==1 & vaxtime>3 & vaxtime<4) | status==0) %>% 
  
  #if vaxtime>4, label as unvaccinated
  mutate(statuslandmark=ifelse(vaxtime>=landmark+7&status==1,0,status),
         #if vaxtime is >4, need to censor at vaccination, so event label, for those 
         eventlandmark=ifelse(vaxtime>=landmark+7&status==1,0,event),
         timelandmark=ifelse(vaxtime>=landmark+7&status==1,vaxtime,time)-landmark,
         neversubset=ifelse(statuslandmark==1|status==0,1,0))

Yall = Surv(setfirstweek$timelandmark, 
         setfirstweek$eventlandmark == 1)

plot(survfit(Yall ~ setfirstweek$statuslandmark),lty = c("solid", "dashed"), col = c("black", "grey"), xlab = "Survival Time In Days", ylab = "Survival Probabilities",xlim=c(0,500),ylim=c(0.10,1))
legend("topright", c("No vaccination", "Vaccination"), lty = c("solid", "dashed"), col = c("black", "grey"))


Ynever = Surv(setfirstweek[setfirstweek$neversubset==1,]$timelandmark, 
         setfirstweek[setfirstweek$neversubset==1,]$eventlandmark == 1)

plot(survfit(Ynever ~ setfirstweek[setfirstweek$neversubset==1,]$statuslandmark),lty = c("solid", "dashed"), col = c("black", "grey"), xlab = "Survival Time In Days", ylab = "Survival Probabilities",xlim=c(0,500),ylim=c(0.10,1))
legend("topright", c("No vaccination", "Vaccination"), lty = c("solid", "dashed"), col = c("black", "grey"))



########### GROK SUGGESTION TIME SERIES ---------------------------------------


library(tidyverse)
library(ggplot2)
set.seed(NULL)

# Step 1: Simulate base data (in days)
n <- 1000000
time <- rexp(n, 1/(20 * 30))        # Mean 720 days
exposure_time <- runif(n, 24, 120)
exposed <- rbinom(n,1,0.9)

cohort <- data.frame(id = 1:n, time, exposure_time,exposed, time) %>% mutate(status=1)

# Step 2: 7-day windows (0 to 84 days)
days <- seq(0, 365, by = 7)

# Step 3: Longitudinal data - Excluding future exposed (drop only post-exposure)
long_data_excl <- cohort %>%
  group_by(id) %>%
  reframe(
    week_start = days[days <= min(time, 365)],
    exposure_time = exposure_time,
    time = time,
    status = status,
    exposed = exposed
  ) %>%
  mutate(
    at_risk = week_start <= time,  # At risk until event/censoring
    event = ifelse(week_start == floor(time / 7) * 7 & status == 1, 1, 0), 
    # Keep until exposure_time, then exclude (like censoring)
    include = (week_start >= exposure_time+3 & exposed==1) | exposed==0 # include if exposed before landmark or never exposed
  ) %>%
  filter(include)

# Step 4: Longitudinal data - Including future exposed (full cohort)
long_data_incl <- cohort %>%
  group_by(id) %>%
  reframe(
    week_start = days[days <= min(time, 365)],
    exposure_time = exposure_time,
    exposed=exposed,
    time = time,
    status = status
  ) %>%
  mutate(
    at_risk = week_start <= time,
    event = ifelse(week_start == floor(time / 7) * 7 & status == 1, 1, 0)
  )

# Step 5: Summarize overall mortality
summary_excl <- long_data_excl %>%
  group_by(week_start) %>%
  summarise(
    n_at_risk = sum(at_risk),
    events = sum(event),
    prop_mortality = events / n_at_risk,
    .groups = "drop"
  ) %>%
  mutate(analysis = "Excluding Future Exposed")

summary_incl <- long_data_incl %>%
  group_by(week_start) %>%
  summarise(
    n_at_risk = sum(at_risk),
    events = sum(event),
    prop_mortality = events / n_at_risk,
    .groups = "drop"
  ) %>%
  mutate(analysis = "Including Future Exposed")

# Combine
combined_summary <- bind_rows(summary_excl, summary_incl)

# Step 6: ggplot
ggplot(combined_summary, aes(x = week_start, y = prop_mortality, color = analysis)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Overall Mortality by 7-Day Windows: Excluding vs. Including Future Exposed",
    x = "Start of Week (days)",
    y = "Proportion of Events (Mortality Rate)",
    color = "Analysis Type"
  ) +
  scale_x_continuous(breaks = seq(0, 365, by = 28)) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()


