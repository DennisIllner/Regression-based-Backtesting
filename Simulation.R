library(rugarch)
library(doParallel)
library(rugarch)
library(dplyr)
library(sandwich)
library(car)
library(doParallel)
library(foreach)
library(tidyverse)

# Defining function for different MZ tests
sim_func_grach <- function(
    MCreps=100, core.max = 7, k=.95,tt){
  
  MZ.test.mean <- function(Proxy, forecast){
    
    ######## Classic MZ test
     MZ.fit <- lm(Proxy~forecast)
        pvalN <- car::linearHypothesis(MZ.fit, c("(Intercept)=0", "forecast=1"), vcov. = sandwich::NeweyWest(MZ.fit))$`Pr(>F)`[2]
        pvalH <- car::linearHypothesis(MZ.fit, c("(Intercept)=0", "forecast=1"), vcov. = sandwich::vcovHAC(MZ.fit))$`Pr(>F)`[2]
    
    ######## log-Tranformations
        alphadj<--log(2)-(-digamma(1))
        logMZ.fit<- lm(log(Proxy) ~ log(forecast))
        logpvalN <- car::linearHypothesis(logMZ.fit, c(paste0("(Intercept)=",alphadj), "log(forecast)=1"), vcov. = sandwich::NeweyWest(logMZ.fit))$`Pr(>F)`[2]
        logpvalH <- car::linearHypothesis(logMZ.fit, c(paste0("(Intercept)=",alphadj), "log(forecast)=1"), vcov. = sandwich::vcovHAC(logMZ.fit))$`Pr(>F)`[2]
    # without correction
        logpvalNnc <- car::linearHypothesis(logMZ.fit, c("(Intercept)=0", "log(forecast)=1"), vcov. = sandwich::NeweyWest(logMZ.fit))$`Pr(>F)`[2]
        logpvalHnc <- car::linearHypothesis(logMZ.fit, c("(Intercept)=0", "log(forecast)=1"), vcov. = sandwich::vcovHAC(logMZ.fit))$`Pr(>F)`[2]
   
    ######## Absolute Returns on  forecast^0.5
        betaadj<-(2/pi)^0.5
        absMZ.fit <-lm(abs(Proxy)~sqrt(forecast))
        abspvalN <- car::linearHypothesis(absMZ.fit, c("(Intercept)=0", paste0("sqrt(forecast)=",betaadj)), vcov. = sandwich::NeweyWest(absMZ.fit))$`Pr(>F)`[2]
        abspvalH <- car::linearHypothesis(absMZ.fit, c("(Intercept)=0", paste0("sqrt(forecast)=",betaadj)), vcov. = sandwich::vcovHAC(absMZ.fit))$`Pr(>F)`[2]
    # withoutcorrection
        abspvalNnc <- car::linearHypothesis(absMZ.fit, c("(Intercept)=0", "sqrt(forecast)=1"), vcov. = sandwich::NeweyWest(absMZ.fit))$`Pr(>F)`[2]
        abspvalHnc <- car::linearHypothesis(absMZ.fit, c("(Intercept)=0", "sqrt(forecast)=1"), vcov. = sandwich::vcovHAC(absMZ.fit))$`Pr(>F)`[2]
    
    ######## F-Tests:
        FMZ.fit<-lm((Proxy-forecast)~forecast)
    # MZ-F
        globalFpvalN<-car::linearHypothesis(FMZ.fit, c("forecast=0"),.vcov = sandwich::NeweyWest(FMZ.fit))$`Pr(>F)`[2]
        globalFpvalH<-car::linearHypothesis(FMZ.fit, c("forecast=0"),.vcov = sandwich::vcovHAC(FMZ.fit))$`Pr(>F)`[2]
    # MZ-F* (Wald with c(0,0)):
        waldFpvalN<-car::linearHypothesis(FMZ.fit, c("(Intercept)=0","forecast=0"),.vcov = sandwich::NeweyWest(FMZ.fit))$`Pr(>F)`[2]
        waldFpvalH<-car::linearHypothesis(FMZ.fit, c("(Intercept)=0","forecast=0"),.vcov = sandwich::vcovHAC(FMZ.fit))$`Pr(>F)`[2]
    
    # p-values
    return(list(pValN       =pvalN,
                pValH       =pvalH,
             logpValN       =logpvalN,
             logpValH       =logpvalH,
           logpValNnc       =logpvalNnc,
           logpValHnc       =logpvalHnc,
             abspValN       =abspvalN,
             abspValH       =abspvalH,
           abspValNnc       =abspvalNnc,
           abspValHnc       =abspvalHnc,
         globalFpValN       =globalFpvalN,
         globalFpValH       =globalFpvalH,
           waldFpValN       =waldFpvalN,
           waldFpValH       =waldFpvalH))
    
  }
  
  # Cluster configuration (parallelization)
  
  cl <- makeCluster(min(parallel::detectCores()-1, core.max) )
  registerDoParallel(cl)
  
  t0 <- Sys.time()
  MCsim <- foreach(
    i_MC = 1:MCreps,
    .errorhandling = "pass",
    .packages=c(
      "dplyr", "tibble", "lmtest", "sandwich", "rugarch"
    )
  )%dopar%{
    
    # Defining DGP 
    sgarch.spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
      mean.model = list(armaOrder = c(0,0), include.mean = T),
      distribution.model = "norm", 
      fixed.pars = list(mu = 0, omega=.05, alpha1=.1, beta1=.85)
    )
    
    data <- data.frame("tt" = 1:tt)
    
    path.sgarch = ugarchpath(sgarch.spec, n.sim= tt, n.start=5000, m.sim = 1)
    
    # extract returns of DGP
    data$r.sg <-  path.sgarch@path$seriesSim  
    
    # compute squared return
    data$r2.sg <- data$r.sg^2
    
    # extract volatility of DGP
    data$sig.sg <- path.sgarch@path$sigmaSim
   
    # compute variance of DGP
    data$sig2.sg <- data$sig.sg^2
    
    # Defining forecast model
    X <- NULL
    X[1] <-  data$sig2.sg[1]
    X[2] <-  data$sig2.sg[2]
    X[3] <-  data$sig2.sg[3]
    for (i in 4:nrow(data)) {
      X[i] <- (
        (1-k) + (.85/.95) * k * X[i-1] + (.1/.95)* k *(data$r.sg[i-1])^2 
      )
    }
    data$X <- X
    
    
    # Realized variance calculation
    
    #start 
    gen_rv_component <- function(m, tt) {
      
      t_save <- NULL
      lambda <- 78 / m
      for (t in 1:tt) {
        i_save <- NULL
        for (i in 1:m) {
          j_save <- NULL
          for (j in  (lambda*(i-1)+1):(lambda*i)) {
            xi_it <- rnorm( n=1,mean = 0, sd = sqrt(1/78) )
            j_save <-rbind(j_save,xi_it)
          }
          inner_<-sum(j_save)^2
          i_save <-rbind(i_save,(inner_))
        }
        i_save
        
        t_save <- rbind(t_save,sum(i_save))     
      }
      t_save
      return(t_save)
    }
    
    
    t_save1<- gen_rv_component(m =1, tt )
    RV1 <- data$sig2.sg*t_save1
    
    t_save13<- gen_rv_component(m =13, tt )
    RV13 <- data$sig2.sg*t_save13
    
    t_save78<- gen_rv_component(m =78, tt )
    RV78 <- data$sig2.sg*t_save78
    
    #end
    
    # saving results
    tibble(
      k=k,
      tt = tt,
           R2N = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$pValN,
           R2H = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$pValH,
        logR2N = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$logpValN,
        logR2H = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$logpValH,
      logR2Nnc = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$logpValNnc,
      logR2Hnc = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$logpValHnc,
        absR2N = MZ.test.mean(Proxy=data$r.sg,  forecast = X)$abspValN,
        absR2H = MZ.test.mean(Proxy=data$r.sg,  forecast = X)$abspValH,
      absR2Nnc = MZ.test.mean(Proxy=data$r.sg,  forecast = X)$abspValNnc,
      absR2Hnc = MZ.test.mean(Proxy=data$r.sg,  forecast = X)$abspValHnc,
      globFR2N = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$globalFpValN,
      globFR2H = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$globalFpValH,
      waldFR2N = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$waldFpValN,
      waldFR2H = MZ.test.mean(Proxy=data$r2.sg, forecast = X)$waldFpValH,
      
          RV1N = MZ.test.mean(Proxy=RV1, forecast = X)$pValN,
          RV1H = MZ.test.mean(Proxy=RV1, forecast = X)$pValH,
       logRV1N = MZ.test.mean(Proxy=RV1, forecast = X)$logpValN,
       logRV1H = MZ.test.mean(Proxy=RV1, forecast = X)$logpValH,
     logRV1Nnc = MZ.test.mean(Proxy=RV1, forecast = X)$logpValNnc,
     logRV1Hnc = MZ.test.mean(Proxy=RV1, forecast = X)$logpValHnc,
       absRV1N = MZ.test.mean(Proxy=RV1, forecast = X)$logpValN,
       absRV1H = MZ.test.mean(Proxy=RV1, forecast = X)$logpValH,
     absRV1Nnc = MZ.test.mean(Proxy=RV1, forecast = X)$logpValNnc,
     absRV1Hnc = MZ.test.mean(Proxy=RV1, forecast = X)$logpValHnc,
     globFRV1N = MZ.test.mean(Proxy=RV1, forecast = X)$globalFpValN,
     globFRV1H = MZ.test.mean(Proxy=RV1, forecast = X)$globalFpValH,
     waldFRV1N = MZ.test.mean(Proxy=RV1, forecast = X)$waldFpValN,
     waldFRV1H = MZ.test.mean(Proxy=RV1, forecast = X)$waldFpValH,
      
         RV13N = MZ.test.mean(Proxy=RV13, forecast = X)$pValN,
         RV13H = MZ.test.mean(Proxy=RV13, forecast = X)$pValH,
      logRV13N = MZ.test.mean(Proxy=RV13, forecast = X)$logpValN,
      logRV13H = MZ.test.mean(Proxy=RV13, forecast = X)$logpValH,
      absRV13N = MZ.test.mean(Proxy=RV13, forecast = X)$logpValN,
      absRV13H = MZ.test.mean(Proxy=RV13, forecast = X)$logpValH,
    globFRV13N = MZ.test.mean(Proxy=RV13, forecast = X)$globalFpValN,
    globFRV13H = MZ.test.mean(Proxy=RV13, forecast = X)$globalFpValH,
    waldFRV13N = MZ.test.mean(Proxy=RV13, forecast = X)$waldFpValN,
    waldFRV13H = MZ.test.mean(Proxy=RV13, forecast = X)$waldFpValH,
      
         RV78N = MZ.test.mean(Proxy=RV78, forecast = X)$pValN,
         RV78H = MZ.test.mean(Proxy=RV78, forecast = X)$pValH,
      logRV78N = MZ.test.mean(Proxy=RV78, forecast = X)$logpValN,
      logRV78H = MZ.test.mean(Proxy=RV78, forecast = X)$logpValH,
      absRV78N = MZ.test.mean(Proxy=RV78, forecast = X)$logpValN,
      absRV78H = MZ.test.mean(Proxy=RV78, forecast = X)$logpValH,
    globFRV78N = MZ.test.mean(Proxy=RV78, forecast = X)$globalFpValN,
    globFRV78H = MZ.test.mean(Proxy=RV78, forecast = X)$globalFpValH,
    waldFRV78N = MZ.test.mean(Proxy=RV78, forecast = X)$waldFpValN,
    waldFRV78H = MZ.test.mean(Proxy=RV78, forecast = X)$waldFpValH
      
    )
  }
  t1 <- Sys.time(); print(t1-t0)
  stopCluster(cl)
  
  do.call("rbind",MCsim)
}
results <- sim_func_grach(MCreps = 1, core.max = 7, k = 0.95, tt = 5) #teset zum durchlauf


####################################
# Defining values for misspecification
KK <- c(0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94, 0.95,0.96,0.97,0.98,0.99, 1)
# Defining sample size
TT=c(100,250,500,1000)

# Main simulation
res_df <- tibble()
for(tt in TT){
  for(K in KK){
    df <- sim_func_grach(
      tt=tt,
      k=K,
      MCreps=1000, core.max=7
    )
    res_df <- rbind(res_df,df)
  }
}
simgarch <- res_df

# changing data structure
Rejections1<-simgarch %>% 
  
  
  pivot_longer(
    cols=c(R2N,  R2H,   logR2N,  logR2H,   logR2Nnc,  logR2Hnc,   absR2N,     absR2H,     absR2Nnc,   absR2Hnc,   globFR2N,   globFR2H,   waldFR2N,   waldFR2H,
           RV1N, RV1H,  logRV1N, logRV1H,  logRV1Nnc, logRV1Hnc,  absRV1N,    absRV1H,    absRV1Nnc,  absRV1Hnc,  globFRV1N,  globFRV1H,  waldFRV1N,  waldFRV1H,
           RV13N,RV13H, logRV13N,logRV13H, absRV13N,  absRV13H,   globFRV13N, globFRV13H, waldFRV13N, waldFRV13H,
           RV78N,RV78H, logRV78N,logRV78H, absRV78N,  absRV78H,   globFRV78N, globFRV78H, waldFRV78N, waldFRV78H),
    names_to = "type", 
    values_to = "values"
  ) %>% 
  filter(k==0.95) %>%
  group_by(tt,type) %>%
  summarise(
    rrates = mean(as.numeric(values < .1))
  )


############Code for tables###################

Rejections2 <- Rejections1 %>%
  pivot_wider(
    names_from = type,
    values_from = rrates)
# sizes with HAC and Newey-West correction
sizeTableHAC<-Rejections2 %>% select(R2H,RV1H,RV13H,RV78H,globFR2H,globFRV1H,globFRV13H,globFRV78H,waldFR2H,waldFRV1H,waldFRV13H,waldFRV78H)
sizeTableNewey<- Rejections2 %>% select(R2N,RV1N,RV13N,RV78N,globFR2N,globFRV1N,globFRV13N,globFRV78N,waldFR2N,waldFRV1N,waldFRV13N,waldFRV78N)

# sizes of regressions on transformations with HAC and Newey-West correction
sizetabletransformationsHAC<- Rejections2 %>% select( logR2H,  logR2Hnc, logRV1H, logRV1Hnc, absR2H, absR2Hnc, absRV1H, absRV1Hnc, logRV13H, logRV78H, absRV13H, absRV78H )
sizetabletransformationsNewey<- Rejections2 %>% select( logR2N,  logR2Nnc, logRV1N, logRV1Nnc, absR2N, absR2Nnc, absRV1N, absRV1Nnc, logRV13N, logRV78N, absRV13N, absRV78N)

# calculate threshold where size is 10% for each sample size and test
sa_rates <- 
  simgarch %>% 
  pivot_longer(
    cols= c( R2N,R2H,       logR2N,  logR2H,    absR2N,   absR2H,   globFR2N,   globFR2H,   waldFR2N,   waldFR2H,
           RV1N,RV1H,     logRV1N, logRV1H,   absRV1N,  absRV1H,  globFRV1N,  globFRV1H,  waldFRV1N,  waldFRV1H,
         RV13N,RV13H,   logRV13N,logRV13H,  absRV13N, absRV13H, globFRV13N, globFRV13H, waldFRV13N, waldFRV13H,
         RV78N,RV78H,   logRV78N,logRV78H,  absRV78N, absRV78H, globFRV78N, globFRV78H, waldFRV78N, waldFRV78H),
    names_to = "type", 
    values_to = "values"
  ) %>% 
  filter(k==.95) %>% 
  group_by(tt,k,type) %>% 
  summarise(
    SArates = quantile(values, prob = .1)
  ) %>% 
  ungroup() %>% 
  group_by(tt,k,type) %>% 
  mutate(group_id = cur_group_id())%>% ungroup() %>% 
  select(group_id,SArates) 
sa_rates

# calculate size-adjusted power rates
sa_power_rates <- 
  simgarch %>% 
  pivot_longer(
    cols=c(R2N,R2H,   logR2N,  logR2H,    absR2N,   absR2H,   globFR2N,   globFR2H,   waldFR2N,   waldFR2H,
           RV1N,RV1H,     logRV1N, logRV1H,   absRV1N,  absRV1H,  globFRV1N,  globFRV1H,  waldFRV1N,  waldFRV1H,
           RV13N,RV13H,   logRV13N,logRV13H,  absRV13N, absRV13H, globFRV13N, globFRV13H, waldFRV13N, waldFRV13H,
           RV78N,RV78H,   logRV78N,logRV78H,  absRV78N, absRV78H, globFRV78N, globFRV78H, waldFRV78N, waldFRV78H),
    names_to = "type", 
    values_to = "values"
  ) %>% 
  group_by(tt,type) %>% 
  mutate(group_id = cur_group_id()) %>% 
  left_join(.,sa_rates ) %>% 
  group_by(tt,k,type) %>% 
  summarise(
    rejr = mean(as.numeric(values < SArates))
  )
sa_power_rates



############ Code for graphical illustrations##############################

# Create powerfunctions for MZ-Classic, MZ-F and MZ-F* with T= {100, 500, 1000} and RV1,RV13 and RV78
datalong_bigplot1 <- sa_power_rates %>%
  filter(tt %in% c(100,500,1000),type %in% c("RV1H", "RV13H", "RV78H", "globFRV1H", "globFRV13H", "globFRV78H", "waldFRV1H", "waldFRV13H", "waldFRV78H")) %>%
  mutate(MZtypes = case_when(
    type %in% c("RV1H", "RV13H", "RV78H") ~ "MZ-Classic",
    type %in% c("globFRV1H", "globFRV13H", "globFRV78H") ~ "MZ-F",
    type %in% c("waldFRV1H", "waldFRV13H", "waldFRV78H") ~ "MZ-F*",
  ))%>%
  mutate(color = case_when(
    type %in% c("RV1H", "globFRV1H", "waldFRV1H") ~ "RV1",
    type %in% c("RV13H", "globFRV13H", "waldFRV13H") ~ "RV13",
    type %in% c("RV78H", "globFRV78H", "waldFRV78H") ~ "RV78",
  ))%>%
  mutate(marker = case_when(
    k == 0.98 & type %in% c("RV1H", "globFRV1H", "waldFRV1H") ~ "point",
    k == 0.98 & type %in% c("RV13H", "globFRV13H", "waldFRV13H") ~ "square",
    k == 0.98 & type %in% c("RV78H", "globFRV78H", "waldFRV78H") ~ "triangle",
    TRUE ~ "Standard"
  ))
library(ggplot2)
library(RColorBrewer)

color_scheme <- brewer.pal(3, "Dark2")


points_data <- datalong_bigplot1 %>%
  filter(marker == "point")
square_data <- datalong_bigplot1 %>%
  filter(marker == "square")
triangle_data <- datalong_bigplot1 %>%
  filter(marker == "triangle")



powerfunctionsRV<-
  ggplot(datalong_bigplot1, aes(x = k, y = rejr, color = color)) +
  geom_line() +
  geom_point(data = points_data) +
  geom_point(data = square_data, shape = 15) +
  geom_point(data = triangle_data, shape = 17) +
  facet_grid(MZtypes ~ tt, switch = "y") +
  ylim(0, 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0.95,  color = "darkgrey") +
  scale_color_manual(values = color_scheme)+
  labs(x = "k", y = "size adj. power")+
  theme_bw() +
  theme(legend.position = "none") 

ggsave(filename = "C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/powerfunctionsRV.png", plot = powerfunctionsRV, dpi = 600, width = 8, height = 5)




####################################
# Create powerfunctions for MZ-Classic, MZ-F and MZ-F* with T= {500, 1000} and RV1 and squared return

datalong_bigplot1RV1R2 <- sa_power_rates %>%
  filter(tt %in% c(500,1000),type %in% c("RV1H", "R2H", "globFRV1H", "globFR2H", "waldFRV1H", "waldFR2H")) %>%
  mutate(MZtypes = case_when(
    type %in% c("RV1H", "R2H") ~ "MZ-Classic",
    type %in% c("globFRV1H", "globFR2H") ~ "MZ-F",
    type %in% c("waldFRV1H", "waldFR2H") ~ "MZ-F*",
    
  ))%>%
  mutate(color = case_when(
    type %in% c("RV1H", "globFRV1H", "waldFRV1H") ~ "RV1",
    type %in% c("R2H", "globFR2H", "waldFR2H") ~ "R2",
  ))%>%
  mutate(marker = case_when(
    k == 0.98 & type %in% c("RV1H", "globFRV1H", "waldFRV1H") ~ "point",
    k == 0.98 & type %in% c("R2H", "globFR2H", "waldFR2H") ~ "square",
    TRUE ~ "Standard"
  ))

library(viridis)
color_scheme <- viridis_pal(option = "D")(3)


points_data <- datalong_bigplot1RV1R2 %>%
  filter(marker == "point")
square_data <- datalong_bigplot1RV1R2 %>%
  filter(marker == "square")




powerfunctionsR2RV<-
  ggplot(datalong_bigplot1RV1R2, aes(x = k, y = rejr, color = color)) +
  geom_line() +
  geom_point(data = points_data ) +
  geom_point(data = square_data, shape = 15 ) +
  facet_grid(MZtypes ~ tt, switch = "y") +
  ylim(0, 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0.95,  color = "darkgrey") +
  scale_color_manual(values = color_scheme, 
                     labels = c("RV1" = "RV1",
                                "RV2" = "RV2"),
                     name = NULL,
                     guide = guide_legend(override.aes = list(shape = c(15, 16))))+
  labs( x = "k", y = "size adj. power", shape = "marker") +
  theme_bw() +
  theme(strip.placement = "outside",
        legend.position = "none")

ggsave(filename = "C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/powerfunctionsR2RV.png", plot = powerfunctionsR2RV, dpi = 600, width = 8, height = 5)



##################################
# Create powerfunctions of MZ regressions on log transformations

datalong_biglog <- sa_power_rates %>%
  filter(tt %in% c(1000),type %in% c("logRV1H", "logR2H", "RV1H", "R2H" )) %>%
  mutate(MZtypes = case_when(
    type %in% c("logR2H", "R2H") ~ "squared return",
    type %in% c("logRV1H","RV1H") ~ "RV1",
  ))%>%
  mutate(color = case_when(
    type %in% c("logRV1H","logR2H")~"log",
    type %in% c("RV1H", "R2H")  ~ "nolog",
  ))%>%
  mutate(marker = case_when(
    k == 0.98 & type %in% c("RV1H", "R2H") ~ "point",
    k == 0.98 & type %in% c("logRV1H", "logR2H") ~ "square",
    TRUE ~ "Standard"
  ))

library(RColorBrewer)


color_scheme <- viridis_pal(option = "D")(3)


points_data <- datalong_biglog %>%
  filter(marker == "point")
square_data <- datalong_biglog %>%
  filter(marker == "square")



powerfunctionslog<-
  ggplot(datalong_biglog, aes(x = k, y = rejr, color = color)) +
  geom_line() +
  geom_point(data = points_data ) +
  geom_point(data = square_data, shape = 15 ) +
  facet_grid(MZtypes~., switch = "y") +
  ylim(0, 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0.95,  color = "darkgrey") +
  scale_color_manual(values = color_scheme, 
                     name = NULL)+
  labs( x = "k", y = "size adj. power", shape = "marker") +
  theme_bw() +
  theme(strip.placement = "outside",
        legend.position = "none")

ggsave(filename = "C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/powerfunctionslog.png", plot = powerfunctionslog, dpi = 600, width = 6, height = 6)
######################################################################
# Create powerfunctions of MZ regressions on forecast^0.5

datalong_bigabs <- sa_power_rates %>%
  filter(tt %in% c(1000),type %in% c("absRV1H", "absR2H", "RV1H", "R2H" )) %>%
  mutate(MZtypes = case_when(
    type %in% c("absR2H", "R2H") ~ "squared return",
    type %in% c("absRV1H","RV1H") ~ "RV1",
  ))%>%
  mutate(color = case_when(
    type %in% c("absRV1H","absR2H")~"abs",
    type %in% c("RV1H", "R2H")  ~ "noabs",
  ))%>%
  mutate(marker = case_when(
    k == 0.98 & type %in% c("RV1H", "R2H") ~ "point",
    k == 0.98 & type %in% c("absRV1H", "absR2H") ~ "square",
    TRUE ~ "Standard"
  ))

color_scheme <- viridis_pal(option = "D")(3)

points_data <- datalong_bigabs %>%
  filter(marker == "point")
square_data <- datalong_bigabs %>%
  filter(marker == "square")



powerfunctionsabs<-
  ggplot(datalong_bigabs, aes(x = k, y = rejr, color = color)) +
  geom_line() +
  geom_point(data = points_data ) +
  geom_point(data = square_data, shape = 15 ) +
  facet_grid(MZtypes~ ., switch = "y") +
  ylim(0, 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0.95,  color = "darkgrey") +
  scale_color_manual(values = color_scheme, 
                     name = NULL)+
  labs( x = "k", y = "size adj. power", shape = "marker") +
  theme_bw() +
  theme(strip.placement = "outside",
        legend.position = "none")

ggsave(filename = "C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/powerfunctionsabs.png", plot = powerfunctionsabs, dpi = 600, width = 6, height = 6)

########################################################################
# Create powerfunctions for MZ-Classic, MZ-F and MZ-F* with T={100, 1000} and RV1,RV13 and RV 78 with other dimensions

datalong_bigplot2 <- sa_power_rates %>%
  filter(tt %in% c(100,1000),type %in% c("RV1H", "RV13H", "RV78H", "globFRV1H", "globFRV13H", "globFRV78H", "waldFRV1H", "waldFRV13H", "waldFRV78H")) %>%
  mutate(RVtypes = case_when(
    type %in% c("RV1H", "globFRV1H", "waldFRV1H") ~ "RV1",
    type %in% c("RV13H","globFRV13H","waldFRV13H") ~ "RV13",
    type %in% c("RV78H","globFRV78H","waldFRV78H") ~ "RV78",
  ))%>%
  mutate(color = case_when(
    type %in% c("RV1H", "RV13H", "RV78H") ~ "MZ-Classic",
    type %in% c("globFRV1H", "globFRV13H", "globFRV78H") ~ "MZ-F",
    type %in% c("waldFRV1H", "waldFRV13H", "waldFRV78H") ~ "MZ-F*",
  ))%>%
  mutate(marker = case_when(
    k == 0.97 & type %in% c("RV1H", "RV13H", "RV78H") ~ "point",
    k == 0.98 & type %in% c("globFRV1H", "globFRV13H", "globFRV78H") ~ "square",
    k == 0.99 & type %in% c("waldRV1H", "waldRV13H", "waldFRV78H") ~ "triangle",
    TRUE ~ "Standard"
  ))

library(RColorBrewer)


color_scheme <- brewer.pal(3, "Dark2")


points_data <- datalong_bigplot2 %>%
  filter(marker == "point")
square_data <- datalong_bigplot2 %>%
  filter(marker == "square")
triangle_data <- datalong_bigplot2 %>%
  filter(marker == "triangle")




powerfunctionstest<-
  ggplot(datalong_bigplot2, aes(x = k, y = rejr, color = color)) +
  geom_line() +
  geom_point(data = points_data ) +
  geom_point(data = square_data, shape = 15 ) +
  geom_point(data = triangle_data, shape = 17  ) +
  
  facet_grid(RVtypes ~ tt, switch = "y") +
  ylim(0, 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = color_scheme, 
                     labels = c("MZ-Classic" = "MZ_Classic",
                                "MZ-F" = "MZ-F",
                                "MZ-F*" = "MZ-F*"),
                     name = NULL,
                     guide = guide_legend(override.aes = list(shape = c(16, 15, 17))))+
  
  labs( x = "k", y = "size adj. power", shape = "marker") +
  theme_bw() +
  theme(strip.placement = "outside",
        legend.position = "none") 

ggsave(filename = "C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/powerfunctionstest.png", plot = powerfunctionstest, dpi = 600, width = 8, height = 5)


#############################################################
#Create powerfunctions for comparison of Newey-West and HAC standard error correction with T={100, 1000}

datalong_bigplotnewey<- sa_power_rates %>% filter(type %in% c("RV1N","RV13N","RV78N","RV1H","RV13H","RV78H"))


datalong_bigplotnewey <- sa_power_rates %>%
  filter(tt %in% c(100,1000),type %in% c("RV1N","RV13N","RV78N","RV1H","RV13H","RV78H")) %>%
  mutate(MZtypes = case_when(
    type %in% c("RV1H", "RV1N" ) ~ "RV1",
    type %in% c("RV13H", "RV13N") ~ "RV13",
    type %in% c("RV78H", "RV78N") ~ "RV78",
  ))%>%
  mutate(color = case_when(
    type %in% c("RV1H", "RV13H","RV78H") ~ "Andrews",
    type %in% c("RV1N", "RV13N","RV78N") ~ "Newey-West",
  ))%>%
  mutate(marker = case_when(
    k == 0.98 & type %in% c("RV1H", "RV13H","RV78H") ~ "point",
    k == 0.98 & type %in% c("RV1N", "RV13N","RV78N") ~ "square",
    TRUE ~ "Standard"
  ))

color_scheme <- viridis_pal(option = "D")(3)


points_data <- datalong_bigplotnewey %>%
  filter(marker == "point")
square_data <- datalong_bigplotnewey %>%
  filter(marker == "square")


powerfunctionsRVnewey<-
  ggplot(datalong_bigplotnewey, aes(x = k, y = rejr, color = color)) +
  geom_line() +
  geom_point(data = points_data) +
  geom_point(data = square_data, shape = 15) +
  facet_grid(tt ~ MZtypes, switch = "y") +
  ylim(0, 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0.95,  color = "darkgrey") +
  scale_color_manual(values = color_scheme)+
  labs(x = "k", y = "size adj. power")+
  theme_bw() +
  theme(legend.position = "none") 

ggsave(filename = "C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/powerfunctionsRVnewey.png", plot = powerfunctionsRVnewey, dpi = 600, width = 8, height = 5)












