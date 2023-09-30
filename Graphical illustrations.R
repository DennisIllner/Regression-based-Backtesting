library(rugarch)
library(dplyr)
library(sandwich)
library(car)
library(ggplot2)
library(gridExtra)
library(patchwork)

# specify DGP
sgarch.spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = T),
  distribution.model = "norm", 
  fixed.pars = list(mu = 0, omega=.05, alpha1=.1, beta1=.85)
)
# specify sample size 
tt=250

# generate returns
path.sgarch = ugarchpath(sgarch.spec, n.sim= tt, n.start=6000, m.sim = 1)

# extract returns
r.sg <-  path.sgarch@path$seriesSim  

# compute squared return
r2.sg <- r.sg^2   

# extract volatility of DGP     
sig.sg <- path.sgarch@path$sigmaSim 

# compute variance
sig2.sg <- sig.sg^2                   


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
        #(lambda*(i-1)+1)
        #(lambda*i) 
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

#Computation of Realized Variances
t_save1<- gen_rv_component(m =1, tt )
RV1 <- sig2.sg*t_save1

t_save13<- gen_rv_component(m =13, tt )
RV13 <- sig2.sg*t_save13

t_save78<- gen_rv_component(m =78, tt )
RV78 <- sig2.sg*t_save78

#end

# Defining forecast model
calculate_X <- function(k) {
  X <- numeric(tt)
  X[1:3] <- sig2.sg[1:3]
  
  for (i in 4:tt) {
    X[i] <- (1 - k) + (.85 / .95) * k * X[i - 1] + (.1 / .95) * k * (r.sg[i - 1])^2
  }
  
  return(X)
}

# computation of forecasts
X0.95<-calculate_X(0.95)
X0.80<-calculate_X(0.8)
X1.00<-calculate_X(1.05)



# MZ regression
MZreg <- lm(r2.sg ~ X0.95)
MZreg %>% summary



# Plots for MZ regression with correct and misspecified forecast model
plot1 <- ggplot(data.frame(r2.sg, X0.95), aes(x = X0.95, y = r2.sg)) +
  geom_point(color = "dodgerblue") +
  geom_abline(intercept = coef(lm(r2.sg ~ X0.95))[1], slope = coef(lm(r2.sg ~ X0.95))[2], color = "darkorange") +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(x = "Forecast", y = "squared return") +
  theme_minimal() +
  labs(x = "Forecast", y = "squared return")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  theme(
    panel.background = element_rect(fill = NA),  
    plot.background = element_rect(fill = NA),   
    panel.grid.major = element_line(color = "gray", size = 0.5)
  )

plot2 <- ggplot(data.frame(r2.sg, X0.80), aes(x = X0.80, y = r2.sg)) +
  geom_point(color = "dodgerblue") +
  geom_abline(intercept = coef(lm(r2.sg ~ X0.80))[1], slope = coef(lm(r2.sg ~ X0.80))[2], color = "darkorange") +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(x = "Forecast", y = "squared return") +
  theme_minimal() +
  labs(x = "Forecast", y = "squared return")+
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  theme(
    panel.background = element_rect(fill = NA),  
    plot.background = element_rect(fill = NA),   
    panel.grid.major = element_line(color = "gray", size = 0.5)
  )


MZplot <- grid.arrange(plot1, plot2, ncol = 2, widths = c(1, 1))
ggsave(filename ="C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/MZplot.png", plot = MZplot, dpi = 600, height = 5, width = 10)


  
# Plot for time series of returns and forecasts with varying degrees of persistence  
  Misspecifications <- ggplot() +
      geom_line(aes(x = 1:tt, y = r.sg), color = "darkgrey") +
      geom_line(aes(x = 1:tt, y = X0.95), color = "darkblue", size = 1.3) +
      geom_line(aes(x = 1:tt, y = X1.00), color = "darkorange", size = 0.8) +
      geom_line(aes(x = 1:tt, y = X0.80), color = "purple", size = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_point(data = data.frame(time = 245, x0.95 = X0.95[245]), aes(x = time, y = x0.95), color = "darkblue", size = 4) + 
      geom_point(data = data.frame(time = 245, x1.00 = X1.00[245]), aes(x = time, y = x1.00), color = "darkorange", size = 4, shape=15) +# Punkt für X0.95 bei tt = 200
      geom_point(data = data.frame(time = 245, x0.80 = X0.80[245]), aes(x = time, y = x0.80), color = "purple", size = 4, shape=17) +
      labs(x = "Time", y = "return and forecasts") +
      scale_color_manual(values = c("darkgrey", "darkblue", "darkorange", "purple")) +
      theme_minimal() +
      guides(color = "none")
  ggsave(filename ="C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/Misspecifications.png", plot = Misspecifications, dpi = 600, height = 4, width = 8 )
  
# Plot for time series of returns and different Volatility Proxies 
  TSProxys<-ggplot() +
      geom_line(aes(x = 1:tt, y = r.sg), color = "darkgrey") +
      geom_line(aes(x = 1:tt, y = RV1), color = "darkblue", size = 0.75) +
      geom_line(aes(x = 1:tt, y = RV78), color = "darkorange", size = 0.75) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_point(data = data.frame(time = 245, rv1 = RV1[245]), aes(x = time, y = rv1), color = "darkblue", size = 4) + 
      geom_point(data = data.frame(time = 235, rv78 = RV78[235]), aes(x = time, y = rv78), color = "darkorange", size = 4, shape=15) +
       labs(x = "Time", y = "return and volatility proxys") +
      scale_color_manual(values = c("darkgrey", "darkblue", "orange", "darkorange")) +
      theme_minimal() +
      guides(color = "none")
  ggsave(filename ="C:/Users/illne/Desktop/Uni/6.Semester/Bachelorarbeit/Text/TSProxys.png", plot = TSProxys, dpi = 600, height = 4, width = 8 )
  
  
  
  
  
  
  

    
  
  
  
  
  
 