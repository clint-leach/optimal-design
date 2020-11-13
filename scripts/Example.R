#########################################################################
### Basic Optimal Design Example - Proof of Concept
### Perry Williams & Mevin Hooten
### Computation time with iMac Pro (2017) with 3.0 GHz, 128 GB Mem, 20 cores:
###               18.6 secs : 15 potential designs,   100,000 MCMC samples
###               2.50 mins : 15 potential designs, 1,000,000 MCMC samples
#########################################################################

# Load/Install required packages and functions =================================
library(doParallel)
library(magrittr)
library(plyr)
library(ggplot2)
library(msm)
library(patchwork)

# GLM - probit pprb
source("src/probit.pprb.R")

# GLM - probit link
source("src/probit.reg.mcmc.R")

# Simulate data ================================================================

# Set the RNG seed
set.seed(1701)

# Size of full domain 
N <- 100

# Size of initial sample
n <- 10

# Parameters
betatrue <- matrix(c(0.0, 0.5), nrow = 2)#rnorm(2) %>% matrix(2)

# Generating data
data <- data.frame(loc = 1:N, x = rnorm(N, 0, 3)) %>%
  mutate(mu = rnorm(N, betatrue[1] + betatrue[2] * x, 1),
         p = pnorm(betatrue[1] + betatrue[2] * x),
         y = ifelse(mu >= 0, 1, 0),
         K1 = loc %in% sample(1:N, n),
         lat = rep(1:10, 10),
         lon = rep(1:10, each = 10))

# Design matrix for matrix math
X <- cbind(1, data$x)

# Fit model ====================================================================

# Beta priors
beta.mn <- matrix(0, 2, 1)
beta.var <- 2.25

# MCMC
n.iter <- 100000

chain <- probit.reg.mcmc(subset(data, K1 == 1),
                         beta.mn,
                         beta.var,
                         n.iter)

beta.samp <- data.frame(chain$beta.save)
names(beta.samp) <- c("beta0", "beta1")
beta.samp$iter <- c(1:n.iter)

# Calculate posterior predictive distributions
ppd <- ddply(beta.samp, .(iter), summarise,
             mu = beta0 + beta1 * data$x,
             mu.ppd = rnorm(N, mu, 1),
             yhat = ifelse(mu.ppd >= 0, 1, 0),
             p = pnorm(mu),
             x = data$x,
             loc = c(1:N))

variances <- ddply(ppd, .(loc), summarise, 
                   vary = var(yhat), 
                   varcond = mean(p * (1 - p)) + var(p),
                   x = mean(x), 
                   p = mean(p))

ggplot(variances, aes(x, vary)) + 
  geom_point() + 
  geom_line(aes(x, varcond))

# Identifying optimal design ===================================================

# Number of new samples

n.new <- 1

# Unsampled sites

unsampled <- which(data$K1 == 0)

# All possible combinations of unsampled sites (i.e., all designs)

D <- combn(unsampled, n.new)
M <- dim(D)[2]

# Generating new data
pred <- ddply(ppd, .(loc), summarise,
                p.pred = mean(yhat),
                y.pred = ifelse(p.pred > 0.5, 1, 0),
                score = NA)

# Loop over all designs
score <- rep(NA, M)
registerDoParallel(cores=detectCores())

updates <- foreach(d=1:M, .export=c('probit.pprb')) %do% {

    cat(round(d/M*100),"% - ")

    # Update model with ypred
    y.pred <- subset(pred, loc == D[d])$y.pred
    X.pred <- subset(data, loc == D[d])$x
    updated <- probit.pprb(y.pred, X.pred, chain$beta.save)
    
    # Compute expected total predictive variance
    p <- X %*% t(updated$beta.save) %>%
      pnorm()
    
    vary <- rep(0, N)
    
    for(i in 1:N){
      yhat <- rbinom(n.iter, 1, p[i, ])
      vary[i] <- var(yhat)
    }
    
    # meanvar <- (p * (1 - p)) %>% apply(1, mean)
    # varmean <- p %>% apply(1, var)
    # score[d] <- sum(meanvar + varmean)
    
    score[d] <- sum(vary)
    
    return(updated$beta.save)
}

# Proccess output ==============================================================

# Selecting optimal design
scores <- data[D, ] %>% cbind(score)
dopt <- which.min(scores$score)
xopt <- scores$x[which.min(scores$score)]
opt <- scores[which.min(scores$score), ]

# Processing posterior of p and related functions of uncertainty
ppd %>% 
  mutate(var = p * (1 - p)) %>% 
  ddply(.(loc), summarise, 
        med = median(p),
        low = quantile(p, 0.025),
        high = quantile(p, 0.975),
        range = high - low,
        varp = var(p),
        var = mean(var)) %>% 
  join(scores) -> scores

# Plots ========================================================================

# Initial fit
ppd %>% 
  ggplot(aes(x, p)) + 
  stat_summary(alpha = 0.2, color = "gray",  
               fun.min = function(z) quantile(z, 0.025), 
               fun.max = function(z) quantile(z, 0.975),
               geom = "ribbon") +
  stat_summary(fun = "median", geom = "line") + 
  geom_point(data = subset(data, K1 == 1), aes(x, y, shape = factor(y)), size = 2) + 
  scale_shape_manual(values = c(1, 19), guide = FALSE) + 
  geom_vline(xintercept = xopt, linetype = 2) +
  theme_classic() +
  ggtitle("A") -> fig1a

# Initial sampling locations
data %>% 
  ggplot(aes(lon, lat, fill = x)) + 
  geom_raster() +
  scale_fill_viridis_c(option = "E") + 
  geom_point(aes(lon, lat, shape = factor(y)), data = subset(data, K1), size = 3) +
  scale_shape_manual("y", values = c(1, 19)) + 
  theme_classic() + 
  scale_y_continuous("northing", expand = c(0, 0)) + 
  scale_x_continuous("easting", expand = c(0, 0)) + 
  theme(axis.text = element_blank()) + 
  ggtitle("B") -> fig1b

# Scores by x
scores %>% 
  ggplot(aes(x, score)) + 
  geom_line(alpha = 0.5) +
  geom_point(aes(color = score)) + 
  scale_color_viridis_c(option = "B", guide = FALSE) + 
  geom_vline(xintercept = xopt, linetype = 2) + 
  geom_rug(data = subset(data, K1 == 1), aes(x = x), inherit.aes = FALSE) +
  ylab("design criterion") + 
  theme_classic() + 
  ggtitle("C") -> fig1c

# Map of scores
scores %>% 
  ggplot(aes(lon, lat, fill = score)) + 
  geom_raster() +
  scale_fill_viridis_c("design\ncriterion", option = "B") + 
  geom_point(aes(lon, lat, shape = factor(y)), data = subset(data, K1), size = 3, inherit.aes = FALSE) +
  scale_shape_manual(values = c(1, 19), guide = FALSE) + 
  geom_point(aes(lon, lat), data = subset(scores, loc == D[dopt]), color = "white", shape = 4, size = 3) + 
  theme_classic() + 
  scale_y_continuous("northing", expand = c(0, 0)) + 
  scale_x_continuous("easting", expand = c(0, 0)) + 
  theme(axis.text = element_blank()) + 
  ggtitle("D") -> fig1d

pdf(file = "output/fig1.pdf", width = 7, height = 6)

fig1a + fig1b + fig1c + fig1d + plot_layout(ncol = 2)

dev.off()
