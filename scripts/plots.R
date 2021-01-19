library(magrittr)
library(plyr)
library(ggplot2)
library(patchwork)
library(reshape2)

# Reading in output from julia code
data <- readRDS("output/data.rds") %>% 
  arrange(sample(1:100, replace = FALSE)) %>% 
  mutate(lat = rep(1:10, 10),
         lon = rep(1:10, each = 10),
         y = as.numeric(y))

post <- readRDS("output/chain.rds") %>% 
  melt(varnames = c("loc", "iter"), value.name = "p") %>% 
  ddply(.(loc), summarise,
        med = median(p),
        lower = quantile(p, 0.025),
        upper = quantile(p, 0.975)) %>% 
  join(data)

xopt <- data$x[which.min(data$score)]
dopt <- data$loc[which.min(data$score)]

# Plots ========================================================================

# Initial fit
post %>% 
  ggplot(aes(x, med)) + 
  geom_line() + 
  geom_ribbon(aes(min = lower, max = upper), alpha = 0.2) + 
  geom_point(data = subset(data, K1 == 1), aes(x, as.numeric(y), shape = factor(y)), size = 2) + 
  scale_shape_manual(values = c(1, 19), guide = FALSE) + 
  geom_vline(xintercept = xopt, linetype = 2) +
  theme_classic() +
  ylab("p") +
  ggtitle("a") -> fig1a

# Initial sampling locations
data %>% 
  ggplot(aes(lon, lat, fill = x)) + 
  geom_raster() +
  scale_fill_distiller(type = "div", palette = "BrBG") + 
  geom_point(aes(lon, lat, shape = factor(y)), data = subset(data, K1 == 1), size = 3) +
  scale_shape_manual("y", values = c(1, 19)) + 
  theme_classic() + 
  scale_y_continuous("northing", expand = c(0, 0)) + 
  scale_x_continuous("easting", expand = c(0, 0)) + 
  theme(axis.text = element_blank()) + 
  ggtitle("b") -> fig1b

# Scores by x
data %>% 
  ggplot(aes(x, score)) + 
  geom_line(alpha = 0.5) +
  geom_point(aes(color = score)) + 
  scale_color_viridis_c(option = "B", guide = FALSE) + 
  geom_vline(xintercept = xopt, linetype = 2) + 
  geom_rug(data = subset(data, K1 == 1), aes(x = x), inherit.aes = FALSE) +
  ylab("design criterion") + 
  theme_classic() + 
  ggtitle("c") -> fig1c

# Map of scores
data %>% 
  ggplot(aes(lon, lat, fill = score)) + 
  geom_raster() +
  scale_fill_viridis_c("design\ncriterion", option = "B") + 
  geom_point(aes(lon, lat, shape = factor(y)), data = subset(data, K1 == 1), size = 3, inherit.aes = FALSE) +
  scale_shape_manual(values = c(1, 19), guide = FALSE) + 
  geom_point(aes(lon, lat), data = subset(data, loc == dopt), color = "white", shape = 4, size = 3) + 
  theme_classic() + 
  scale_y_continuous("northing", expand = c(0, 0)) + 
  scale_x_continuous("easting", expand = c(0, 0)) + 
  theme(axis.text = element_blank()) + 
  ggtitle("d") -> fig1d

pdf(file = "output/fig1.pdf", width = 7, height = 6)

fig1a + fig1b + fig1c + fig1d + plot_layout(ncol = 2)

dev.off()  
