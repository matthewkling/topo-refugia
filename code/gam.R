library(tidyverse)
library(mgcv)

# simulate hump-shaped data, proportion ~ cwd
n <- 1000
d <- data.frame(cwd = rnorm(n)) %>%
      mutate(p = 1 - scales::rescale(abs(cwd) + rnorm(n)))

# fit GAM smoothing spline
fit <- gam(p ~ s(cwd, k=2), data=d, family=binomial(logit))
d$fit <- predict(fit, d, type="response")

# plot
ggplot(d) + 
      geom_point(aes(cwd, p)) +
      geom_line(aes(cwd, fit), color="red")
