
######## test true slope recoverability via logistic regression and LDA #########


library(tidyverse)
library(doParallel)
select <- dplyr::select


# this function simulates a set of presences and absences 
# based on a gaussian climatic niche and a "slope" coefficient 
# describing how northness modifies macroclimate
simulateDataset <- function(n, # number of sites
                            slope, # true slope of macroclimate-northness relationship
                            niche_mean, # species niche optimum
                            niche_sd, # species nice breadth
                            prevalence){ # occupancy rate of suitable sites
      
      macroclimate <- runif(n, 0, 30)
      northness <- rnorm(n, 0, .35)
      northness[abs(northness)>1] <- 0
      microclimate <- macroclimate + northness * slope
      
      suitability <- function(x, niche_mean, niche_sd, prevalence){
            y <- dnorm(x, niche_mean, niche_sd) / 
                  dnorm(niche_mean, niche_mean, niche_sd)
            y <- y / mean(y) * prevalence
            return(y)
      }
      
      suitable <- suitability(microclimate, niche_mean, niche_sd, prevalence)
      occur <- rbinom(n, 1, suitable)
      
      return(data.frame(macroclimate,
                        northness,
                        occur))
}




# get range of sample sizes representing actual FIA dataset
d <- readRDS("e:/fia/topoclimate/data_munged.rds")
d <- d %>% filter(lon < -95) %>% select(plt_cn, subp, genus, species) %>% distinct()
abundances <- group_by(d, genus, species) %>% summarize(n=n())
n <- length(unique(paste(d$plt_cn, d$subp)))
prevalences <- quantile(abundances$n, c(0, .1, .5, .9, 1)) / n

# simulation parameters
params <- expand.grid(n=n,
                      slope=c(1, 3, 5),
                      niche_mean=15,
                      niche_sd=c(1, 3, 5),
                      prevalence=prevalences,
                      quantile=c(.5, .75, .9, .95),
                      trim=c(T, F),
                      rep=1:5) %>% #100
      mutate(i = 1:nrow(.))

cl <- makeCluster(detectCores())
registerDoParallel(cl)
r <- foreach(i = 1:nrow(params), .combine="rbind",
             .packages="tidyverse", .noexport="d") %dopar% {
                   x <- params[i,]
                   
                   df <- simulateDataset(x$n, x$slope, x$niche_mean, x$niche_sd, x$prevalence)
                   
                   if(x$trim) df <- filter(df, abs(northness) < .5)
                   pres <- df %>% filter(occur==1)
                   if(nrow(pres)<100) return(NA)
                   df <- df %>% filter(ecdf(pres$macroclimate)(macroclimate) > x$quantile)
                   pres <- df %>% filter(occur==1)
                   if(nrow(pres)<100) return(NA)
                   
                   fit <- glm(occur ~ macroclimate + northness, data=df,
                              family=binomial(link="logit"))
                   x$fit <- coef(fit)["northness"] / coef(fit)["macroclimate"]
                   return(x)
             }
stopCluster(cl)

vars <- c("slope", "niche_sd", "prevalence", 
          "quantile", "trim")
pd <- ecoclim::pairsData(r %>% mutate(err=abs(fit-slope), 
                                      prevalence=round(prevalence, 3), 
                                      err_ratio=abs(log10(fit/slope))) %>% 
                               as.data.frame(),
                         vars, 
                         c("err", "err_ratio", "rep"), mirror=T) %>%
      group_by(x_var, y_var, x_value, y_value) %>%
      summarize(err=mean(err, na.rm=T),
      err_ratio=mean(err_ratio, na.rm=T)) %>% ungroup() %>%
      mutate(x_value=factor(x_value), y_value=factor(y_value)) %>%
      filter(!is.na(x_value), !is.na(y_value)) %>%
      mutate(x_var=factor(x_var, levels=vars), 
             y_var=factor(y_var, levels=vars))

p <- ggplot(pd, aes(x_value, y_value, fill=err)) + 
      geom_tile() +
      facet_grid(y_var ~ x_var, scales="free") +
      scale_fill_gradientn(colours=c("black", "blue", "red", "yellow"), trans="log10") +
      labs(fill="slope MAE  ", x="parameter value", y="parameter value") +
      theme_minimal() + theme(legend.position="top")
ggsave("figures/niche_simulation_logistic_MAE_heatmap.png", width=5, height=5.5, units="in")



p <- ggplot(pd, aes(x_value, y_value, fill=err_ratio)) + 
      geom_tile() +
      facet_grid(y_var ~ x_var, scales="free") +
      scale_fill_gradientn(colours=c("black", "blue", "red", "yellow")) +
      labs(fill="slope mean abs log10 error ratio\n(<0.1 is within +/- 10% of true slope) ", x="parameter value", y="parameter value") +
      theme_minimal() + theme(legend.position="top")
ggsave("figures/niche_simulation_logistic_MAE_ratio_heatmap.png", width=5, height=5.5, units="in")

