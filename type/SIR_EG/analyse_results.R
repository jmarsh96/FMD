# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(dplyr)
library(ggplot2)

ReturnMedianHistogram <- function(res, cur_parameter, type, truth) {
  p <- res %>% filter(parameter == cur_parameter, type == type) %>% 
    group_by(run) %>% summarise(median = median(value)) %>% 
    ggplot(aes(x=median)) + 
    geom_histogram(aes(y=..density..), fill = "white", color="black") + 
    xlab(cur_parameter) + geom_vline(aes(xintercept = truth), colour="red") + 
    labs(title=paste0("Posterior median estimates for ", cur_parameter))
}


PlotMedians <- function(res_long, type, truth) {
  require(gridExtra)
  p1 <- ReturnMedianHistogram(res, "beta", type, truth[1])
  p2 <- ReturnMedianHistogram(res, "gamma", type, truth[2])
  p3 <- ReturnMedianHistogram(res, "lambda", type, truth[3])
  p4 <- ReturnMedianHistogram(res, "inf_error", type, truth[4])
  grid.arrange(p1,p2,p3,p4)
}


res <- readRDS("results/1/res_long.rds")
truth <- c(0.002, 0.5, 4*10^(-6), 0)
PlotMedians(res,1,truth)





res %>% filter(parameter == "beta") %>% 
  group_by(run) %>% summarise(median = median(value)) %>% 
  ggplot(aes(x=median)) + geom_histogram()



res %>% filter(parameter == "inf_error") %>% 
  group_by(run) %>% summarise(median = median(value)) %>% 
  ggplot(aes(x=median)) + geom_histogram()


specific_run <- 5
res %>% 
  filter(run == specific_run, parameter == "inf_error") %>% 
  ggplot(aes(x=iteration,y=value)) + geom_line()



res %>% filter(parameter == "beta", type == type) %>% 
  group_by(run) %>% summarise(median = median(value)) %>% 
  ggplot(aes(x=median)) + geom_histogram(aes(y=..density..),fill="white",color="black") + xlab("beta") +
  geom_vline(aes(xintercept = 0.002),colour="red") + title("a")





