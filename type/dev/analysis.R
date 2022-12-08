rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

exposure_data <- readRDS("res/exposure_matrix.rds")





pdf("exp_sum_trace.pdf")
par(mfrow=c(3,1))
errors <- numeric(200)
for(i in 1:200) {
  if(!is.na(exposure_data$truth[i])) {
    cur_fs <- exposure_data$final_size[i]
    plot(exposure_data$exposure_matrix[i,],type="l",ylab="Exposure sum",main=paste0("run = ",i, ", final size = ", cur_fs))
    abline(h=exposure_data$truth[i],col=2)
  }
}
dev.off()


relative_error <- (apply(exposure_data$exposure_matrix, 1, median) - exposure_data$truth)/exposure_data$truth

pdf("error_high_beta.pdf")
hist(relative_error, main="high_beta")
plot(exposure_data$final_size,relative_error, xlab="Final size", ylab="Relative error", main="high_beta")
dev.off()
