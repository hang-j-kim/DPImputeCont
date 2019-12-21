# DPImputeCont
- Imputation engine for continuous data using a Dirichlet Process (DP) Gaussian-mixture model, introduced in Kim, H. J., Reiter, J. P., Wang, Q., Cox, L. H. and Karr, A. F. (2014), "Multiple imputation of missing or faulty values under linear constraints," _Journal of Business and Economics Statistics_, 32, 375-386"

- More documentation will be added later.
- Provided that `devtools` package has been installed on the system, one can install and load the library by

```
devtools::install_github("hangjun25/DPImputeCont")
library(DPImputeCont)
```

## Example

Implementing the package while drawing some diagnostic plots with the simulated true data : Having detailed comments which explain some options for the package. Simulation datasets used are `01_Sim_Population.RData` and `02_Sim_Sample.RData`.

```
unlink(".RData")
rm(list = ls(all = TRUE))

############# Fuctions to draw scatter plot #####################

one_panel_fn <- function(i_var, j_var, Input, Output, MAIN, XLIM = NULL, YLIM = NULL) {
    plot(Input[, i_var], Input[, j_var], main = MAIN, xlab = paste("Var", i_var), ylab = paste("Var", 
        j_var), col = "red", pch = 1, cex = 0.4, xlim = XLIM, ylim = YLIM)
    points(Output[, i_var], Output[, j_var], col = "blue", pch = 20, cex = 0.4)
}  # one_panel_fn 

scatter_plot_fn <- function(D_true, Y_imputed, File_Name, i_rep) {
    Plot_Output <- Y_imputed
    if (i_rep == 0) {
        MAIN <- "Blue: Imputed data w/o chain convergence"
    } else {
        MAIN <- paste0("Blue: Imputed values ", i_rep)
    }  # 
    png(file = File_Name, width = 800, height = 700, pointsize = 20)
    par(mfrow = c(2, 2), mai = c(0.8, 0.8, 0.4, 0.4), family = "serif", mgp = c(1.5, 0.5, 0))
    one_panel_fn(i_var = 1, j_var = 2, Input = D_true, Output = Plot_Output, MAIN = "Red: True values", 
        XLIM = c(-4, 13), YLIM = c(0, 14))
    one_panel_fn(i_var = 1, j_var = 3, Input = D_true, Output = Plot_Output, MAIN = MAIN, XLIM = c(-4, 
        13), YLIM = c(-2, 14))
    one_panel_fn(i_var = 2, j_var = 3, Input = D_true, Output = Plot_Output, MAIN = "", XLIM = c(0, 
        14), YLIM = c(-2, 14))
    one_panel_fn(i_var = 2, j_var = 4, Input = D_true, Output = Plot_Output, MAIN = "", XLIM = c(0, 
        14), YLIM = c(-3, 10))
    dev.off()
}  # scatter_plot_fn

############################################################################################################ 

library(DPImputeCont)

sessionInfo()

################## Input data formats

load("01_Sim_Population.RData")
load("02_Sim_Sample.Rdata")

varnames <- dimnames(D_sample)[[2]]

################## Read the data

data_obj <- readData(Y_in = D_sample, RandomSeed = 99)

################## Make and initialize the model

model_obj <- createModel(data_obj, K_mix_comp = 30)
# K_mix_comp is the number of the mixture components.  The default value is 50.  Empirically, K=30
# is enough for most datasets with the number of variables less than 5.  After running the code, we
# can check if K_mix_comp is enough big by summarizing result_obj$draw_no_occ_cluster

Y_imputed <- Y_std <- model_obj$Y_mat
for (i_sample in 1:dim(Y_imputed)[[1]]) {
    Y_imputed[i_sample, ] <- data_obj$mean_Y_input + data_obj$sd_Y_input * Y_std[i_sample, ]
}
scatter_plot_fn(D_true = D_pop, Y_imputed = Y_imputed, File_Name = "12_MI_ScatterPlots_0.png", i_rep = 0)

# run 1 iteration.
model_obj$Iterate()
# model_obj$.where_we_are

# run many iterations
burn <- 500
m_Imp <- 10
thin <- 100
# burn=10; m_Imp=5; thin=10
where_sample_EI <- burn + seq(from = thin, to = (m_Imp * thin), by = thin)

result_obj <- multipleImp(model_obj = model_obj, n_burnin = burn, m_Imp = m_Imp, interval_btw_Imp = thin)

################## Save results and make plots

save(data_obj = data_obj, model_obj = model_obj, result_obj = result_obj, varnames = varnames, where_sample_EI = where_sample_EI, 
    burn = burn, file = "11a_ImputedData.RData")

for (i_Imp in 1:m_Imp) {
    scatter_plot_fn(D_true = D_pop, Y_imputed = result_obj$multiple_Imp[i_Imp, , ], File_Name = paste0("12_MI_ScatterPlots_", 
        i_Imp, ".png"), i_rep = i_Imp)
}  # for (i_EI in 1:m) 

png(file = "12_Imp_ConvCheck.png", width = 800, height = 700, pointsize = 20)

par(mfrow = c(2, 2), mai = c(0.8, 0.8, 0.4, 0.4), family = "serif", mgp = c(1.5, 0.5, 0))  # b l t r

i_var <- 1
plot(result_obj$draw_weighted_mu[, i_var], main = paste0("Weighted mu for variable ", i_var), type = "l", 
    ylab = "Weighted mu", xlab = "Iteration")
# abline(h = mean(result_obj$draw_weighted_mu[,i_var]), lty='dotted', col='green', lwd=2)
abline(h = True_Weighted_Mu[i_var], lty = "dotted", col = "green", lwd = 2)
abline(v = burn, col = "red", lwd = 2, lty = "dotted")
abline(v = where_sample_EI, col = "blue", lty = "dotted", lwd = 1.2)

i_var <- 2
plot(result_obj$draw_weighted_mu[, i_var], main = paste0("Weighted mu for variable ", i_var), type = "l", 
    ylab = "Weighted mu", xlab = "Iteration")
# abline(h = mean(result_obj$draw_weighted_mu[,i_var]), lty='dotted', col='green', lwd=2)
abline(h = True_Weighted_Mu[i_var], lty = "dotted", col = "green", lwd = 2)
abline(v = burn, col = "red", lwd = 2, lty = "dotted")
abline(v = where_sample_EI, col = "blue", lty = "dotted", lwd = 1.2)

plot(result_obj$draw_no_occ_cluster, main = "No. of occupying clusters", type = "l", ylab = "No. of clusters", 
    xlab = "Iteration", ylim = c(-0.5, 20))
# abline(h = mean(result_obj$draw_no_occ_cluster), lty='dotted', col='green', lwd=2)
abline(v = burn, col = "red", lwd = 2, lty = "dotted")
abline(v = where_sample_EI, col = "blue", lty = "dotted", lwd = 1.2)

plot(result_obj$draw_alpha, main = "alpha - Concentration param.", type = "l", ylab = "alpha", xlab = "Iteration")
# abline(h = mean(result_obj$draw_alpha), lty='dotted', col='green', lwd=2)
abline(v = burn, col = "red", lwd = 2, lty = "dotted")
abline(v = where_sample_EI, col = "blue", lty = "dotted", lwd = 1.2)

dev.off()
```


## Implementing the package (used in practice where we do not know the true data)

```
library(DPImputeCont)

sessionInfo()

####################### 

load("02_Sim_Sample.Rdata")
head(D_sample)
varnames <- dimnames(D_sample)[[2]]

data_obj <- readData(Y_in = D_sample, RandomSeed = 99)

model_obj <- createModel(data_obj, K_mix_comp = 30)

burn <- 500
m_Imp <- 10
thin <- 100
result_obj <- multipleImp(model_obj = model_obj, n_burnin = burn, m_Imp = m_Imp, interval_btw_Imp = thin)

save(data_obj = data_obj, model_obj = model_obj, result_obj = result_obj, varnames = varnames, burn = burn, 
    file = "11b_ImputedData.RData")

####################### 

dim(result_obj$multiple_Imp)

head(result_obj$multiple_Imp[1, , ])

head(result_obj$multiple_Imp[m_Imp, , ])
```
