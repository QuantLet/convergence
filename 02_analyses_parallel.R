library(maps)
library(sp) 
library(spdep) 
library(maptools) 
# library(gpclib)
# install.packages("RANN")
library(RANN)
library(lmtest)
library(stringr)
library(dummies)
library(openxlsx)
library(ggplot2)
library(directlabels)
library(fBasics)
library(reshape2)

library(dplyr)
library(xts)
library(readr)
library(sf)
library(tidyr)
library(tibble)

library(RColorBrewer)
library(broom)
library(gridExtra)
library(grid)
library(ggpubr)
library(kSamples)
library(viridis)
library(Hmisc)
library(markovchain)
library(mc2d)


# functions

source("functions/function_compare_two_distributions.R")
source("functions/convergence_functions.R")
source("functions/Wilcox2003_functions.R")

################################################################################
# read data

data_income <- read_csv("data/data_income.csv")

data_exams <- read_csv("data/data_exams.csv")


##############################################################################

# prepare for transition matrices

#--- 3-yearly transitions
# edu       2003-2006, 2006-2009, 2009-2012
# vs income 2009-2012, 2012-2015, 2015-2018

# income 

data_matrix_income3y <-
  data.frame(data_income %>% 
               select(code, 
                      ends_with("2009"),
                      ends_with("2012"),
                      ends_with("2015")) %>% 
               pivot_longer(names_to = "year",
                            values_to = "y_1", 
                            -code) %>% 
               select(-year),
             data_income %>% 
               select(code, 
                      ends_with("2012"),
                      ends_with("2015"),
                      ends_with("2018")) %>% 
               pivot_longer(names_to = "year",
                            values_to = "y", 
                            -code) %>% 
               select(-year, -code))


# exams

data_matrix_exams3y <-
  data.frame(data_exams %>% 
               select(code, 
                      ends_with("2003"),
                      ends_with("2006"),
                      ends_with("2009")) %>% 
               gather(key = "year", 
                      value = "y_1", 
                      -code),
             data_exams %>% 
               select(code, 
                      ends_with("2006"),
                      ends_with("2009"),
                      ends_with("2012")) %>% 
               gather(key = "year_1", 
                      value = "y", -code)
  )

###########################################################################
# matrix 3y

#--------------------------------
# income

(matrix_income3y <- calculate_trans_matrix2(data_matrix_income3y,
                                            ngroups = 5, 
                                            lang = "EN"))

plot_income3y_podreg_matrix <-
  plot(matrix_income3y, "dodgerblue", FALSE, lang = "EN")

plot_income3y_podreg_erg <- 
  plot_erg.tmatrix(matrix_income3y, 
                   use_initial = TRUE,
                   lang = "EN")

grid.arrange(plot_income3y_podreg_matrix, 
             plot_income3y_podreg_erg, 
             ncol = 1, 
             nrow = 2)


#--------------------------------
# exams

(matrix_exams3y <- calculate_trans_matrix2(data_matrix_exams3y,
                                           ngroups = 5, 
                                           lang = "EN"))

plot_exams3y_powiaty_matrix <-
  plot(matrix_exams3y, "dodgerblue", FALSE, lang = "EN")

plot_exams3y_powiaty_erg <- 
  plot_erg.tmatrix(matrix_exams3y, 
                   use_initial = TRUE,
                   lang = "EN")

grid.arrange(plot_exams3y_powiaty_matrix, 
             plot_exams3y_powiaty_erg, 
             ncol = 1, 
             nrow = 2)

#------------------------------------------
# test of parallelism - comparison of matrices

# 3-letnie
compare_tmat3y_podreg <- test_trans_matrices2(matrix_exams3y,
                                              matrix_income3y)

compare_tmat3y_podreg_df <- data.frame(statistic = c(compare_tmat3y_podreg$chi2_stat1, 
                                                     compare_tmat3y_podreg$chi2_stat2),
                               p.value = c(compare_tmat3y_podreg$chi2_pvalue1,  
                                           compare_tmat3y_podreg$chi2_pvalue))

rownames(compare_tmat3y_podreg_df) <-  c("matrix income = matrix exams",
                                        "matrix exams = matrix income")

compare_tmat3y_podreg_df <- rownames_to_column(compare_tmat3y_podreg_df, "variant")

compare_tmat3y_podreg_df


#--------------------------------------
# test of equality of ergodic vectors

#------------------------------------------------------------------
# tests M and B from wilcox (2013) article

compare_ergvec3y <- test_ergodic_vectors_Wilcox2003(matrix_exams3y,
                                                    matrix_income3y, 
                                                    method = "both")

compare_ergvec3y


#-----------------------------------------------------------
### kernels

kernel_income3y <- calculate_cond_kde_adaptive(data_matrix_income3y)

kernel_exams3y <- calculate_cond_kde_adaptive(data_matrix_exams3y)

plot_kernel(kernel.data = kernel_income3y, 
            xlab = paste0("Relative income in year t-3"),
            ylab = paste0("Relative income in year t"),
            lang = "EN",
            cex.main = 2,
            main = "a) relative income, 2009-2018",
            gmin = 50, gmax = 350, gby = 50,
            BW = TRUE)

plot_kernel(kernel.data = kernel_exams3y, 
            xlab = paste0("Relative exam result in year t-3"),
            ylab = paste0("Relative exam result in year t"),
            lang = "EN",
            cex.main = 2,
            main = "b) relative exam result, 2003-2012",
            gmin = 80, gmax = 130, gby = 10,
            BW = TRUE)


#---------------------------------------------------------------------
# test of equality of initial distributions (after standardization)

matrix_exams3y_scaled <- scale(data_matrix_exams3y[, c("y_1", "y")])
matrix_income3y_scaled <- scale(data_matrix_income3y[, c("y_1", "y")])

data_stdy <- rbind(data.frame(variable = "exams", 
                              matrix_exams3y_scaled),
                   data.frame(variable = "income", 
                              matrix_income3y_scaled)
                   )


ggplot(data = data_stdy) + 
  stat_ecdf(aes(x = y_1,
                color = variable),
            size = 2,
            geom = "step") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.text=element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16)) +
  #scale_color_manual(values = c("#4DAF4A", "#E41A1C")) +
  scale_color_grey() +
  labs(x = "standardized value of the variable",
       y = "cumulative distribution function",
       color = element_blank())



# tests 
ks.test(x = data.frame(matrix_exams3y_scaled)$y_1,
        y = data.frame(matrix_income3y_scaled)$y_1)

ad.test(x = data.frame(matrix_exams3y_scaled)$y_1,
        y = data.frame(matrix_income3y_scaled)$y_1)


#-----------------------------------------------------------------
# comparison of ergodic - scaled analogously
# as respective initial !!!!!


# income
ergodic_income3y <- ergodicKDE(kernel_income3y)
ergodic_income3y$value <- as.integer(ergodic_income3y$value)

# rescaling according to initial income
mean_ <- as.numeric(attr(matrix_income3y_scaled, 
                           "scaled:center")[1])
sd_ <- as.numeric(attr(matrix_income3y_scaled,
                        "scaled:scale")[1])

ergodic_income3y$value <- 
  (ergodic_income3y$value - mean_)/sd_

# change densities to frequencies
ergodic_income3y$ergodic <-
  round(ergodic_income3y$ergodic * 
          (nrow(data_matrix_income3y) /
             sum(ergodic_income3y$ergodic, na.rm = TRUE)))


# exams
ergodic_exams3y <- ergodicKDE(kernel_exams3y)
ergodic_exams3y$value <- as.integer(ergodic_exams3y$value)

# przeskalowuję tak jak rozkład początkowy exams
mean_ <- as.numeric(attr(matrix_exams3y_scaled, 
                           "scaled:center")[1])
sd_ <- as.numeric(attr(matrix_exams3y_scaled,
                        "scaled:scale")[1])

ergodic_exams3y$value <- 
  (ergodic_exams3y$value - mean_)/sd_

# change densities to frequencies
ergodic_exams3y$ergodic <-
  round(ergodic_exams3y$ergodic * 
          (nrow(data_matrix_exams3y) /
             sum(ergodic_exams3y$ergodic, na.rm = TRUE)))

#-----------------------------------------------------------
# tests for ergodic

# powiaty
ks.test(x = rep(ergodic_income3y$value,
                ergodic_income3y$ergodic),
        y = rep(ergodic_exams3y$value, 
                ergodic_exams3y$ergodic))
# ad
kSamples::ad.test(x = rep(ergodic_income3y$value,
                          ergodic_income3y$ergodic),
                  y = rep(ergodic_exams3y$value, 
                          ergodic_exams3y$ergodic))




# comparison or ergodic densities

all_dens <- rbind(data.frame(variable = "income",
                             ergodic_income3y),
                  data.frame(variable = "exams",
                             ergodic_exams3y))


# empirical CDFs

ecdf_ergodic_income3y <- Hmisc::Ecdf(x = ergodic_income3y$value, 
                                     weights = ergodic_income3y$ergodic)

ecdf_ergodic_exams3y <- Hmisc::Ecdf(x = ergodic_exams3y$value, 
                                    weights = ergodic_exams3y$ergodic)

all_ecdf <- rbind(data.frame(variable = "exams",
                             x = ecdf_ergodic_exams3y$x,
                             y = ecdf_ergodic_exams3y$y),
                  data.frame(variable = "income",
                             x = ecdf_ergodic_income3y$x,
                             y = ecdf_ergodic_income3y$y)
                  )


ggplot(data = all_ecdf %>% 
             arrange(variable)) +
    geom_line(aes(x = x, y = y, 
                  group = variable,
                  col = variable),
              size = 2) +
    theme_bw() +
    theme(legend.position = "bottom", 
          legend.text=element_text(size = 20),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16)) +
    #scale_color_manual(values = c("#4DAF4A", "#E41A1C")) +
    scale_color_grey() +
    labs(x = "standardized value of the variable", 
         y = "cumulative distribution function", 
         color = element_blank()) 