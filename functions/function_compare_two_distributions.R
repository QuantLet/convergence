compare_two_distributions <- function(input_data1,
                                      input_data2,
                                      cond = "y_1", 
                                      var = "y",
                                      ngroups = 5,
                                      lang = "EN") {
  if (lang == "EN") {
    group_ = "group "
  }
  if (lang == "PL") {
    group_ = "grupa "
  }
  
  n_units = nrow(input_data1)
  data1_y_1 <- input_data1[, cond]
  data1_y <- input_data1[, var]
  data2_y_1 <- input_data2[, cond]
  data2_y <- input_data2[, var]
  
  quantiles_ <- seq(0, 1, length.out = (ngroups + 1))
  
  borders_1 <- as.matrix(quantile(data1_y_1, quantiles_ , na.rm = TRUE))
  borders_2 <- as.matrix(quantile(data2_y_1, quantiles_ , na.rm = TRUE))
  
  borders_1[1] <- 0
  borders_2[1] <- 0
  borders_1[length(borders_1)] <- 1000
  borders_2[length(borders_2)] <- 1000
  
  borders_labels1 <- rep(" ", ngroups)
  borders_labels2 <- rep(" ", ngroups)
  borders_labels1[1] <- paste0("<=", round(borders_1[2], 1),"%")
  borders_labels2[1] <- paste0("<=", round(borders_2[2], 1),"%") 
  
  for(i in 2:(length(borders_labels1) - 1)) {
    borders_labels1[i] <- paste0("(", round(borders_1[i],1),"%, ",
                                round(borders_1[i+1],1),"%]")
    borders_labels2[i] <- paste0("(", round(borders_2[i],1),"%, ",
                                 round(borders_2[i+1],1),"%]") 
  }
  borders_labels1[length(borders_labels1)] <- 
    paste0(">", round(borders_1[length(borders_labels1)], 1),"%")
  
  borders_labels2[length(borders_labels2)] <- 
    paste0(">", round(borders_2[length(borders_labels2)], 1),"%") 
  
  data1_g_1 <- cut(data1_y_1, 
             breaks = as.numeric(borders_1),
             labels = paste0(group_, 1:(length(borders_1) - 1)),
             right = TRUE)
  
  data2_g_1 <- cut(data2_y_1, 
                   breaks = as.numeric(borders_2),
                   labels = paste0(group_, 1:(length(borders_2) - 1)),
                   right = TRUE)
  
  data1_g <- cut(data1_y, 
           breaks = as.numeric(borders_1),
           labels = paste0(group_, 1:(length(borders_1) - 1)),
           right = TRUE)
  
  data2_g <- cut(data2_y, 
                 breaks = as.numeric(borders_2),
                 labels = paste0(group_, 1:(length(borders_2) - 1)),
                 right = TRUE)
  
  
  table_g_1 <- round(100 * table(data1_g_1, data2_g_1) / n_units, 2)
  table_g <- round(100 * table(data1_g, data2_g) / n_units, 2)
  
  table_diff <- table_g - table_g_1
  
  # print(table_g_1)
  # print(table_g)
  
  return(list(borders_labels1 = borders_labels1,
              borders_labels2 = borders_labels2,
              table_g_1 = table_g_1, 
              table_g = table_g,
              table_diff = table_diff))
  
  
}
