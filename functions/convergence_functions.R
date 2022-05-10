#------------------------------------------------------------
# univariate Gausian Kernel
# input:
# x - numeric vector
# output:
# value of Gaussian in point x

kernel_univ <- function(x) { 1 / sqrt(2 * pi) * exp(-0.5 * x**2)}

# generalize to use different kernel functions !!!
# package kedd !!!

#------------------------------------------------------------
# bivariate kernel
# Input:
# x1   n1xm1-matrix, containing first variable points at which kernel function
#	is to be evaluated.
# x2   n2xm2-matrix, conformable with x1 containing second variable points at
#	which the kernel function is to be evaluated.
# Output:
#     nxm-matrix, conformable with x1 and x2 with value of kernel function at
#	implied points. 

kernel_biv <- function(x1, x2) { kernel_univ(x1) * kernel_univ(x2) }

#------------------------------------------------------------
# joint distribution kernel

# input 
# y - numeric vector of observations of variable y
# x - numeric vector of observations of variable x
# hwy - a single value or vector of length nrow(y) of optimal bandwidth for y points
# hwx - a single value or vector of length nrow(x) of optimal bandwidth for x points
# xgrid - vector of gridpoints on x dimension for which f(y, x) will be calculated
# xgrid - vector of gridpoints on y dimension for which f(y, x) will be calculated
# output:  
# matrix_ L x m zawierajaca wartosci f(y,x) dla punktow (x,y) z xgrid, ygrid

# zmieniæ:
# - dodaj automatyczne obliczanie optimal bandwidth 
# niech xgrid, ygrid nie bêdzie matrix_¹

kernel_joint <- function (y, x, hwy, hwx, ygrid, xgrid) {
  L <- ifelse(is.matrix(ygrid), nrow(ygrid), length(ygrid))
  m <- ifelse(is.matrix(xgrid), ncol(xgrid), length(xgrid))
  if (length(hwy) == 1) hwy <- rep(hwy, length(y))
  if (length(hwx) == 1) hwx <- rep(hwx, length(x))
  fjy_x <- matrix(0, L, m)
  y1 <- matrix(rep(y, L), ncol = L, byrow = F)
  y_grid1 <- matrix(rep(ygrid, nrow(y1)), ncol = L, byrow = T)
  y_y <- y_grid1 - y1
  x1 <- matrix(rep(x, m), ncol = m, byrow = F)
  x_grid1 <- matrix(rep(xgrid, nrow(x1)), ncol = L, byrow = T)
  x_x <- x_grid1 - x1
  fjy_x <- (t(kernel_univ(y_y / hwy) / hwy) %*% 
              (kernel_univ(x_x / hwx) / hwx) ) / L
  return (fjy_x)
}

#------------------------------------------------------------
# conditional kernel distribution
# input 
# y - numeric vector of observations of variable y
# x - numeric vector of observations of variable x
# hwx - optimal bandwidth for x (a single value or vector of length nrow(x))
# hwy - optimal bandwidth for y (a single value or vector of length nrow(y))
# xgrid - vector of gridpoints on x dimension for which f(y|x) will be calculated
# xgrid - vector of gridpoints on y dimension for which f(y|x) will be calculated
# Output
# matrix_ L x m zawierajaca wartosci f(y|x) dla punktow (x,y) z xgrid, ygrid

kernel_cond <- function (y, x, hwy, hwx, ygrid, xgrid) {
  L <- ifelse(is.matrix(ygrid), nrow(ygrid), length(ygrid))
  m <- ifelse(is.matrix(xgrid), ncol(xgrid), length(xgrid))
  if (length(hwy) == 1) hwy <- array(hwy, length(y))
  if (length(hwx) == 1) hwx <- array(hwx, length(x))
  
  fjy_x <- matrix(0, L, m)
  
  y1 <- matrix(rep(y, L), ncol = L, byrow = F)
  y_grid1 <- matrix(rep(ygrid, nrow(y1)), ncol = L, byrow = T)
  y_y <- y_grid1 - y1
  x1 <- matrix(rep(x, m), ncol = m, byrow = F)
  x_grid1 <- matrix(rep(xgrid, nrow(x1)), ncol = L, byrow = T)
  x_x <- x_grid1 - x1
  fjy_x <-  
    t( (t(kernel_univ(y_y / as.vector(hwy)) / as.vector(hwy)) %*% 
          (kernel_univ(x_x / as.vector(hwx)) / as.vector(hwx)) ) )/
    colSums( kernel_univ(x_x / as.vector(hwx)) / as.vector(hwx) )
  
  fjy_x <- t(fjy_x)
  
  return (fjy_x)
}

optimal_bandwidth <- function(x) {
  1.06 * min(sd(x, na.rm = T), IQR(x, na.rm = T) / 1.34) * length(x)**(-0.2)
}
  

#------------------------------------------------------------
# conditional distribution kernel - simple one step estimation
# input:
# input.data - data.frame including a vector of a considered variable (y)
#              and a conditioning variable (x) for f(y|x)
# var - the name of numeric column with values of the considered variable
# cond - the name of numeric column with values of the conditioning variable

# zmien
# - dodaj argument ngrids

calculate_cond_kde <- function (input.data, cond = "y_1", var = "y",
                                  grid_by = 0.05, grid_min = NA,
								  grid_max = NA) {
  input.data_ <- input.data[, c(cond, var)]
  
  scale_ <- max(abs(min(input.data_)), 
                abs(max(input.data_))
                )
  scale_ <- 10**(nchar(ceiling(scale_)) - 1)
  
  input.data_ <- input.data_ / scale_
  
  if(is.na(grid_min)) grid_min <- floor(10 * min(input.data_)) / 10
  if(is.na(grid_max))  grid_max <- ceiling(10 * max(input.data_)) / 10
  
  # jesli parametr ngrids pusty, zrob co 0.05
  ngrids <- round(((grid_max - grid_min) / grid_by + 1), 0)
  
  n_ <- nrow(input.data_)
  # initially all weights = 1
  w <- array(1, n_) 
  y <- input.data_[, var]
  y_1 <- input.data_[, cond]
  
  kernel.grids <- seq(grid_min, grid_max, by = grid_by)
  
  # optimal bandwidth
  # Silverman's rule od thumb - uogólnij !!!!
  # to bardziej elegancko zapisz
  data.opthy <- optimal_bandwidth(y)
  data.opthx <- optimal_bandwidth(y_1)
  
  fy_x1 <- kernel_cond(y, y_1, data.opthy, data.opthx, kernel.grids, kernel.grids)
  
  kernel.c1 <- rep(kernel.grids, ngrids)
  kernel.c2 <- sort(kernel.c1)
  kernel.c3 <- as.numeric(fy_x1)
  
  kernel.out <- data.frame(y = kernel.c1,
                           x = kernel.c2,
                           density = kernel.c3)
  
  return(kernel.out)
}

#------------------------------------------------------------
# conditional distribution kernel - two step adaptive estimation
# input:
# input.data - data.frame including a vector of a considered variable (y)
#              and a conditioning variable (x) for f(y|x)
# var - the name of numeric column with values of the considered variable
# cond - the name of numeric column with values of the conditioning variable

# zmien
# - dodaj argument ngrids

calculate_cond_kde_adaptive <- function (input.data, cond = "y_1", var = "y",
                                  grid_by = 0.05, grid_min = NA, grid_max = NA) {
  input.data_ <- input.data[, c(cond, var)]
  
  scale_ <- max(abs(min(input.data_)), 
                abs(max(input.data_))
                )
  scale_ <- 10**(nchar(ceiling(scale_)) - 1)
  
  input.data_ <- input.data_ / scale_
  
	if(is.na(grid_min)) grid_min <- floor(10 * min(input.data_)) / 10
	if(is.na(grid_max)) grid_max <- ceiling(10 * max(input.data_)) / 10
	
	# jesli parametr ngrids pusty, zrob co 0.05
	ngrids <- round(((grid_max - grid_min) / grid_by + 1), 0)
	
	n_ <- nrow(input.data_)
	# initially all weights = 1
	w <- array(1, n_) 
	y <- input.data_[, var]
	y_1 <- input.data_[, cond]
	
	kernel.grids <- seq(grid_min, grid_max, by = grid_by)

	# optimal bandwidth
	# Silverman's rule od thumb - uogólnij !!!!
	# to bardziej elegancko zapisz
	data.opthy <- optimal_bandwidth(y)
	data.opthx <- optimal_bandwidth(y_1)
	
	fy_x2 <- kernel_joint(y = y, x = y_1, 
	                      hwy = data.opthy, hwx = data.opthx,
	                      ygrid = y, xgrid = y_1)

	fy_x3 <- diag(fy_x2) # tworzy wektor z elementow na przekatnej 
	
	fg <- exp(mean(log(fy_x3))) # geometric mean of weights
	w2 <- w * ((fy_x3 / fg)**(-0.5)) # adaptive weights
	data.opthy <- w2 * data.opthy
	data.opthx <- w2 * data.opthx

	fy_x1 <- kernel_cond(y, y_1, data.opthy, data.opthx, kernel.grids, kernel.grids)

	kernel.c1 <- rep(kernel.grids, ngrids)
	kernel.c2 <- sort(kernel.c1)
	kernel.c3 <- as.numeric(fy_x1)
	
	kernel.out <- data.frame(y = kernel.c1,
	                         x = kernel.c2,
	                         density = kernel.c3)

	return(kernel.out)
}

#------------------------------------------------------------
# calculating transition matrix
# input:
# input.data - data.frame including a vector of a considered variable (y)
#              and a conditioning variable (x) for f(y|x)
# var - the name of numeric column with values of the considered variable
# cond - the name of numeric column with values of the condition

calculate_trans_matrix <- function(input.data,
                                   start = "y_1", end = "y",
                                   ngroups = 5,
								   lang = "EN") {
  require(markovchain)
  
  if (lang == "EN") {
       group_ = "group "
	   ergodic_ = "ergodic"
	   }
  if (lang == "PL") {
       group_ = "grupa "
	   ergodic_ = "ergodyczny"
	   }
  
  y_1 <- input.data[, start]
  y <- input.data[, end]
  quantiles_ <- seq(0, 1, length.out = (ngroups + 1))
  borders_ <- as.matrix(quantile(y_1, quantiles_ , na.rm = TRUE))
  
  borders_[1] <- 0
  borders_[length(borders_)] <- 1000
  
  borders.labels <- rep(" ", ngroups)
  borders.labels[1] <- paste0("<=", round(borders_[2], 1),"%") 
  
  for(i in 2:(length(borders.labels) - 1)) {
    borders.labels[i] <- paste0("(", round(borders_[i],1),"%, ",
                                round(borders_[i+1],1),"%]") 
    }
  borders.labels[length(borders.labels)] <- 
    paste0(">", round(borders_[length(borders.labels)], 1),"%") 
  
  g_1 <- cut(y_1, 
             breaks = as.numeric(borders_),
             labels = paste0(group_, 1:(length(borders_) - 1)),
             right = TRUE)
  
  g <- cut(y, 
           breaks = as.numeric(borders_),
           labels = paste0(group_, 1:(length(borders_) - 1)),
           right = TRUE)
  
  initial_ <- round(100 * table(g_1) / sum(table(g_1)), 1)
  
  trans.matrix <- as.matrix(table(g_1, g))
  
  if(ncol(trans.matrix) < nrow(trans.matrix)) {
    cols_m_ <- matrix(0, nrow = nrow(trans.matrix),
                      ncol = nrow(trans.matrix) - ncol(trans.matrix))
    colnames(cols_m_) <- paste0(group_, 
                                c(1:nrow(trans.matrix)))[(1 - as.numeric(1:ngroups %in% gsub(group_, "", colnames(trans.matrix)))) == 1]
    trans.matrix <- cbind(trans.matrix,cols_m_)
    trans.matrix <- trans.matrix[, sort(colnames(trans.matrix))]
    }
  
  if(ncol(trans.matrix) > nrow(trans.matrix)) {
    rows_m_ <- matrix(0, nrow = ncol(trans.matrix) - nrow(trans.matrix),
                      ncol = ncol(trans.matrix))
    rownames(rows_m_) <- paste0(group_, 
                                c(1:ncol(trans.matrix)))[(1 - as.numeric(1:ngroups %in% gsub(group_, "", rownames(trans.matrix)))) == 1]
    trans.matrix <- rbind(trans.matrix, rows_m_)
    trans.matrix <- trans.matrix[sort(rownames(trans.matrix)),]
    }
  
  for (nazwa_col in paste0(group_, 1:ngroups)) {
    if(! nazwa_col %in% colnames(trans.matrix)) {
      trans.matrix <- cbind(trans.matrix, rep(0, nrow(trans.matrix)))
      colnames(trans.matrix)[ncol(trans.matrix)] <- nazwa_col
      trans.matrix <- rbind(trans.matrix, rep(0, ncol(trans.matrix)))
      rownames(trans.matrix)[nrow(trans.matrix)] <- nazwa_col
    }
    }
  trans.matrix <- trans.matrix[, sort(colnames(trans.matrix))]
  trans.matrix <- trans.matrix[sort(rownames(trans.matrix)), ]
  
  count <- as.numeric(rowSums(trans.matrix, na.rm = T))
  
  trans.matrix <- 100 * trans.matrix / replace(count, which(count == 0), 1)
  
  rownames(trans.matrix) <- paste0(rownames(trans.matrix),
                                   " (", count, ")")
  colnames(trans.matrix) <- paste0(colnames(trans.matrix),
                                        " \n ", borders.labels)
  
  # adding an ergodic vector
  
  trans.matrix.MC <- new("markovchain", 
                         transitionMatrix = matrix(trans.matrix/100, 
                                                   nrow = ngroups, 
                                                   byrow = F)
                         )
  erg_ <- round(100 * steadyStates(trans.matrix.MC), 2)
  
  trans.matrix = rbind(trans.matrix, erg_)
  
  rownames(trans.matrix)[ngroups + 1] <- ergodic_
  
  round(trans.matrix, 1)
}

#------------------------------------------------------------
# calculating transition matrix
# input:
# input.data - data.frame including a vector of a considered variable (y)
#              and a conditioning variable (x) for f(y|x)
# var - the name of numeric column with values of the considered variable
# cond - the name of numeric column with values of the condition
# ngroups - number of groups in a transition matrix
# lang - language of labels ("EN" - default or "PL")

calculate_trans_matrix2 <- function(input.data,
                                    cond = "y_1", var = "y",
                                    ngroups = 5,
                                    lang = "EN") {
  require(markovchain)
  
  if (lang == "EN") {
    group_ = "group "
    ergodic_ = "ergodic"
  }
  if (lang == "PL") {
    group_ = "grupa "
    ergodic_ = "ergodyczny"
  }
  
  n_units = nrow(input.data)
  y_1 <- input.data[, cond]
  y <- input.data[, var]
  quantiles_ <- seq(0, 1, length.out = (ngroups + 1))
  borders_ <- as.matrix(quantile(y_1, quantiles_ , na.rm = TRUE))
  
  borders_[1] <- 0
  borders_[length(borders_)] <- 1000
  
  borders.labels <- rep(" ", ngroups)
  borders.labels[1] <- paste0("<=", round(borders_[2], 1),"%") 
  
  for(i in 2:(length(borders.labels) - 1)) {
    borders.labels[i] <- paste0("(", round(borders_[i],1),"%, ",
                                round(borders_[i+1],1),"%]") 
  }
  borders.labels[length(borders.labels)] <- 
    paste0(">", round(borders_[length(borders.labels)], 1),"%") 
  
  g_1 <- cut(y_1, 
             breaks = as.numeric(borders_),
             labels = paste0(group_, 1:(length(borders_) - 1)),
             right = TRUE)
  
  g <- cut(y, 
           breaks = as.numeric(borders_),
           labels = paste0(group_, 1:(length(borders_) - 1)),
           right = TRUE)
  
  initial_ <- round(100 * table(g_1) / sum(table(g_1)), 1)
  
  trans.matrix <- as.matrix(table(g_1, g))
  
  if(ncol(trans.matrix) < nrow(trans.matrix)) {
    cols_m_ <- matrix(0, nrow = nrow(trans.matrix),
                      ncol = nrow(trans.matrix) - ncol(trans.matrix))
    colnames(cols_m_) <- paste0(group_, 
                                c(1:nrow(trans.matrix)))[(1 - as.numeric(1:ngroups %in% gsub(group_, "", colnames(trans.matrix)))) == 1]
    trans.matrix <- cbind(trans.matrix,cols_m_)
    trans.matrix <- trans.matrix[, sort(colnames(trans.matrix))]
  }
  
  if(ncol(trans.matrix) > nrow(trans.matrix)) {
    rows_m_ <- matrix(0, nrow = ncol(trans.matrix) - nrow(trans.matrix),
                      ncol = ncol(trans.matrix))
    rownames(rows_m_) <- paste0(group_, 
                                c(1:ncol(trans.matrix)))[(1 - as.numeric(1:ngroups %in% gsub(group_, "", rownames(trans.matrix)))) == 1]
    trans.matrix <- rbind(trans.matrix, rows_m_)
    trans.matrix <- trans.matrix[sort(rownames(trans.matrix)),]
  }
  
  for (nazwa_col in paste0(group_, 1:ngroups)) {
    if(! nazwa_col %in% colnames(trans.matrix)) {
      trans.matrix <- cbind(trans.matrix, rep(0, nrow(trans.matrix)))
      colnames(trans.matrix)[ncol(trans.matrix)] <- nazwa_col
      trans.matrix <- rbind(trans.matrix, rep(0, ncol(trans.matrix)))
      rownames(trans.matrix)[nrow(trans.matrix)] <- nazwa_col
    }
  }
  trans.matrix <- trans.matrix[, sort(colnames(trans.matrix))]
  trans.matrix <- trans.matrix[sort(rownames(trans.matrix)), ]
  
  initial_counts <- as.numeric(rowSums(trans.matrix, na.rm = T))
  final_counts <- as.numeric(colSums(trans.matrix, na.rm = T))
  
  trans.matrix <- 100 * trans.matrix / replace(initial_counts, which(initial_counts == 0), 1)
  
  # adding an ergodic vector
  
  trans.matrix.MC <- new("markovchain", 
                         transitionMatrix = matrix(trans.matrix/100, 
                                                   nrow = ngroups, 
                                                   byrow = F)
  )
  erg_ <- round(100 * steadyStates(trans.matrix.MC), 2)
  ergodic_counts <- round(0.01 * erg_ * n_units)
  
  if(sum(ergodic_counts) != n_units) 
    ergodic_counts[which.max(ergodic_counts)] <- ergodic_counts[which.max(ergodic_counts)] + n_units - sum(ergodic_counts)
  
  trans.matrix <- round(trans.matrix, 2)
  
  half_life = -log(2)/log(eigen(trans.matrix/100)$values[2])
  
  trans.matrix_return <- list(trans_probs = trans.matrix,
                              ergodic = erg_,
                              n_units = n_units,
                              n_groups = ngroups,
                              group_borders = borders.labels,
                              initial_counts = initial_counts,
                              final_counts = final_counts,
                              ergodic_counts = ergodic_counts,
							  average_persistence_prob = mean(diag(trans.matrix)),
							  half_life = half_life)
  
  structure(trans.matrix_return, class = "tmatrix")
}

#-----------------------------------------------------------
# metoda plot dla tmatrix  
plot.tmatrix <- function(tmatrix, 
                         gradient.high = "darkred",
                         show.zeros = T,
                         lang = "EN")
{
  require(ggplot2)
  
  if (lang == "EN") {
    ylab_ <- "group in initial period"
    xlab_ <- "group in final period"
    options(OutDec =  ".")
  }
  
  if (lang == "PL") {
    ylab_ <- "grupa w okresie pocz¹tkowym"
    xlab_ <- "grupa w okresie koñcowym"
	options(OutDec =  ",")
  }
  
  ngroups <- tmatrix$n_groups
  matrix.df <- data.frame(tmatrix$trans_probs)
  cell.labels <- paste0(sprintf("%2.1f", matrix.df$Freq), "%")
  
  if(show.zeros == F) cell.labels <- gsub("^0.0%$","", cell.labels, perl = T)
  
  if(lang == "PL") cell.labels <- gsub(".",",", cell.labels, fixed = T)
  
  plot1 <- ggplot(data =  matrix.df, 
                  aes(x = g, y = g_1)) +
    geom_tile(aes(fill = Freq), colour = "darkgray") +
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          panel.background = element_rect(fill = 'white', colour = 'white')) +
    scale_x_discrete(labels = paste0(1:ngroups, "\n", tmatrix$group_borders)) +
    scale_y_discrete(labels = paste0(1:ngroups, " (", tmatrix$initial_counts, ")")) +
    ylab(ylab_) +
    xlab(xlab_) +
    scale_fill_gradient2(low = "white", high = gradient.high, guide = F) +
    annotate("text", x = as.numeric(matrix.df$g), 
             y = as.numeric(matrix.df$g_1), 
             label = cell.labels, 
             color = "black", 
             fontface = 1,
             size = 8)
  
  return(plot1)
}
 

#------------------------------------------------------------
# zamienic na metodê plot() dla trans.matrix???

plot_trans_matrix <- function(matrix_, 
                              gradient.high = "darkred",
                              show.zeros = T,
							  lang = "EN")
{
  require(ggplot2)
  require(reshape2)
  
  if (lang == "EN") {
     ylab_ <- "group in initial period"
	 xlab_ <- "group in final period"
	 options(OutDec =  ".")
  }
  
  if (lang == "PL") {
     ylab_ <- "grupa w okresie pocz¹tkowym"
	 xlab_ <- "grupa w okresie koñcowym"
	 options(OutDec =  ",")
  }
  
  ngroups <- ncol(matrix_)
  
  matrix.df <- data.frame(group_y_1 = 1:ngroups,
                           matrix_[-(ngroups+1),], stringsAsFactors = F)
  
  row.names(matrix.df) <- NULL
  names(matrix.df)[-1] <- 1:ngroups
  
  matrix.plot <- melt(matrix.df, id = "group_y_1")
  
  names(matrix.plot)[2:3] <- c("group_y", "p")
  
  cell.labels <- paste0(sprintf("%2.1f", matrix.plot$p), "%")
  
  if(show.zeros == F) cell.labels <- gsub("^0.0%$","", cell.labels, perl = T)
  
  if(lang == "PL") cell.labels <- gsub(".",",", cell.labels, fixed = T)
  
  plot1 <- ggplot(data =  matrix.plot, 
                  aes(x = group_y, y = group_y_1)) +
    geom_tile(aes(fill = p), colour = "darkgray") +
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 16),
          panel.background = element_rect(fill = 'white', colour = 'white')) +
    scale_x_discrete(labels = colnames(matrix_)) +
    scale_y_discrete(limit = gsub(" (","\n(", rownames(matrix_)[1:ngroups], fixed = T)) +
    ylab(ylab_) +
    xlab(xlab_) +
    scale_fill_gradient2(low = "white", high = gradient.high, guide = F) +
    annotate("text", x = as.numeric(matrix.plot$group_y), 
             y = as.numeric(matrix.plot$group_y_1), 
             label = cell.labels, 
             color = "black", 
             fontface = 1,
             size = 8)
  
  return(plot1)
}



#------------------------------------------------------------
# plot_trans_matrix_erg

plot_trans_matrix_erg <- function(matrix_, lang = "EN", horiz = F) {
  
  require(stringr)
  require(scales)
  
  if (lang == "EN") {
    initial_text <- "initial"
	ergodic_text <- "ergodic"
	prob_txt <- "probability (%)"
	xlab_txt <- "group"
	x_name <- "distribution"
	options(OutDec =  ".")
  }
  
   if (lang == "PL") {
    initial_text <- "pocz¹tkowy"
	ergodic_text <- "ergodyczny"
	prob_txt <- "prawdopodobieñstwo, %"
	xlab_txt <- "grupa"
	x_name <- "rozk³ad"
	options(OutDec =  ",")
  }
  
  ngroups <- ncol(matrix_)
  initial_ <- str_extract_all(row.names(matrix_), "\\([^()]+\\)", simplify = T)[1:ngroups]
  initial_ <- as.numeric(substring(initial_, 2, nchar(initial_) - 1))
  initial_ <- round(100 * initial_/sum(initial_), 1)
  
  ergodic_ <- as.numeric(matrix_[ngroups + 1,])
  
  comparison_ <- rbind(initial_ , ergodic_)
  rownames(comparison_) <- c(initial_text, ergodic_text)
  colnames(comparison_) <- colnames(matrix_)
  
  comparison_ <- data.frame(t(comparison_))
  comparison_$labels <- rownames(comparison_)
  rownames(comparison_) <- NULL
  
  # zeby porzadek etykiet (przedzia³ow) byl wlasciwy
  levels_ <- comparison_$labels  
  comparison_ <- melt(comparison_, id = "labels")
  
  comparison_$labels <- factor(comparison_$labels, ordered = T,
                               levels = levels_)
  
  labels_ <- paste0(sprintf("%2.1f", comparison_$value),"%")
  if(lang == "PL") labels_ <- gsub(".", ",", labels_, fixed = T)
  
  p <- ggplot(comparison_, 
         aes(x = labels, y = value, fill = variable)) +  
    geom_bar(position = "dodge", 
             stat = "identity") + 
    labs(y = prob_txt, size = 2, 
         fill = x_name) + 
    scale_fill_brewer(palette = "Set1") +
    theme_bw(base_size = 20) +
    theme(legend.position = 'bottom', 
	      legend.text = element_text(size = 20),
		  axis.title.x = element_blank()) +
    geom_text(aes(y = value+1,
                  label = labels_),
              position = position_dodge(width = 0.9), 
              vjust = 0, hjust = 0.5, size = 5)
			  
   if (horiz == T) p <- p + coord_flip()
   
   p  
}


#-----------------------------------------------------------
# plot_erg <- function() UseMethod("plot_erg")

plot_erg.tmatrix <- function(tmatrix, 
                             lang = "EN", 
                             use_initial = T,
                             use_final = F,
                             use_ergodic = T,
                             scale.counts = F,
                             horiz = F) {
  if(use_initial + use_final + use_ergodic == 0) {
    warning("At least one of the dimensions (initial, final, ergodic) has to be plotted!")
    stop }
  
  require(scales)
  require(tidyr)
  
  if (lang == "EN") {
    initial_text <- "initial"
    final_text <- "final"
    ergodic_text <- "ergodic"
    ergodic1_text <- "ergodic no 1"
    ergodic2_text <- "ergodic no 2"
    prob_txt <- ifelse(scale.counts == F, "probability (%)", "count")
    xlab_txt <- "group"
    x_name <- "distribution"
    options(OutDec =  ".")
  }
  
  if (lang == "PL") {
    initial_text <- "pocz¹tkowy"
    final_text <- "koñcowy"
    ergodic_text <- "ergodyczny"
    ergodic1_text <- "ergodyczny nr 1"
    ergodic2_text <- "ergodyczny nr 2"
    prob_txt <- ifelse(scale.counts == F, "prawdopodobieñstwo (%)", "liczebnoœci")
    xlab_txt <- "grupa"
    x_name <- "rozk³ad"
    options(OutDec =  ",")
  }
  
  comparison_ <- data.frame(group = as.factor(1:tmatrix$n_groups),
                            initial = tmatrix$initial_counts, 
                            final = tmatrix$final_counts, 
                            ergodic = t(tmatrix$ergodic_counts))
  
  names(comparison_) <- gsub(".", "", names(comparison_),
                             fixed = TRUE)
  
  if(use_initial == F) comparison_$initial <- NULL
  if(use_final == F) comparison_$final <- NULL
  if(use_ergodic == F) comparison_[,grep("ergodic", 
                                        names(comparison_))] <- NULL
  
  if(use_initial + use_final + use_ergodic > 1) {
    if(scale.counts == F) comparison_[, -1] <- round(100*comparison_[, -1]/ sum(comparison_[, 2], na.rm = T), 1)
    comparison_ <- tidyr::gather(comparison_, distribution, value, -group) } else {
      if(scale.counts == F) comparison_[, -1] <- round(100*comparison_[, -1]/ sum(comparison_[, -1], na.rm = T), 1)
      comparison_$distribution <- names(comparison_)[2]
      names(comparison_)[2] <- "value"
    }
  
  # zeby porzadek etykiet (przedzia³ow) byl wlasciwy
  comparison_$distribution <- as.factor(sapply(as.list(comparison_$distribution),
                                               function(x) get(paste0(x, "_text"))))
  
  labels_ <- paste0(sprintf("%2.1f", comparison_$value),"%")
  if(lang == "PL") labels_ <- gsub(".", ",", labels_, fixed = T)
  
  p <- ggplot(comparison_, 
              aes(x = group, 
                  y = value, fill = distribution)) +  
    geom_bar(position = "dodge", 
             stat = "identity") + 
    labs(x = xlab_txt, y = prob_txt, size = 2, 
         fill = x_name) + 
    scale_fill_brewer(palette = "Set1") +
    theme_bw(base_size = 20) +
    theme(legend.position = 'bottom', legend.text = element_text(size = 20)) 
  
  if (scale.counts == F) { p <- p +
    geom_text(aes(y = value + 1,
                  label = labels_),
              position = position_dodge(width = 0.9), 
              vjust = 0, hjust = 0.5, size = 5) } else
              {p <- p +
                geom_text(aes(y = ifelse(value == min(value), value * 1.2, 0.9 * value),
                              label = value),
                          position = position_dodge(width = 0.9), 
                          vjust = 0, hjust = 0.5, size = 5)}
  
  if (horiz == T) p <- p + coord_flip()
  
  p  
}



#------------------------------------------------------------
# plot kernel

plot_kernel <- function(kernel.data, gmin, gmax, gby, xlab, ylab, 
                        main = "", cex.main = 1,
						use_palette = colorRampPalette(rev(brewer.pal(9, "Spectral"))),
                        BW = F, nlevels = 10, lang = "EN") { 
  
  require(RColorBrewer)
  require(grDevices)
  
  kernel_ <- reshape(kernel.data, idvar = "x", timevar = "y", direction = "wide")
  if(lang == "EN") main.scale = "density" else main.scale = "gêstoœæ"
  
  kernel2_ <- as.matrix(kernel_[, 2:ncol(kernel_)])
  x <- as.matrix(100 * kernel_[,1])
  
  if (BW) {
    contour(x, x, kernel2_, xlab = xlab, ylab = ylab, nlevels = nlevels,
            main = main, cex.main = cex.main, cex.lab = 2, 
			cex.axis = 2, drawlabels = F, las = 1, 
            asp = 1)
    abline(h = seq(gmin, gmax, by = gby), col = "gray", lty = "dashed", lwd = 2)
    abline(v = seq(gmin, gmax, by = gby), col = "gray", lty = "dashed", lwd = 2)
    abline(a = 0, b = 1, col = "black", lty = 4, lwd = 4)
  } else {
    filled.contour(x, x, kernel2_, 
                   xlab = xlab, 
                   ylab = ylab, 
                   key.title = title(main = main.scale),
				   main=main,
				   color.palette = use_palette, 
                 asp = 1, 
				 plot.title = title(main = main, cex.main = cex.main),
                 plot.axes = {
                   axis(1, seq(gmin, gmax, by = gby),
                        cex.axis = 1.3,
                        cex = 1.3,
                        cex.lab = 1.3);
                   axis(2, seq(gmin, gmax, by = gby),
                        cex.axis = 1.3,
                        cex = 1.3,
                        cex.lab = 1.3); 
                   abline(h = seq(gmin, gmax, by = gby),
                          col = "white", lty = "dashed", lwd = 2); 
                   abline(v = seq(gmin, gmax, by = gby),
                          col = "white", lty = "dashed", lwd = 2);
                   abline(a = 0, b = 1, col = "black", lty = 4, lwd = 4)
                 })
    }
  }


plot_kernel_cde <- function(kernel.data, gmin, gmax, gby, xlab, ylab, 
                            main = "",
                            use_palette = colorRampPalette(rev(brewer.pal(9, "Spectral"))),
                            BW = F, nlevels = 10, lang = "EN") { 
  
  require(RColorBrewer)
  require(grDevices)
  
  kernel_x <- kernel.data$x
  
  kernel_ <- matrix(kernel.data$z, nrow = length(kernel_x))
  
  if(lang == "EN") main.scale = "density" else main.scale = "gêstoœæ"
  
  if (BW) {
    contour(kernel_x, kernel_x, kernel_, xlab = xlab, ylab = ylab, nlevels = nlevels,
            main = main, cex.lab = 2, cex.axis = 2, drawlabels = F,
            asp = 1)
    abline(h = seq(gmin, gmax, by = gby), col = "gray", lty = "dashed", lwd = 2)
    abline(v = seq(gmin, gmax, by = gby), col = "gray", lty = "dashed", lwd = 2)
    abline(a = 0, b = 1, col = "black", lty = 4,lwd = 4)
  } else {
    filled.contour(kernel_x, kernel_x, kernel_, 
                   xlab = xlab, 
                   ylab = ylab, 
                   key.title = title(main=main.scale),
                   main=main,
                   color.palette = use_palette, 
                   asp = 1, 
                   plot.axes = {
                     axis(1, seq(gmin, gmax, by = gby),
                          cex.axis = 1.3,
                          cex = 1.3,
                          cex.lab = 1.3);
                     axis(2, seq(gmin, gmax, by = gby),
                          cex.axis = 1.3,
                          cex = 1.3,
                          cex.lab = 1.3); 
                     abline(h = seq(gmin, gmax, by = gby),
                            col = "white", lty = "dashed", lwd = 2); 
                     abline(v = seq(gmin, gmax, by = gby),
                            col = "white", lty = "dashed", lwd = 2);
                     abline(a = 0, b = 1, col = "black", lty = 4, lwd = 4)
                   })
  }
}

  
  
#------------------------------------------------------------


# graficzna konwergencja beta

KonwergencjaBetaGraf <- function(zbior_danych, zmienna_wzrost, zmienna_warpocz, zmienna_jednostka,odstep=0.2,
                               xlab="Wzglêdny PKB na mieszkañca w okresie pocz¹tkowym",
                               ylab="Tempo wzrostu", 
                               main="Analiza konwergencji typu beta") {
  
  plot(zmienna_warpocz,zmienna_wzrost,pch=20,cex=2.5,xlab=xlab,ylab=ylab,cex.axis=1.3,cex.lab=1.3)
  abline(lm(zmienna_wzrost~zmienna_warpocz), lwd=2)
  title(main=main,cex=2)
  text(zmienna_warpocz, zmienna_wzrost+odstep,zmienna_jednostka, cex=1.2)
  abline(h=0, lty="dashed"); abline (v=100, lty="dashed")
}

# regresja dla konwergencji beta

KonwergencjaBetaReg<-function(zbior_danych, zmienna_wzrost, zmienna_warpocz) {
  summary(lm(zmienna_wzrost~zmienna_warpocz))
}



# wspó³czynnik zmiennoœci
cv <- function(x)(100*sd(x)/mean(x)) 

# wspó³czynnik zmiennoœci wa¿ony
wt.cv <- function(x,w)(100*wt.sd(x,w)/wt.mean(x,w)) 



#------------------------------------------------------------
ergodicKDE <- function(kernel_)
{
  require(markovchain)
  
  tmA <- matrix(kernel_$density,
                nrow = sqrt(nrow(kernel_)),
                byrow = T)
  
  rowSumy <- rowSums(tmA)
  colSumy <- colSums(tmA)
  
  groups_ <- as.character(unique(kernel_$x))
  
  for (i in 1:ncol(tmA))
    tmA[i,] <- tmA[i,] / rowSumy[i]
  
  dtmcA <- new("markovchain",
               transitionMatrix = tmA,
               states = groups_)
  
  erg_ <- as.numeric(100*steadyStates(dtmcA))
  
  ergodic.dist <- data.frame(value = 100*as.numeric(groups_),
                             ergodic = erg_)
  
  return(ergodic.dist)  
}

ergodicKDE_cde <- function(kernel_)
{
  require(markovchain)
  
  tmA <- matrix(kernel_$z,
                nrow = length(kernel_$x),
                byrow = T)
  
  rowSumy <- rowSums(tmA)
  colSumy <- colSums(tmA)
  
  groups_ <- as.character(unique(kernel_$x))
  
  for (i in 1:ncol(tmA))
    tmA[i,] <- tmA[i,] / rowSumy[i]
  
  dtmcA <- new("markovchain",
               transitionMatrix = tmA,
               states = groups_)
  
  erg_ <- as.numeric(100*steadyStates(dtmcA))
  
  ergodic.dist <- data.frame(value = 100*as.numeric(groups_),
                             ergodic = erg_)
  
  return(ergodic.dist)  
}



#------------------------------------------------------------
# testowe !!!!

# !!!! test równoœci macierzy i wektorów ergodycznych

# matrix1 <- macierz_PKB
# matrix2 <- macierz_trwanie

test_trans_matrices <- function(matrix1, matrix2) {
  n_groups1 <- ncol(matrix1)
  n_groups2 <- ncol(matrix2)
  
  stopifnot(n_groups1==n_groups2)
    
  non_zero1 <- which(matrix1 != 0)
  non_zero2 <- which(matrix2 != 0)
  
  ni1 <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", rownames(matrix1)[1:n_groups1]))
  ni2 <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", rownames(matrix2)[1:n_groups2]))
  
  # dla matrix1 danej egzogenicznie
  chi2_stat1 <- sum((((matrix2/100 - matrix1/100)**2) * ni1)[non_zero1]/
                      (matrix1[non_zero1]/100), na.rm = T)
  
  # dla matrix2 danej egzogenicznie
  chi2_stat2 <- sum((((matrix1/100 - matrix2/100)**2) * ni2)[non_zero2]/
                      (matrix2[non_zero2]/100), na.rm = T)
  
  chi2_pvalue1 <- pchisq(chi2_stat1, length(non_zero1), lower.tail = FALSE)
  chi2_pvalue2 <- pchisq(chi2_stat2, length(non_zero2), lower.tail = FALSE) 
  
  return (list(chi2_stat1 = chi2_stat1,
              chi2_pvalue1 = chi2_pvalue1,
              chi2_stat2 = chi2_stat2,
              chi2_pvalue = chi2_pvalue2
              )
          )
  }

# wariant 2
  
test_trans_matrices2 <- function(matrix1, matrix2) {
  n_groups1 <- matrix1$n_groups
  n_groups2 <- matrix2$n_groups
  
  stopifnot(n_groups1 == n_groups2)
  
  non_zero1 <- which(matrix1$trans_probs != 0)
  non_zero2 <- which(matrix2$trans_probs != 0)
  
  ni1 <- matrix1$initial_counts
  ni2 <- matrix2$initial_counts
  
  # dla matrix1 danej egzogenicznie
  chi2_stat1 <- sum((((matrix2$trans_probs/100 - matrix1$trans_probs/100)**2) * ni1)[non_zero1]/
                      (matrix1$trans_probs[non_zero1]/100), na.rm = T)
  
  # dla matrix2 danej egzogenicznie
  chi2_stat2 <- sum((((matrix1$trans_probs/100 - matrix2$trans_probs/100)**2) * ni2)[non_zero2]/
                      (matrix2$trans_probs[non_zero2]/100), na.rm = T)
  
  chi2_pvalue1 <- pchisq(chi2_stat1, length(non_zero1), lower.tail = FALSE)
  chi2_pvalue2 <- pchisq(chi2_stat2, length(non_zero2), lower.tail = FALSE) 
  
  return (list(chi2_stat1 = chi2_stat1,
               chi2_pvalue1 = chi2_pvalue1,
               chi2_stat2 = chi2_stat2,
               chi2_pvalue = chi2_pvalue2)
          )
}

 
  

#------------------------------------------------------------
# testowe !!!



test_ergodic_vectors <- function(matrix1, matrix2) {
  
  n_groups1 <- ncol(matrix1)
  n_groups2 <- ncol(matrix2)
  
  stopifnot(n_groups1 == n_groups2)
  
  ergodic1 <- matrix1[n_groups1 + 1,]/100
  ergodic2 <- matrix2[n_groups2 + 1,]/100
  
  n_units <- sum(as.numeric(gsub(".*\\((.*)\\).*", "\\1", rownames(matrix1)[1:n_groups1])))
  
  chi_stat1 <- chisq.test(x = round(n_units*ergodic1),
                          p = ergodic2,
                          rescale.p = T,
                          simulate.p.value = T)
  
  chi_stat2 <- chisq.test(x = round(n_units*ergodic2),
                          p = ergodic1,
                          rescale.p = T,
                          simulate.p.value = T)
  
  return(list(chi_stat1, chi_stat2))
  
}

test_ergodic_vectors_Wilcox2003 <- function(matrix1, matrix2, method = "M", B = 1000) {
  
  if (!method %in% c("B", "M", "both") ) {warning("Provide correct method: \"M\", \"B\" or \"both\"!")
    stop}
  
  source("Wilcox2003_functions.R")
  
  n_groups1 <- matrix1$n_groups
  n_groups2 <- matrix2$n_groups
  
  stopifnot(n_groups1 == n_groups2)
  
  vec1 <- rep(1:n_groups1, matrix1$ergodic_counts)
  vec2 <- rep(1:n_groups2, matrix2$ergodic_counts)
  
  if (method != "B") result_M <- data.frame(binband(vec1, vec2, plotit = F))
  if (method != "M") {
    result_B <- disc2comSK(vec1, vec2, nboot = B, SEED = F)
    result_B <- data.frame(test = result_B$test, p.value = result_B$p.value)
  }
  
  if (method == "B") result_ <- list(B = result_B)
  if (method == "M") result_ <- list(M = result_M)
  if (method == "both") result_ <- list(M = result_M, B = result_B)
  
  return(result_)
}

#-----------------------------------------------

ergodicKDE <- function(kernel_) {
  require(markovchain)
  
  tmA <- matrix(kernel_$density,
                nrow = sqrt(nrow(kernel_)),
                byrow = T)
  
  rowSumy <- rowSums(tmA)
  colSumy <- colSums(tmA)
  
  groups_ <- unique(kernel_$x)
  
  for (i in 1:ncol(tmA))
    tmA[i,] <- tmA[i,] / rowSumy[i]
  
  dtmcA <- new("markovchain",
               transitionMatrix = tmA,
               states = as.character(groups_))
  
  erg_ <- as.numeric(steadyStates(dtmcA)) / (100*(max(groups_) - min(groups_))/length(groups_))
  
  ergodic.dist <- data.frame(value = 100 * as.numeric(groups_),
                             ergodic = erg_)
  
  return(ergodic.dist)  
}
