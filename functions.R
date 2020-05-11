
# Compute vector mean count from vector window
vmc <- function(vm.win){
  mean(abs(vm.win - mean(vm.win)))
}

## Function to compute local maxima
## source: https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
localMaxima <- function(x) {
  #cat("x is ",x, "\n")
  #y <- diff(c(-.Machine$integer.max, x)) > 0L
  y <- diff(c(-Inf, x)) > 0L
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

## Function which cut x vector at indices given by x.cut.idx,
## approximate cut parts into common length and average 
## parts point-wise into one vector
cut.and.avg <- function(x, x.cut.idx,length_out = 200){
  x.cut.idx.l <- length(x.cut.idx)
  mat.out <- matrix(NA, nrow = x.cut.idx.l-1, ncol = length_out)
  for (i in 1:(length(x.cut.idx)-1)){
    xp <- x[x.cut.idx[i]:x.cut.idx[i+1]]
    xpa <- approx(seq(0, 1, length.out = length(xp)), xp, seq(0, 1, length.out = length_out))$y
    mat.out[i, ] <- as.numeric(scale(xpa, scale = TRUE, center = TRUE))
  }
  out <- apply(mat.out, 2, mean)
  as.numeric(scale(out, scale = TRUE, center = TRUE))
}

out.plot1 <- function(val, out, fs = 30, xlims = c(0,30), save_plots = T, ind = this_ind){
  plot_name = paste0("pattern_eval_", this_ind, ".pdf")
  yrange <- c(-1, 1) * max(abs(val))
  y.h <- 0
  plt <- ggplot()
  for (i in 1:nrow(out)){
    cat("I is ", i, "\n")
    tau1_i <- out[i, "tau_i"]
    tau2_i <- tau1_i + out[i, "T_i"] - 1
    tau1_i <- tau1_i/fs
    tau2_i <- tau2_i/fs
    plt <- 
      plt + 
      geom_vline(xintercept = tau1_i, color = "red") + 
      geom_vline(xintercept = tau2_i, color = "red") + 
      annotate(
        "rect",
        fill = "pink", 
        alpha = 0.3,
        xmin = tau1_i, 
        xmax = tau2_i, 
        ymin = yrange[1],
        ymax = yrange[2]
      )
  }
  geom_line.df <- data.frame(x = seq(0, by = 1/fs, length.out = length(val)), y = val)
  plt <- 
    plt + 
    geom_line(data = geom_line.df, 
              aes(x = x, y = y), 
              color = "black", 
              size = 0.3) + 
    theme_bw(base_size = 9) + 
    labs(x = "Time [s]", y = "Black line: x",
         title = "Black line: signal x\nRed vertical lines: start and end points of identified pattern occurrence\nRed shaded area: area corresponding to identified pattern occurrence") +
    xlim(xlims)
  #plot(plt)
  if(save_plots == T){
    save_plot(plot_name, plt, base_width = 12, base_height = 3)
  }
}
  