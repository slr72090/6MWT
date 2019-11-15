## Analyze six min walk test data ## ----------------------------------------------
## Sylvia Ranjeva, 10/2019

devtools::install_github("martakarass/adept")

## Load packages
require(dplyr)
require(reshape2)
require(lubridate)
require(gridExtra)
require(ggplot2)
require(cowplot)
require(magrittr)
require(adept)

## Function scripts
source("functions.R")

## LOAD SAMPLE DATA ------------------------------------------------------------------------------------------

# IF running on a high-performance computing cluster: will incorporate later as volume of data increases --------------------------------------------------------
#args = commandArgs(trailingOnly=TRUE)
#subjID = as.numeric(args[1])
#filename <- paste0("walk_",subjID,".csv")

#------------------------------------------------------------------------------------------------------------
n_files = 10
filenames = paste0("./6MWT/0",c(1:n_files),"_pocket.csv")
df <- data.frame()
for(i in c(1:n_files)){
  filename <- filenames[i]
  if(file.exists(filename)){
  this_df <- read.csv(filename) %>% 
  mutate(date_time= as.POSIXct(Time, format = "%Y-%m-%dT%H:%M:%SZ"),
         t_rel = Seconds - Seconds[1],
         vm = sqrt(X^2 + Y^2 + Z^2))
  df <- rbind(df, this_df)
  }
}

## ANALYZE AND PLOT RAW DATA ----------------------------------------------------------------------------------

## Plotting specs  ----------------------------
max_time_display =90 #number of seconds to plot on x axis 
textSize = 10
save_plots = T
source("plot_themes.R")

## Plot raw accelerometry  ------------------------------
dfm_raw <- df %>% 
  select(c(ID, t_rel, X,Y,Z)) %>% 
  melt(id.vars = c("ID", "t_rel"))

p_raw <- ggplot(dfm_raw,aes(x = t_rel, y = value, color = variable)) + 
  geom_line() +
  xlim(0,max_time_display) + 
  xlab("Relative time") + ylab("Acceleration") + 
  labs(color = "Accelerometer axis") + 
  facet_grid(ID~.) + 
  plot_themes

## Plot vector magnitude -----------------------------------------------
dfm_vm <- df %>% 
  select(c(ID, t_rel, vm))

p_vm <- ggplot(dfm_vm,aes(x = t_rel, y = vm)) + 
  geom_line() +
  xlim(0,max_time_display) + 
  xlab("Relative time") + ylab("Vector Magnitude") + 
  facet_grid(ID~.) + 
  plot_themes 

## Combined plot ------------------------------------------------
p_combined_raw = plot_grid(p_raw, p_vm, ncol = 1)
if(save_plots){
  save_plot("raw_data.pdf",p_combined_raw, base_width = 12, base_height = 6)
}

## Calculate vector mean counts --------------------------------------------------
# Set specs -----------------------------------------
n_sec = 2 #Length of averaging window 
sampling_freq = 30 #30hz
this_ind = "dr" #specify individual 
resting_thresh = 0.3 #count vmcs below this value as "resting"
length_vec = nrow(df)
cutoff = 2 #chop n seconds off of the start and end 

## Compute vmc vector in n-second windows ------------
vm <- df$vm
win.vl <- sampling_freq * n_sec
rn.seq <- seq(1, to = length(vm), by = win.vl)
vmc.vec <- sapply(rn.seq, function(rn.i){
  vm.win.idx <- rn.i : (rn.i + win.vl - 1)
  vm.win <- vm[vm.win.idx]
  vmc(vm.win)
})

vmc.df <- data.frame(vmc = vmc.vec,
                     t_rel = df$t_rel[rn.seq],
                     rn_seq = rn.seq,
                     ID = df$ID[rn.seq]) 
vmc.df.rest <- 
  vmc.df %>% 
  filter(ID == this_ind) %>%
  mutate(resting = ifelse(vmc < resting_thresh, 1, 0),
         vmc_tau_i = rn_seq) %>%  #- length_vec) %>%
  select( t_rel, resting, vmc_tau_i)

vmc.df <- vmc.df %>% left_join(vmc.df.rest, by = "t_rel") %>% 
  filter(t_rel < max(df$t_rel) - cutoff, 
         t_rel > min(df$t_rel) + cutoff)

## Plot vmc 
p_vmc <- ggplot(vmc.df, aes(x = t_rel/60, y = vmc)) +
  facet_grid(ID ~.) + 
  geom_tile(aes(fill = factor(resting)),  height = Inf, alpha = 0.1) +
  geom_line(size = 0.3) +
  labs(x = "Exercise time [min]", y = "Vector magnitue count",
       title = paste0("Vector magnitude count (vmc) computed over ", n_sec, "second-length windows of (vm)"),
       fill = "Resting: ") + 
  scale_fill_manual(values = c("white", "blue")) + 
  plot_themes

if(save_plots){
  save_plot("vector_mean_counts.pdf", p_vmc, base_width = 12, base_height = 3)
}

## GENERATE TEMPLATE(S) ## -----------------------------------------------------------------------------------------------

# Smooth the accelerometry time series 
window_length = 0.1 # smoothing window length in seconds 
df$vm_smoothed1 <- windowSmooth(x = df$vm, x.fs = sampling_freq, W = window_length)

## Obtain individual-specific subset of data  
sub.ind <- df %>% filter(ID == this_ind)

## Identify local maxima from individual data, plot, and hand-pick representative peaks 
length_plot = 200 # Number of indices to plot 
x1 <- sub.ind[, "vm_smoothed1"]
x1 <- x1[!is.na(x1)]
x1.locMax <- localMaxima(x1)

df_x1 = data.frame(x = c(1:length_plot), y = x1[1:length_plot])

## Plot 1: All maxima 
p1 <- ggplot(df_x1, aes(x = x, y = y)) + geom_line() + 
  geom_vline(xintercept = x1.locMax, color = "red") + 
  xlim(0, length_plot) + 
  xlab("Index") + 
  ylab("") + 
  labs(title = "(vm) local maxima") + 
  plot_themes

## Manually identify stride peaks
stride_peak_vec = c(1, 4, 8, 13, 18)

## Plot 2: Manually identified maxima
p2 <- ggplot(df_x1, aes(x = x, y = y)) + geom_line() + 
  geom_vline(xintercept = x1.locMax[stride_peak_vec], color = "red") + 
  xlim(0, length_plot) + 
  xlab("Index") + 
  ylab("") + 
  labs(title = "(vm) local maxima subset") +
  plot_themes

if(save_plots){
  p_combined <- plot_grid(p1,p2)
  save_plot("stride_peaks.pdf", p_combined, base_width = 12, base_height = 3)
}

## Generate template by selecting n representative strides, interpolating smoothed vm trajectory for each, and averaging/scaling across the n strides 
template.ind1 <- cut.and.avg(x1, x1.locMax[stride_peak_vec])
df_template <- data.frame(x = c(1:length(template.ind1)), y = template.ind1)

p_template <- ggplot(df_template, aes(x = x, y = y)) + 
  geom_line(col = "red") + 
  labs(title = "Template") + 
  xlab("Index") + 
  ylab("") + 
  plot_themes 

if(save_plots){
  save_plot("template.pdf", p_template, base_width = 8, base_height = 3)
}

## SEGMENT DATA ACCORDING TO TEMPLATE ----------------------------------------------------------------------------

# Output segmented data 
out.ind1 <- segmentPattern(x = df$vm,
                          x.fs = sampling_freq,
                          template = template.ind1,
                          x.adept.ma.W = 0.1,
                          finetune = "maxima",
                          finetune.maxima.nbh.W = 0.3,
                          pattern.dur.seq = seq(0.5, 1.8, length.out = 50),
                          similarity.measure = "cor",
                          similarity.measure.thresh = 0.8,
                          compute.template.idx = TRUE,
                          run.parallel = F)

hist(out.ind1$sim_i, 100)

## Visualize segemetntation
if(save_plots){
  out.plot1(val = df$vm, out = out.ind1)
}

## Separate and overlay walking strides 
n_sec = 3
x.ind <- df$vm

# Merge output with raw data and eliminate periods where pt was resting 
plt.df.ind <- 
  out.ind1 %>% 
  merge(vmc.df.rest) %>%
  filter(tau_i >= vmc_tau_i, tau_i < vmc_tau_i + n_sec* sampling_freq, resting == 0) %>%
  mutate(ind = "ind1")

#plt.df <- plt.df.ind

## For data frame #1 (raw vm segments)
stride.acc.vec.ind <- numeric()
stride.tau_i.vec.ind <- numeric()
stride.idx.vec.ind <- numeric()
## For data frame #2 (scaled vm segments)
stride_S.acc.vec.ind <- numeric()
stride_S.tau_i.vec.ind <- numeric()
stride_S.phase.vec.ind <- numeric()

for (i in 1:nrow(plt.df.ind)){
  out.i <- plt.df.ind[i, ]
  x.ind.i <- x.ind[out.i$tau_i : (out.i$tau_i + out.i$T_i - 1)]
  x.ind.i.len <- length(x.ind.i)
  if (var(x.ind.i) < 1e-3) next
  ## For data frame #1 
  stride.acc.vec.ind   <- c(stride.acc.vec.ind, x.ind.i)
  stride.tau_i.vec.ind <- c(stride.tau_i.vec.ind, rep(out.i$tau_i, x.ind.i.len))
  stride.idx.vec.ind <- c(stride.idx.vec.ind, 1:x.ind.i.len)
  ## For data frame #2
  x.ind.i_S <- approx(x = seq(0, 1, length.out = length(x.ind.i)),
                     y = x.ind.i,
                     xout = seq(0, 1, length.out = 200))$y
  x.ind.i_S <- as.numeric(scale(x.ind.i_S))
  stride_S.acc.vec.ind <- c(stride_S.acc.vec.ind, x.ind.i_S)
  stride_S.tau_i.vec.ind <- c(stride_S.tau_i.vec.ind, rep(out.i$tau_i, 200))
  stride_S.phase.vec.ind <- c(stride_S.phase.vec.ind, seq(0, 1, length.out = 200))
}

## data frame #1 
stride.df.ind <- data.frame(acc = stride.acc.vec.ind, 
                           tau_i = stride.tau_i.vec.ind,
                           idx = stride.idx.vec.ind)
## data frame #2
stride_S.df.ind <- data.frame(acc = stride_S.acc.vec.ind, 
                             tau_i = stride_S.tau_i.vec.ind,
                             phase = stride_S.phase.vec.ind)

## Plot segmented walking strides
plt1 <- 
  stride.df.ind %>%
  ggplot(aes(x = idx/sampling_freq, y = acc, group = tau_i)) + 
  geom_line(alpha = 0.1) + 
  theme_bw(base_size = 8) + 
  labs(x = "Stride pattern duration [s]", y = "Vector magnitude [g]",
       title = "Segmented walking strides")
plt2 <- 
  stride_S.df.ind %>%
  ggplot(aes(x = phase, y = acc, group = tau_i)) + 
  geom_line(alpha = 0.1) + 
  theme_bw(base_size = 8) + 
  labs(x = "Stride pattern phase", y = "Vector magnitude (scaled) [g]",
       title = "Segmented walking strides, aligned and scaled") 

if(save_plots){
  p_combined <- plot_grid(plt1, plt2)
  save_plot("segmented_strides.pdf", p_combined, base_width = 12, base_height = 3)
}

## CORRELATION CLUSTERING OF SEGMENTED WALKING STRIDES ## ----------------------------------------------------------------------------------

## Compute strides distance martrix
stride_S.dfdc.ind <- reshape(stride_S.df.ind, timevar = "tau_i", idvar = "phase", direction = "wide")[,-1]
colnames(stride_S.dfdc.ind) = unique(stride_S.df.ind$tau_i)
data.mat.tau_i <- as.numeric(colnames(stride_S.dfdc.ind))
data.mat <-  apply(as.matrix(stride_S.dfdc.ind),2,as.numeric)
D.mat  <- dist(cor(data.mat))

# Compute optimal number of clusters : for now, use 1 cluster
#n_clust = fviz_nbclust(cor(data.mat), cluster::pam, method = "silhouette")

## Get cluster medoids
cluster.k <- 1
medoids.idx <- round(seq(1, ncol(stride_S.dfdc.ind), length.out = cluster.k + 2))
medoids.idx <- medoids.idx[-c(1, medoids.idx + 2)]

## Cluster strides
set.seed(1) #random number generator
pam.out <- cluster::pam(D.mat, cluster.k, diss = TRUE, medoids = medoids.idx)
table(pam.out$clustering)

## Put clustering results into data frame 
data.df <- as.data.frame(t(data.mat))
colnames(data.df) <- seq(0, to = 1, length.out = 200)
data.df$tau_i <- data.mat.tau_i
data.df$cluster <- pam.out$clustering
data.dfm <- melt(data.df, id.vars = c("tau_i", "cluster"))
data.dfm$variable <- as.numeric(as.character(data.dfm$variable))
data.dfm$cluster <- paste0("cluster ", data.dfm$cluster)

data.dfm.agg <- 
  data.dfm %>%
  group_by(variable, cluster) %>% 
  summarise(value = mean(value))

p_clustered <- ggplot(data.dfm, aes(x = variable, y = value, group = tau_i)) + 
  geom_line(alpha = 0.2) + 
  geom_line(data = data.dfm.agg,  aes(x = variable, y = value, group = 1),
            color = "red", size = 1, inherit.aes = FALSE) +  
  facet_grid(cluster ~ .) + 
  theme_bw(base_size = 9) + 
  labs(x = "Stride pattern phase", y = "Vector magnitude (scaled) [g]",
       title = "Segmented walking strides, aligned, scaled, clustered\nRed line: point-wise mean") 

if(save_plots){
  save_plot("clustered_strides.pdf", p_clustered, base_width = 8, base_height = 6)
}

df_strides <- data.dfm %>%
  select(tau_i, cluster)  %>%
  distinct() %>%
  left_join(plt.df.ind, by = "tau_i") 

p_stride_length <-  ggplot(df_strides, aes(x = tau_i / (sampling_freq * 60) , y = T_i / sampling_freq, color = cluster)) + 
  geom_tile(data = vmc.df.rest, 
            aes(x = vmc_tau_i / (sampling_freq * 60), y = 1, fill = factor(resting)), 
            height = Inf, alpha = 0.1, inherit.aes = FALSE) +
  geom_point(alpha = 0.4) + 
  theme_bw(base_size = 9) + 
  labs(x = "Exercise time [min]", y = "Estimated stride duration time [s]",
       color = "Stride assignment: ",
       fill = "Resting: ") + 
  scale_fill_manual(values = c("white", "blue")) + 
  theme(legend.position = "top",
        legend.background = element_rect(fill = "grey90"))  + 
  scale_y_continuous(limits = c(0.5, 1.8))

if(save_plots){
  save_plot("stride_durations.pdf", p_stride_length, base_width = 4, base_height = 3)
}