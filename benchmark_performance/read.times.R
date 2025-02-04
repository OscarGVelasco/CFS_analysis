
PATH="/home/o872o/o872o/cluster_similarity_sc/benchmark/"
df <- data.frame()

## MULTIDATASET
library(dplyr)
library(ggplot2)
df <- data.frame()

for(i in seq(3,56,by=2)){
  #load(file = "/home/o872o/o872o/cluster_similarity_sc/benchmark/time.benchmarking.results.RData")
  load(file = file.path(PATH, paste0("time.3.datasets.CPU.",i,".RData")))
  df <- rbind.data.frame(df, time)
}

load(file = file.path(PATH, paste0("time.3.datasets.CPU.",1,".RData")))
df <- rbind.data.frame(df, c(time_tmp,1), c(time_tmp,1),c(time_tmp,1))

# load(file = file.path(PATH, paste0("time.3.datasets.CPU.",56,".RData")))
# df <- rbind.data.frame(df, c(time,56))
df$clock <- df$elapsed*df$cpus
df$clock <- df$user.self / (df$user.self + df$sys.self) * 100 * df$cpus

df$cpus <- factor(df$cpus)
df.mean <- df %>% group_by(cpus) %>% dplyr::summarize(Mean = mean(elapsed, na.rm=TRUE)/60,)
df.minutes <- df %>% group_by(cpus) %>% dplyr::mutate(minutes = elapsed/60)

library(ggplot2)
pdf(file = "/home/o872o/o872o/cluster_similarity_sc/benchmark/boxplot.mean.times.cpu.processing.3.datasets.750k.cells.pdf",width = 8, height = 4)
g <- ggplot(df, aes(x=cpus, y=elapsed/60)) +
  geom_boxplot(col="#D5BADB") +
  geom_jitter(width = 0.1, col="#86608E") +
  #geom_line(col="blue", linewidth=1) +
  theme_minimal() +
  #annotate('text', x=df$cpus[ order(df$elapsed)][1], y=(min(df$elapsed)/60)+0.6, label = round(min(df$elapsed)/60,digits = 1)) +
  geom_text(data = df.mean, aes(x = cpus, y = Mean+1, label = round(Mean,digits = 1))) +
  #ggplot2::scale_x_continuous(breaks = df$cpus) +
  ggplot2::scale_y_continuous(breaks = seq(0, round(max(df$elapsed/60))+1, by=1)) +
  geom_hline(yintercept=(min(df$elapsed)/60)-0.01, linetype="dashed",
             color = "#F08080", linewidth=0.5) +
  #ylim(c(0,round(max(df$elapsed/60))+1)) +
  ylab("minutes") +
  ggtitle("Computing time for 750k cells divided in 3 datasets \n and with 20 cluster per dataset")
print(g)
dev.off()


ggplot(df, aes(x=cpus, y=clock)) +
  geom_boxplot(col="#D5BADB") +
  geom_jitter(width = 0.1, col="#86608E") +
  #geom_line(col="blue", linewidth=1) +
  theme_minimal() +
  #annotate('text', x=df$cpus[ order(df$elapsed)][1], y=(min(df$elapsed)/60)+0.6, label = round(min(df$elapsed)/60,digits = 1)) +
  geom_text(data = df.mean, aes(x = cpus, y = Mean+1, label = round(Mean,digits = 1))) +
  #ggplot2::scale_x_continuous(breaks = df$cpus) +
  ggplot2::scale_y_continuous(breaks = seq(0, round(max(df$elapsed/60))+1, by=1)) +
  geom_hline(yintercept=(min(df$elapsed)/60)-0.01, linetype="dashed",
             color = "#F08080", linewidth=0.5) +
  #ylim(c(0,round(max(df$elapsed/60))+1)) +
  ylab("minutes") +
  ggtitle("Computing time for 750k cells divided in 3 datasets \n and with 20 cluster per dataset")


###################################
## Size Benchmark

library(dplyr)
library(ggplot2)
df <- data.frame()

for(i in c(5000, 25000, 50000, 100000, 150000, 200000, 250000)){
  
  #load(file = file.path(PATH, paste0("time.3.datasets.SIZE.",i,".1.using.17.CPUs.RData")))
  #df <- rbind.data.frame(df, time)
  load(file = file.path(PATH, paste0("time.3.datasets.SIZE.",i,".2.using.17.CPUs.RData")))
  df <- rbind.data.frame(df, time)
}

df$size2 <- as.integer(df$size*3)
df$size2 <- factor(df$size2)

pdf(file = "/home/o872o/o872o/cluster_similarity_sc/benchmark/boxplot.mean.times.size.cells.processing.3.datasets.17.CPUs.pdf",width = 8, height = 4)
g <- ggplot(df, aes(x=size2, y=elapsed)) +
  geom_boxplot(col="#E5D8BD") +
  geom_jitter(width = 0.1, col="#A65628",size=2) +
  #geom_line(col="blue", linewidth=1) +
  theme_minimal() +
  #annotate('text', x=df$cpus[ order(df$elapsed)][1], y=(min(df$elapsed)/60)+0.6, label = round(min(df$elapsed)/60,digits = 1)) +
  #geom_text(data = df.mean, aes(x = cpus, y = Mean+1, label = round(Mean,digits = 1))) +
  #ggplot2::scale_x_continuous(breaks = df$cpus) +
  ggplot2::scale_y_continuous(breaks = seq(0, round(max(df$elapsed))+20, by=10)) +
  #ylim(c(0,round(max(df$elapsed/60))+1)) +
  ylab("seconds") +
  xlab("number of cells") +
  ggtitle("Computing time for 750k cells divided in 3 datasets \n and with 20 cluster per dataset")
print(g)
dev.off()
