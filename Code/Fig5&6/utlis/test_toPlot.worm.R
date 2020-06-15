
ggplot(tmp_df,
       aes(x = as.numeric(tmp_df$dnds),
           y= as.numeric(tmp_df$logexpr))) +
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  xlab(label = "dnds") + 
  ylab(label = "logexpr")

ggplot(toPlot.worm,
       aes(x = as.numeric(toPlot.worm$dNdS),
           y= as.numeric(toPlot.worm$logexpr))) +
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  xlab(label = "dnds") + 
  ylab(label = "logexpr")