library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("caller", "mapper", "subsample", "vcf", "score", "metric", "value"))
res$caller = factor(res$caller, levels=c('pbsv', 'Sniffles', 'SVIM'), labels=c('pbsv', 'Sniffles', 'SVIM'))


res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    filter(subsample=="pooled", vcf==args[2]) %>%
    ggplot(aes(recall, precision, color=caller, pch=caller)) +
      geom_point(size=1.0) +
      scale_shape_manual(values=c(15,16,17)) +
      scale_color_manual(values=c("deepskyblue3", "goldenrod2", "firebrick2")) +
      scale_x_continuous(breaks=seq(0,100,20), minor_breaks=seq(10,90,20), limits=c(0,100)) +
      scale_y_continuous(breaks=seq(0,100,20), minor_breaks=seq(10,90,20), limits=c(0,100)) +
      geom_path() +
      labs(y = "Precision", x = "Recall", color = "Tool", pch = "Tool") +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=20))

ggsave(args[3], width=5.5, height=4)
