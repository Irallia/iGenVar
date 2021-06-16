library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("caller", "mapper", "subsample", "vcf", "score", "metric", "value"))
res$caller = factor(res$caller, levels=c('pbsv', 'Sniffles', 'SVIM'), labels=c('pbsv', 'Sniffles', 'SVIM'))
res$subsample = recode(res$subsample, pooled.subsampled.10 = 10, pooled.subsampled.20 = 20, pooled.subsampled.30 = 30, pooled.subsampled.40 = 40, pooled.subsampled.50 = 50, pooled.subsampled.60 = 60, pooled.subsampled.70 = 70, pooled.subsampled.80 = 80, pooled.subsampled.90 = 90, pooled = 100)

final <- res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    filter(vcf == args[2]) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    group_by(caller, mapper, subsample, vcf) %>%
    filter(f1 == max(f1)) %>%
    summarise(bestf1 = max(f1, na.rm = TRUE), bestscore = first(score))

final%>%
    ggplot(aes(subsample, bestf1, color=caller, pch=caller)) +
      geom_path(, alpha=0.7) +
      geom_point(size=1.0, alpha=0.7) +
      scale_shape_manual(values=c(15,16,17)) +
      scale_color_manual(values=c("deepskyblue3", "goldenrod2", "firebrick2")) +
      scale_x_continuous(breaks=seq(0,100,20), minor_breaks=seq(10,90,20), limits=c(0,100)) +
      scale_y_continuous(breaks=seq(0,100,20), minor_breaks=seq(10,90,20), limits=c(0,100)) +
      labs(y = "Highest F1-score", x = "Dataset subsample", color = "Tool", pch = "Tool") +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=20))
ggsave(args[3], width=5.5, height=4)
final %>% write_tsv(args[4])
