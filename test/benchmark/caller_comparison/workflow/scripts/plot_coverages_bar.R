library(tidyverse)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("caller", "mapper", "coverage", "zygosity", "type", "score", "metric", "value"))
res$type = factor(res$type, levels=c('del','inv','dup','ins'), labels=c('Deletions','Inversions','Tandem Dupl.','Insertions'))
res$caller = factor(res$caller, levels=c('pbsv', 'Sniffles', 'SVIM'), labels=c('pbsv', 'Sniffles', 'SVIM'))
res$coverage = factor(res$coverage, levels=c(5, 15, 30, 45, 60), labels=c('5x', '15x', '30x', '45x', '60x'))

res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    mutate(f1 = 2*precision*recall/(precision+recall)) %>%
    group_by(caller, mapper, coverage, zygosity, type) %>%
    summarise(bestf1 = max(f1, na.rm = TRUE)) %>%
    ggplot(aes(coverage, bestf1, fill=caller, alpha=zygosity)) +
      geom_col(position=position_dodge(), color="black", size=0.1) +
      scale_fill_manual(values=c("deepskyblue3", "goldenrod2", "firebrick2")) +
scale_alpha_discrete(range = c(0.2, 1.0)) +
      facet_wrap(~type) +
      labs(y = "Best F1", x = "Coverage", fill = "Tool", alpha = "Zygosity") +
      lims(y=c(0,100)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))

ggsave(args[2], width=5, height=4)


