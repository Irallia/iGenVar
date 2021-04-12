library(tidyverse)
library(ggrepel)
library(scales)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("parameters", "caller", "mapper", "subsample", "vcf", "score", "metric", "value"))

res2 <- res %>%
  filter(metric %in% c("recall", "precision")) %>%
  pivot_wider(names_from=metric, values_from=value) %>%
  filter(recall!=0 | precision!=0) %>%
  mutate(precision = 100*precision, recall = 100*recall, f1 = 2*precision*recall/(precision+recall)) %>%
  group_by(parameters, caller, mapper, subsample, vcf) %>%
  summarise(best_f1 = max(f1)) %>%
  pivot_wider(names_from=vcf, values_from=best_f1) %>%
  separate(parameters, into=c("pmd", "pdn", "edn", "cmd"), sep="_", remove=F)

ggplot(res2, aes(giab, giab.gt, label=parameters, group=edn, color=edn)) +
  geom_point(aes(size=cmd)) +
  geom_label_repel() +
  geom_path() +
  facet_wrap(vars(subsample), ncol=1) +
  coord_cartesian(xlim=c(90, 93.5), ylim=c(80, 90))

ggsave(args[2], width=20, height=12)
