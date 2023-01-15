#!/usr/bin/env Rscript
#
# Supplementary Figures and statistics
#

library(ggplot2)
library(cowplot)
library(ggrepel)

################################################################################
# Create Output Scaffolding
################################################################################

s <- ifelse(!dir.exists('./fig'), dir.create('./fig'), FALSE)
s <- ifelse(!dir.exists('./fig/supp'), dir.create('./fig/supp'), FALSE)

s <- ifelse(!dir.exists('./calc'), dir.create('./calc'), FALSE)
s <- ifelse(!dir.exists('./calc/supp'), dir.create('./calc/supp'), FALSE)

################################################################################
# Supplmentary Figure -- Pairwise Frequency Correlations (Train; Test; Predict)
################################################################################

df <- read.csv('./calc/merged_substitution_info.csv')
df$location <- sapply(df$substitution, function(x) { 
  as.numeric(substr(x, 2, str_length(x)-1)) })
df$group <- factor(df$group, levels=c('Train', 'Test', 'Predict'), ordered=T)

for (l in list(c('Train', 'Test'), c('Train', 'Predict'), c('Test', 'Predict'))) {
  f1 <- l[[1]]
  f2 <- l[[2]]
  
  common.subs.df <- as.data.frame.matrix(table(df$substitution, df$group))
  common.subs <- rownames(common.subs.df)[common.subs.df[,f1] == 1
                                          & common.subs.df[,f2] == 1]
  f1.sub <- subset(df, group == f1)
  f2.sub <- subset(df, group == f2)
  
  merged.df <- data.frame(
    substitution = common.subs,
    f1 = f1.sub$frequency[match(common.subs, f1.sub$substitution)],
    f2 = f2.sub$frequency[match(common.subs, f2.sub$substitution)]
  )
  merged.df$f1.scaled <- merged.df$f1 / sum(merged.df$f1)
  merged.df$f2.scaled <- merged.df$f2 / sum(merged.df$f2)
  
  ggplot(merged.df, aes(x=log10(f1), y=log10(f2))) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x) +
    theme_cowplot() +
    background_grid()
  ggsave(paste0('./fig/supp/', f1, '_', f2, '_corr.png'), width=8, height=7.5)
  
  print(paste0('=====', f1, ' - ', f2, ' ====='))
  print(cor.test(merged.df$f1, merged.df$f2, method = 'spearman'))
  print('==========')
}





