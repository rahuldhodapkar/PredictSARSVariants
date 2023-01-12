#!/usr/bin/env Rscript
#
# Manuscript Figures
#

library(ggplot2)
library(cowplot)
library(Biostrings)
library(RColorBrewer)
library(stringr)

################################################################################
# Create Output Scaffolding
################################################################################

s <- ifelse(!dir.exists('./fig'), dir.create('./fig'), FALSE)
s <- ifelse(!dir.exists('./calc'), dir.create('./calc'), FALSE)

################################################################################
# Load Data
################################################################################

df <- read.csv('./calc/variant_meta_df.csv')
df$location <- sapply(df$Mutation, function(x) { 
  as.numeric(substr(x, 2, str_length(x)-1)) })

data('BLOSUM80')
data('PAM30')

################################################################################
# Build Plots
################################################################################

###############################################
## Show Density Plot of Substitution Locations
###############################################

df.sub <- subset(df, location >= 438 & location <= 506)

plot.df <- data.frame(
  group = c(rep('Train', sum(df.sub$in_train == 'True')),
            rep('Test', sum(df.sub$in_test == 'True')),
            rep('Predict', sum(df.sub$in_predict == 'True'))),
  location = c(
    df.sub$location[df.sub$in_train == 'True'],
    df.sub$location[df.sub$in_test == 'True'],
    df.sub$location[df.sub$in_predict == 'True']
  )
)

ggplot(plot.df, aes(x=location, fill=group)) + 
  geom_density(alpha=0.4, adjust=0.1) +
  scale_fill_manual(values=c('#648fff', '#785ef0', '#dc267f')) +
  xlim(c(438, 506)) +
  facet_grid(group ~ .) +
  theme_cowplot() + 
  background_grid(minor='x')


ggplot(plot.df, aes(x=location, fill=group)) + 
  geom_histogram(alpha=0.4, position='dodge') +
  scale_fill_manual(values=c('#648fff', '#785ef0', '#dc267f')) +
  theme_cowplot()      
            
###############################################
## Compare using PAM and BLOSUM
###############################################

set.seed(42)

rand.blosum.scores <- c()
rand.pam.scores <- c()
for (i in 1:sum(df$in_predict == 'True')) {
  ix <- sample(1:25, 2)
  rand.blosum.scores <- c(rand.blosum.scores, BLOSUM80[ix[1], ix[2]])
  rand.pam.scores <-c(rand.pam.scores, PAM30[ix[1], ix[2]])
}

plot.df <- data.frame(
  group = c(rep('Train', sum(df$in_train == 'True')),
            rep('Test', sum(df$in_test == 'True')),
            rep('Predict', sum(df$in_predict == 'True')),
            rep('Random', sum(df$in_predict == 'True'))),
  blosum80 = c(
    df$BLOSUM80[df$in_train == 'True'],
    df$BLOSUM80[df$in_test == 'True'],
    df$BLOSUM80[df$in_predict == 'True'],
    rand.blosum.scores
  ),
  pam30 = c(
    df$PAM30[df$in_train == 'True'],
    df$PAM30[df$in_test == 'True'],
    df$PAM30[df$in_predict == 'True'],
    rand.pam.scores
  )
)

plot.df$group <- factor(plot.df$group,
                        levels=c('Random', 'Train', 'Test', 'Predict'),
                        ordered=T)

ggplot(plot.df, aes(fill=group, y=blosum80)) + 
  geom_boxplot() +
  scale_fill_manual(values=c('#dddddd', '#648fff', '#785ef0', '#dc267f')) +
  theme_cowplot()
ggsave('./fig/blosum80_boxplot.png', width=8, height=4)


# statistics
wilcox.test(blosum80 ~ group, data=subset(plot.df, group %in% c('Predict', 'Random')))
wilcox.test(blosum80 ~ group, data=subset(plot.df, group %in% c('Train', 'Random')))
wilcox.test(blosum80 ~ group, data=subset(plot.df, group %in% c('Test', 'Random')))

ggplot(plot.df, aes(fill=group, y=pam30)) + 
  geom_boxplot() +
  scale_fill_manual(values=c('#dddddd', '#648fff', '#785ef0', '#dc267f')) +
  theme_cowplot()

wilcox.test(pam30 ~ group, data=subset(plot.df, group %in% c('Predict', 'Random')))
wilcox.test(pam30 ~ group, data=subset(plot.df, group %in% c('Train', 'Random')))
wilcox.test(pam30 ~ group, data=subset(plot.df, group %in% c('Test', 'Random')))
ggsave('./fig/pam30_boxplot.png', width=8, height=4)


