#!/usr/bin/env Rscript
#
# Manuscript Figures
#

library(ggplot2)
library(cowplot)
library(Biostrings)
library(RColorBrewer)
library(stringr)
library(ggrepel)

################################################################################
# Create Output Scaffolding
################################################################################

s <- ifelse(!dir.exists('./fig'), dir.create('./fig'), FALSE)
s <- ifelse(!dir.exists('./calc'), dir.create('./calc'), FALSE)

################################################################################
# Load Data
################################################################################

data('BLOSUM80')
data('BLOSUM100')
data('PAM30')
data('PAM250')

################################################################################
# Build Plots
################################################################################

###############################################
## Show Density Plot of Substitution Locations
###############################################

df <- read.csv('./calc/merged_substitution_info.csv')
df$location <- sapply(df$substitution, function(x) { 
  as.numeric(substr(x, 2, str_length(x)-1)) })
df$group <- factor(df$group, levels=c('Train', 'Test', 'Predict'), ordered=T)

plt.mutation.density <- ggplot(df, aes(x=location, fill=group, weight=log10(frequency))) + 
  geom_density(adjust=0.1) +
  scale_fill_manual(values=c('#648fff', '#785ef0', '#dc267f')) +
  facet_grid(group ~ .) +
  theme_cowplot() + 
  scale_x_continuous(
    breaks=seq(from=440, to=506, by=5),
    minor_breaks=seq(from=438, to=506, by=1),
    limits=c(438, 506)
  ) +
  background_grid(minor='x')
ggsave('./fig/mutation_density_plot.png', width=12, height=5)

ggplot(df, aes(x=location, fill=group, weight=log10(frequency))) + 
  geom_density(alpha=0.4, adjust=0.1) +
  geom_histogram(aes(y = stat(density)), color='#000000', alpha=0.4, binwidth = 1) +
  scale_fill_manual(values=c('#648fff', '#785ef0', '#dc267f')) +
  facet_grid(group ~ .) +
  theme_cowplot() + 
  scale_x_continuous(
    breaks=seq(from=440, to=506, by=5),
    minor_breaks=seq(from=438, to=506, by=1),
    limits=c(438, 506)
  ) +
  background_grid(minor='x')
ggsave('./fig/mutation_histogram_plot.png', width=12, height=5)


###############################################
## Compare using PAM and BLOSUM
###############################################
df <- read.csv('./calc/variant_meta_df.csv')
df$location <- sapply(df$Mutation, function(x) { 
  as.numeric(substr(x, 2, str_length(x)-1)) })

set.seed(42)

rand.blosum80.scores <- c()
rand.pam30.scores <- c()
for (i in 1:sum(df$in_predict == 'True')) {
  ix <- sample(1:24, 2)
  rand.blosum80.scores <- c(rand.blosum80.scores, BLOSUM80[ix[1], ix[2]])
  rand.pam30.scores <-c(rand.pam30.scores, PAM30[ix[1], ix[2]])
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
    rand.blosum80.scores
  ),
  pam30 = c(
    df$PAM30[df$in_train == 'True'],
    df$PAM30[df$in_test == 'True'],
    df$PAM30[df$in_predict == 'True'],
    rand.pam30.scores
  )
)

plot.df$group <- factor(plot.df$group,
                        levels=c('Random', 'Train', 'Test', 'Predict'),
                        ordered=T)

#
# Make Plots
#

ggplot(plot.df, aes(fill=group, x=group, y=blosum80)) + 
  geom_boxplot() +
  scale_fill_manual(values=c('#dddddd', '#648fff', '#785ef0', '#dc267f')) +
  theme_cowplot()
ggsave('./fig/blosum80_boxplot.png', width=5, height=4)

wilcox.test(blosum80 ~ group, data=subset(plot.df, group %in% c('Predict', 'Random')))
wilcox.test(blosum80 ~ group, data=subset(plot.df, group %in% c('Train', 'Random')))
wilcox.test(blosum80 ~ group, data=subset(plot.df, group %in% c('Test', 'Random')))

ggplot(plot.df, aes(fill=group, x=group, y=pam30)) + 
  geom_boxplot() +
  scale_fill_manual(values=c('#dddddd', '#648fff', '#785ef0', '#dc267f')) +
  theme_cowplot()
ggsave('./fig/pam30_boxplot.png', width=5, height=4)

wilcox.test(pam30 ~ group, data=subset(plot.df, group %in% c('Predict', 'Random')))
wilcox.test(pam30 ~ group, data=subset(plot.df, group %in% c('Train', 'Random')))
wilcox.test(pam30 ~ group, data=subset(plot.df, group %in% c('Test', 'Random')))

################################################################################
# Create Pie / Donut Charts
################################################################################
df <- read.csv('./calc/merged_substitution_info.csv')
df$location <- sapply(df$substitution, function(x) { 
  as.numeric(substr(x, 2, str_length(x)-1)) })

sub.venn.df <- as.data.frame.matrix(table(df$substitution, df$group))

sub.venn.df[sub.venn.df$Predict == 1
            & sub.venn.df$Train == 0
            & sub.venn.df$Test == 1,]
sub.venn.df$ID <- paste0(sub.venn.df$Predict, sub.venn.df$Train, sub.venn.df$Test)

predict.df <- subset(df, group == 'Predict')
predict.df$ID <- sub.venn.df$ID[match(predict.df$substitution, rownames(sub.venn.df))]
predict.df$ID <- factor(predict.df$ID, levels=c('101', '110', '111', '100'), 
                        ordered=T)

plt.pred.identities <- ggplot(predict.df, aes(x=location, fill=ID)) +
  geom_histogram(binwidth = 1, color='#000000', size=0.5) +
  scale_fill_manual(values=c('#785ef0', '#648fff', '#fe6100', '#dc267f')) +
  scale_x_continuous(
    breaks=seq(from=440, to=506, by=5),
    minor_breaks=seq(from=438, to=506, by=1),
    limits=c(438, 506)
  ) +
  theme_cowplot() + 
  background_grid(minor='xy')
ggsave('./fig/prediction_identity_plot.png', width=12, height=2)

plt.merged <- plot_grid(
  plt.mutation.density + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
  plt.pred.identities,
  ncol=1, align='v', axis='rl', rel_heights = c(3.4,1))
ggsave('./fig/merged_substitution_plots.png', width=12, height=7.3)

################################################################################
# Create Predicted Mutation Burden Plot
################################################################################
pred.7df4 <- read.csv('./data/7DF4_StoACE2.tsv', sep='\t', comment.char='#')
pred.8d8q.2130 <- read.csv('./data/8D8Q_Sto2130.tsv', sep='\t', comment.char='#')
pred.8d8q.2196 <- read.csv('./data/8D8Q_Sto2196.tsv', sep='\t', comment.char='#')

colnames(pred.7df4) <- paste0(colnames(pred.7df4), '_7df4')
colnames(pred.8d8q.2130) <- paste0(colnames(pred.8d8q.2130), '_8d8q_2130')
colnames(pred.8d8q.2196) <- paste0(colnames(pred.8d8q.2196), '_8d8q_2196')

merged.df <- merge(pred.7df4, pred.8d8q.2130, by.x='Mutation_7df4', by.y='Mutation_8d8q_2130')
merged.df <- merge(merged.df, pred.8d8q.2196, by.x='Mutation_7df4', by.y='Mutation_8d8q_2196')

merged.df$DDG_8d8q_avg <- sapply(
  1:nrow(merged.df),
  function(i) {
     mean(merged.df$DDG_8d8q_2130[[i]], merged.df$DDG_8d8q_2196[[i]])
  })

merged.df$ID <- sub.venn.df$ID[match(merged.df$Mutation_7df4, rownames(sub.venn.df))]
merged.df$ID <- factor(merged.df$ID, levels=c('101', '110', '111', '100'), 
                        ordered=T)

ggplot(merged.df, aes(x=-DDG_7df4, y=-DDG_8d8q_avg, fill=ID)) +
  scale_fill_manual(values=c('#785ef0', '#648fff', '#fe6100', '#dc267f')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color='#000000', shape=21, size=3, alpha=0.5) +
  theme_cowplot() +
  xlim(values=c(-3,3)) + 
  ylim(values=c(-1.5,1.5)) +
  background_grid(minor='xy') +
  theme(
    axis.line = element_blank()
  ) +
  geom_text_repel(aes(label=Mutation_7df4), max.overlaps = Inf)
ggsave('./fig/predicted_binding_energies_labeled.png', width=9, height=8.5)


ggplot(merged.df, aes(x=-DDG_7df4, y=-DDG_8d8q_avg, fill=ID)) +
  scale_fill_manual(values=c('#785ef0', '#648fff', '#fe6100', '#dc267f')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color='#000000', shape=21, size=5, alpha=0.5) +
  theme_cowplot() +
  xlim(values=c(-3,3)) + 
  ylim(values=c(-1.5,1.5)) +
  background_grid(minor='xy') +
  theme(
    axis.line = element_blank()
  )
ggsave('./fig/predicted_binding_energies.png', width=9, height=8.5)


ggplot(merged.df, aes(x=-DDG_7df4, y=-DDG_8d8q_avg, fill=ID)) +
  scale_fill_manual(values=c('#785ef0', '#648fff', '#fe6100', '#dc267f')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(color='#000000', shape=21, size=5, alpha=0.5) +
  theme_cowplot() +
  xlim(values=c(-3,3)) + 
  ylim(values=c(-1.5,1.5)) +
  background_grid(minor='xy') +
  theme(
    axis.line = element_blank()
  ) +
  geom_text_repel(aes(label=Mutation_7df4), max.overlaps = Inf, nudge_x = 2, nudge_y = -2,
                  data=merged.df[merged.df$DDG_7df4 < 0 & merged.df$DDG_8d8q_avg > 0,])
ggsave('./fig/predicted_binding_energies_keylabel.png', width=9, height=8.5)


