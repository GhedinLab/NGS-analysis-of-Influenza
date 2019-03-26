library(ggplot2)
# h1 <- read.csv(file="H1N1.plotme.txt", header=TRUE, sep=",")

# h3 <- read.csv(file="h3n2.plotme.txt", header=TRUE, sep=",")
# # H1N1.plotme.txt  h3n2.plotme.txt
# p <- ggplot(h1, aes(x=group, y=value)) + 
# geom_boxplot() +
# geom_jitter(shape=16, position=position_jitter(0.2)) +
# ylim(0,25)

# ggsave('h1.pdf',units='in',width=4,height=6)


# p <- ggplot(h3, aes(x=group, y=value)) + 
# geom_boxplot() +
# geom_jitter(shape=16, position=position_jitter(0.2)) +
# ylim(0,25)

# ggsave('h3.pdf',units='in',width=4,height=6)



h1 <- read.csv(file="H1.H3.redemux.reseq.plotme.txt", header=TRUE, sep=",")

# h3 <- read.csv(file="h3n2.plotme.txt", header=TRUE, sep=",")
# H1N1.plotme.txt  h3n2.plotme.txt
p <- ggplot(h1, aes(x=group, y=value)) + 
geom_boxplot() +
facet_grid(strain + type ~ group, space = 'free', scales='free' ) +
geom_jitter(shape=16, position=position_jitter(0.2)) +
ylim(0,25)
ggsave('all.pdf',units='in',width=4,height=16,  useDingbats=FALSE)
