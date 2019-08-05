# Mae
# July 22
# Red size histogram


df <- read.csv("data/red_urchin_sizes.csv")

hist(df$size)

df$cuts <- cut(df$size, 5)

ggplot(df)+
  geom_bar(aes(x = cuts))

table(df$cuts)
range(df$size)

df$cuts <- cut(df$size, breaks = c(0, 20, 40, 60, 80, 100, 122))


# 20-40, 40-60, 60-80, 80-100, 100-120