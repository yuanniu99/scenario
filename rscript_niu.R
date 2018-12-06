library(car)
library(ggplot2)
library(dplyr)


# Import data
df.effc <- read.csv('data/efficacy.csv')
df.rand <- read.csv('data/randomization.csv')
df.subj <- read.csv('data/subject.csv')

# Merge
df <- merge(merge(df.effc, df.rand, by = 'subject'), df.subj, by = 'subject')
nrow(df)
colnames((df))
summary(df)

# dropput
hist(df$duration)
max(df$duration)
table(df$duration)
print('dropout rate:')
sum(df$duration < 365)/sum(df$duration == 365)

  ## remove dropout
df <- df[df$duration == 365,]

# descriptive and missing
summary(df)
cor(df[!is.na(df$mucus.viscosity),c('nosebleeds', 'previous.year', 'mucus.viscosity')])
# +++++++++++++++++++++++ Exploratary +++++++++++++++++++++++
# 1. treatment
# plot of nosebleeds by arm
boxplot(nosebleeds ~ arm, data=df, main="Treatment Effect on Nosebleed", 
        xlab="Treatment Level", ylab="Nosebleed Rate")

qplot(sample = nosebleeds, data = df, color=arm)

# ck homogeneity of variance between groups
leveneTest(nosebleeds ~ arm , df)

summary(aov(glm(nosebleeds ~ arm, data=df, family = 'poisson')))
df %>% group_by(arm) %>% summarize(mean_eff = mean(nosebleeds, na.rm = T))

# 2. mucus viscosity
boxplot(mucus.viscosity ~ arm, df)

df[!is.na(df$mucus.viscosity),] %>% 
  ggplot() +
  aes(x = mucus.viscosity, y = nosebleeds) +
  geom_point() +
  geom_smooth(method = 'lm')

summary(glm(nosebleeds ~ mucus.viscosity, df, family = 'poisson'))

## interaction of mucus viscosity and arm
summary(aov(glm(nosebleeds ~ arm * mucus.viscosity, df, family = 'poisson')))

df1 <- df[!is.na(df$mucus.viscosity),]
x <- df1$mucus.viscosity
df1$viscosity_3group <- 
  case_when(x > mean(x) + sd(x) ~ 'high',
            x < mean(x) + sd(x) & x > mean(x) - sd(x) ~ 'average',
            x < mean(x) - sd(x) ~ 'low')
count(df1, viscosity_3group)

df1 %>% ggplot() +
  aes(x = arm, y = nosebleeds, group = viscosity_3group, color = viscosity_3group) +
  geom_point(color = 'grey') +
  geom_smooth (method = 'lm', se = FALSE, fullrange = TRUE)

# 3. country
table(df$country, df$arm) 

## difference of mean value of nosebleeds across countries
df %>% group_by(country) %>% 
  summarise(mean_eff = mean(nosebleeds)) %>%
  ggplot(aes(x = country, y = mean_eff, group=1)) +
  geom_point() + 
  geom_line() +
  labs(title = 'Mean Nosebleed Rate by Country', y = '# Severe Nosebleeds Annually', x = 'Country')


fit2 <- aov(glm(nosebleeds ~ arm * country, data=df, family = 'poisson'))
summary(fit2)

interaction.plot(x.factor = df$country,
                 trace.factor = df$arm,
                 response = df$nosebleeds,
                 type = "b",
                 col = c('red', 'blue'))
## understand country difference of nosebleed level
posthoc <- TukeyHSD(fit2, "country", conf.level = 0.95) 
print(posthoc)
plot(posthoc)

# 4. effect of paper tissue
df %>% group_by(tissue.use) %>%
  summarize(mean_nosebleed_rate = mean(nosebleeds)) %>%
  ggplot()+
  aes(x = tissue.use, y = mean_nosebleed_rate, group = 1) +
  geom_line() +
  geom_point()   

## interaction effect of tissue and arm
summary(aov(glm(nosebleeds ~ arm * tissue.use, df, family = 'poisson')))

interaction.plot(x.factor = df$tissue.use,
                 trace.factor = df$arm,
                 response = df$nosebleeds,
                 
                 type = "b",
                 col = c('red', 'blue'))

# 5. previous year nosebleed level
boxplot(previous.year ~ arm, data=df)
df%>%group_by(arm) %>% summarize(meaneff = mean(previous.year))
 
df %>% ggplot() +
  aes(x = previous.year, y = nosebleeds) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(x = '# of Severe Nosebleed in Previous Year', y = 'Nosebleed Rate')

summary(aov(glm(nosebleeds ~ arm * previous.year, df, family = 'poisson')))
interaction.plot(x.factor = as.factor(df$previous.year),
                 trace.factor = df$arm,
                 response = df$nosebleeds,
                 
                 type = "b",
                 col = c('red', 'blue'))


# 6. eye color
table(df$eye.colour)
summary(aov(nosebleeds ~ eye.colour, df))

# +++++++++++++++++++++++++++ Model +++++++++++++++++++++++++++++++

# poisson regression to model count of nosebleed, using country as blocking factor and previous year condition as covariate
mod1 <- glm(nosebleeds ~ 1 + arm * mucus.viscosity + country + previous.year, data = df1, family = 'poisson')
summary(mod1)

  ## goodness of fit
pchisq(mod1$deviance, mod1$df.residual, lower.tail = F)

mod1$deviance/mod1$df.residual

  ## residual
dev.off()
par(mfrow=c(1,2))
res <- residuals(mod1, type="deviance")
plot(exp(predict(mod1)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)

# treatment effect
anov1 <- Anova(aov(mod1), type = 'III')
print(anov1)

# predicted
result <- data.frame(df1, mod1$fitted)
