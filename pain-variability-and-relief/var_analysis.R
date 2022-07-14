library(ggplot2)
library(patchwork)
library(performance)
library(boot)
library(ggeffects)
library(pwr)
library(emmeans)

#####
# Load and prepare data
pl1 <- read.csv("data/placeboI_ratings_updated.csv")
pl2 <- read.csv("data/placeboII_ratings.csv")[,-2]

pl1.groups <- rio::import("data/placeboI_groups.xlsx")
colnames(pl1.groups)[ncol(pl1.groups)] <- "group"
pl1.groups$group <- factor(pl1.groups$group, levels = c(3,1), labels = c("No treatment", "Placebo"))
pl2.groups <- rio::import("data/placeboII_groups.xlsx")
colnames(pl2.groups)[1] <- "ID"
pl2.groups$group <- factor(pl2.groups$group, levels = c(0,1,2), labels = c("No treatment", "Placebo", "Drug"))


pl1_subs <- na.omit(rio::import("data/placebo_lists.xlsx")[,1])
pl2_subs <- na.omit(rio::import("data/placebo_lists.xlsx")[,2])

#####
# Define functions

rsp <- function(object) {
  df <- object$df.residual
  object <- summary(object)
  ts <- coef(object)[,3] 
  rsq <- object$r.squared
  
  sign(ts) * sqrt(ts^2 * (1-rsq) / df)
}

p <- function(object) {
  df <- object$df.residual
  object <- summary(object)
  ts <- coef(object)[,3] 
  
  pt(ts, df)
}

rsp.boot <- function(data, ids, model.formula) {
  model <- lm(model.formula, data[ids,])
  atanh(rsp(model))
}

adjusted_response <- function(model, variable, by = NA) {
  if(!is.na(by)) {
    mf <- model.frame(model)
    mf[,variable] <- mf[,variable] - mean(mf[,variable])
    model <- update(model, data = mf)
    
    if(is.numeric(variable))
      variable <- colnames(model.frame)[variable]
    indices <- which(grepl(variable,names(coef(model))))
    betas <- coef(model)[indices]
    x <- model.matrix(model)[,indices]
    ybar <- mean(model$model[,1])
    
    ybar + x %*% betas + resid(model)
  } else {
    beta <- coef(model)[variable]
    x <- model$model[,variable]
    ybar <- mean(model$model[,1])
    
    ybar + beta * (x - mean(x)) + resid(model)
  }
}

rsp.loo <- function(model) {
  df <- model.frame(model)
  
  rsps <- lapply(1:nrow(df), function(i) {
    model.tmp <- lm(formula(model), df[-i,])
    rsp(model.tmp)
  })
  do.call(rbind, rsps)
}

extract_ratings <- function(ratings, groups) {
  ratings.out <- lapply(unique(ratings$subject), function(sub) {
    
    tmp <- na.omit(subset(ratings, subject == sub))[,c(2,3)]
    group <- groups[groups$ID == sub,"group"]
    if(length(group) == 1) {
      
      #tmp <- aggregate(pain ~ floor(as.numeric(day)), tmp, mean)
      colnames(tmp) <- c("day","pain")
      pre <- mean(tmp[tmp$day < 7,"pain"])
      preSD <- sd(tmp[tmp$day < 7,"pain"])
      last_day <- as.numeric(tail(tmp$day,1))
      post <- mean(tmp[tmp$day > last_day-7,"pain"])
      data.frame(subject = sub,
                 group = group,
                 preSD = preSD,
                 pre = pre,
                 post = post)
    }
  })
  do.call(rbind.data.frame, ratings.out)
}

#####
# Get ratings and variabilities

pl1.ratings <- na.omit(extract_ratings(pl1, pl1.groups))
pl2.ratings <- na.omit(extract_ratings(subset(pl2, subject!="NPL041" & subject != "NPL108"), pl2.groups))
contrasts(pl1.ratings$group) <- MASS::contr.sdif(2)
contrasts(pl2.ratings$group) <- MASS::contr.sdif(3)

#####
# Demographics
pl1.included <- subset(pl1.groups, ID %in% pl1.ratings$subject)
aggregate(Age ~ group, pl1.included, function(x) round(cbind(mean(x),sd(x))))
aggregate(Sex ~ group, pl1.included, function(x) cbind(sum(x==2),round(mean(x == 2),2)*100))


pl2.demos <- rio::import("data/placeboII_dmg.xlsx", sheet=3)
pl2.included <- subset(pl2.demos, SID %in% pl2.ratings$subject)
aggregate(Age_years ~ pl2.ratings$group, pl2.included, function(x) round(cbind(mean(x),sd(x))))
aggregate(Sex ~ pl2.ratings$group, pl2.included, function(x) cbind(sum(x=="Female"),round(mean(x == "Female"),2)*100))

###### 
# Get change in R^2 for adding preSD as a linear predictor (no interaction)

round(summary(lm(post ~ pre + group + preSD, pl1.ratings))$r.squared - 
        summary(lm(post ~ pre + group, pl1.ratings))$r.squared, 2)
round(summary(lm(post ~ pre + group + preSD, pl2.ratings))$r.squared - 
        summary(lm(post ~ pre + group, pl2.ratings))$r.squared, 2)


#####
# Get partial R^2s

# placebo 1
model1 <- lm(post ~ pre + group * preSD, pl1.ratings)
model1.c <- lm(post ~ pre + group * I(preSD-mean(preSD)), pl1.ratings)

# print beta coefs (table 1)
round(cbind(coef(model1.c),confint(model1.c)), 1)
check_model(model1.c)

#print semipartial correlations (table 1)
rsp(model1.c)
set.seed(0)
boot.out <- boot(data = pl1.ratings, 
                 statistic = rsp.boot, 
                 R = 1000, 
                 model.formula = formula(post ~ pre + as.factor(group) * I(preSD-mean(preSD)), pl1.ratings))
tanh(boot.ci(boot.out, index = 4, type="bca")$bca[,4:5])
tanh(boot.ci(boot.out, index = 5, type="bca")$bca[,4:5])

#power analysis for groups in study -> should hit 0.8 or above
p(model1.c)
pwr.r.test(n = 20, r = -0.22, sig.level= 0.026)
pwr.r.test(n = 20, r = 0.22, sig.level= 0.969)



# within group, no pre covariate (supp table 1)
lapply(unique(pl1.ratings$group), function(grp) {
  tmp <- subset(pl1.ratings, group == grp)
  model <- lm(post-pre ~ 1 + I(preSD-mean(preSD)), tmp)
  
  boot.out <- boot(data = tmp, 
                   statistic = rsp.boot, 
                   R = 1000, 
                   model.formula = formula(post-pre  ~ 1 +  I(preSD-mean(preSD)), tmp))
  
  list(grp,
       round(cbind(coef(model),confint(model))[2,],1),
       round(c(rsp(model)[2],
               tanh(boot.ci(boot.out, index = 2, type="bca")$bca[,4:5])),2))
})


# within group, pre covariate (supp table 1)
lapply(unique(pl1.ratings$group), function(grp) {
  tmp <- subset(pl1.ratings, group == grp)
  model <- lm(post-pre ~ pre + I(preSD-mean(preSD)), tmp)
  
  boot.out <- boot(data = tmp, 
                   statistic = rsp.boot, 
                   R = 1000, 
                   model.formula = formula(post-pre ~ pre +  I(preSD-mean(preSD)), tmp))
  
  list(grp,
       round(cbind(coef(model),confint(model))[3,],1),
       round(c(rsp(model)[3],
               tanh(boot.ci(boot.out, index = 3, type="bca")$bca[,4:5])),2))
})

predicted.model1 <- ggemmeans(model1, c("preSD","group"))


#figure 1
placebo1.plot <- 
  ggplot() +
  geom_point(data = model1$model,
             aes(x = preSD, 
                 y = adjusted_response(model1,"preSD",T),
                 color = group),
             alpha = 0.5) +
  geom_ribbon(data = predicted.model1,
              aes(x = x,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = group),
              alpha = 0.1) +
  geom_line(data = predicted.model1,
            aes(x = x,
                y = predicted,
                color = group),
            size = 1) +
  ylab("Post-intervention pain") +
  xlab(expression(SD[baseline])) +
  labs(color = "Group", fill = "Group") +
  coord_cartesian(ylim=c(0,10)) +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", n = 3, type="continuous")) +
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1", n = 3, type="continuous")) +
  theme_minimal() +
  guides(fill = "none", color = "none") +
  theme(panel.background = element_rect(fill = "#f9f9f9",color = "#ffffff"),
        strip.text = element_text(size = 8),
        legend.position = "none") +
  ggtitle("Placebo I")

# loo analysis
loo1 <- rbind(data.frame(group = "No treatment", 
                         subject = 1:length(rsp.loo(model1)[,4]), 
                         effect = rsp.loo(model1)[,4]),
              data.frame(group = "Placebo", 
                         subject = 1:length(rsp.loo(model1)[,5]), 
                         effect = rsp.loo(model1)[,5]))
loo1.plot <-
  ggplot() +
  geom_hline(yintercept = 0,
             color = "grey") +
  geom_point(data = loo1,
             aes(y = effect,
                 x = subject),
             size = 0.7) +
  geom_point(data = subset(loo1, subject == 50),
             aes(y = effect,
                 x = subject),
             color = "red") +
  ylab(expression(r[sp])) +
  xlab("Patient #") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",color = "#ffffff"),
        strip.text = element_text(size = 8),
        legend.position = "none") +
  facet_wrap(. ~ group) +
  ggtitle("Placebo I")



#placebo 2  
model2 <- lm(post ~  pre + group * preSD, pl2.ratings)
model2.c <- lm(post ~  pre + group * I(preSD-mean(preSD)), pl2.ratings)
check_model(model2.c)

#print beta coefs (table 1)
round(cbind(coef(model2.c),confint(model2.c)),1)

#no pre covariate
lapply(unique(pl2.ratings$group), function(grp) {
  tmp <- subset(pl2.ratings, group == grp)
  model <- lm(post-pre ~ 1 + I(preSD-mean(preSD)), tmp)
  
  boot.out <- boot(data = tmp, 
                   statistic = rsp.boot, 
                   R = 1000, 
                   model.formula = formula(post-pre ~ 1 +  I(preSD-mean(preSD)), tmp))
  
  list(grp,
       round(cbind(coef(model),confint(model))[2,],1),
       round(c(rsp(model)[2],
               tanh(boot.ci(boot.out, index = 2, type="bca")$bca[,4:5])),2))
})


# with pre covariate
lapply(unique(pl2.ratings$group), function(grp) {
  tmp <- subset(pl2.ratings, group == grp)
  model <- lm(post-pre ~ pre + I(preSD-mean(preSD)), tmp)
  
  boot.out <- boot(data = tmp, 
                   statistic = rsp.boot, 
                   R = 1000, 
                   model.formula = formula(post ~ pre +  I(preSD-mean(preSD)), tmp))
  
  list(grp,
       round(cbind(coef(model),confint(model))[3,],1),
       round(c(rsp(model)[3],
               tanh(boot.ci(boot.out, index = 3, type="bca")$bca[,4:5])),2))
})

#print semi-partial correlations (table 1)
round(rsp(model2),2)
par(mfrow=c(1,3))
plot((rsp.loo(model2)[,5]), ylim=c(-.25,.25))
plot((rsp.loo(model2)[,6]), ylim=c(-.25,.25))
plot((rsp.loo(model2)[,7]), ylim=c(-.25,.25))

# print confidence intervals for semi-partial correlations (table 1)
boot.out <- boot(data = pl2.ratings, 
                 statistic = rsp.boot, 
                 R = 1000, 
                 model.formula = formula(post ~ pre + as.factor(group) * I(preSD-mean(preSD)), pl1.ratings))
tanh(boot.ci(boot.out, index = 5, type="bca")$bca[,4:5])
tanh(boot.ci(boot.out, index = 6, type="bca")$bca[,4:5])
tanh(boot.ci(boot.out, index = 7, type="bca")$bca[,4:5])

predicted.model2 <- ggemmeans(model2, c("preSD","group"))

# 
placebo2.plot <- 
  ggplot() +
  geom_point(data = model2$model,
             aes(x = preSD, 
                 y = adjusted_response(model2,"preSD",T),
                 color = group),
             alpha = 0.5) +
  geom_ribbon(data = predicted.model2,
              aes(x = x,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = group),
              alpha = 0.1) +
  geom_line(data = predicted.model2,
            aes(x = x,
                y = predicted,
                color = group),
            size = 1) +
  ylab("Post-intervention pain") +
  xlab(expression(SD[baseline])) +
  labs(color = "Group", fill = "Group") +
  coord_cartesian(ylim=c(0,10)) +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", n = 3, type="continuous")) +
  scale_fill_manual( values = wesanderson::wes_palette("Zissou1", n = 3, type="continuous")) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",color = "#ffffff"),
        strip.text = element_text(size = 8),
        legend.position = "bottom") +
  ggtitle("Placebo II")

placebo1.plot + placebo2.plot + plot_layout(guides = "collect") &
  theme(legend.position='bottom')

#ggsave("figures/figure1.pdf", device=cairo_pdf, width = 7, height = 3.5)
ggsave("figures/figure1.png", width = 7, height = 3.5)

loo2 <- rbind(data.frame(group = "No treatment", 
                         subject = 1:length(rsp.loo(model2)[,5]), 
                         effect = rsp.loo(model2)[,5]),
              data.frame(group = "Placebo", 
                         subject = 1:length(rsp.loo(model2)[,6]), 
                         effect = rsp.loo(model2)[,6]),
              data.frame(group = "Drug", 
                         subject = 1:length(rsp.loo(model2)[,7]), 
                         effect = rsp.loo(model2)[,7]))
loo2$group <- factor(loo2$group, levels = unique(loo2$group))
loo2.plot <-
  ggplot() +
  geom_hline(yintercept = 0,
             color = "grey") +
  geom_point(data = loo2,
             aes(y = effect,
                 x = subject),
             size = 0.7) +
  geom_point(data = subset(loo2, subject == 14),
             aes(y = effect,
                 x = subject),
             color = "red") +
  ylab(expression(r[sp])) +
  xlab("Patient #") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",color = "#ffffff"),
        strip.text = element_text(size = 8),
        legend.position = "none") +
  facet_wrap(. ~ group) +
  ggtitle("Placebo II")

loo1.plot + loo2.plot + plot_layout(widths = c(2,3))

#ggsave("figures/supp.figure1.pdf", device=cairo_pdf, width = 10, height = 2.5)
ggsave("figures/supp.figure1.png", width = 10, height = 2.5)



combined <- rbind(cbind(study = "Placebo I", pl1.ratings),
                  cbind(study = "Placebo II", pl2.ratings))
combined$delta <- combined$post - combined$pre
ggplot(data = combined,
       aes(x = preSD, y = delta, color = group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="lm",
              aes(fill = group),
              alpha = 0.1) +
  facet_grid(group ~ study) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",color = "#ffffff"),
        strip.text = element_text(size = 8),
        legend.position = "none") +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", n = 3, type="continuous")) +
  scale_fill_manual( values = wesanderson::wes_palette("Zissou1", n = 3, type="continuous")) +
  xlab(expression(SD[baseline])) +
  ylab("Change in pain")
ggsave("figures/supp.figure2.png", width = 5, height = 5)

