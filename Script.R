## What we miss on the road: underreporting in roadkill surveys depending on animal size and observer attentiveness

# Load libraries ------------------------------
library(lme4)
library(performance)
library(ggplot2)
library(DHARMa)
library(readr)
library(plyr)
library(dplyr)


# Import and clean data -----------------------
data <- read_csv("data.csv")
data <- data[data$Test %in% c("AG","CA","Birieux"),] # relevant experiments only, A_Loire and A_Ain are test runs of the uninformed method
data <- data[!is.na(data$Detect), ] # NA means participant did not answer
data$dummy <- 1

## Species in french vernacular;
## Naif => uninformed participant; 
## TV == road annual average Traffic Volume;
## AG (general assembly); CA (board of directors); Birieux (location of informed experiment);
## Seat: C (driver), P (passenger);
## Phase: Informed experiment 2 phases (mid-point stop to list species)
## Bodysize == approx measure of taxidermy volume


# Convert variables to factors
data$Sp <- as.factor(data$Sp)
data$Taxa <- as.factor(data$Taxa)
data$Obs <- as.factor(data$Obs)
data$Phase <- as.factor(data$Phase)
data$Test <- as.factor(data$Test)
data$Naive <- as.factor(data$Naive)
data$ObsID <- as.factor(paste(data$Obs, "_", data$Test, sep = ""))
data$Seat <- as.factor(data$Seat)

summary(as.factor(data$ObsID))


# Journey time histograms for informed ---------------------
hist(as.numeric(data$Journey_time1)/60,main="Journey time for phase 1 (mn)")
hist(as.numeric(data$Journey_time2)/60,main="Journey time for phase 2 (mn)")


## Full detection model (Uninformed == naive vs Informed) -------
summary(data$Test)

probs <- aggregate(Detect ~ interaction(Sp, Naive), data = data, FUN = sum)
probs$Obs <- aggregate(dummy ~ interaction(Sp, Naive), data = data, FUN = sum)[, 2]
probs$Score <- probs$Detect / probs$Obs

body_mass <- aggregate(Bmkg ~ interaction(Sp, Naive), data = data, FUN = mean)
probs$BodyMass <- body_mass$Bmkg
probs$SE <- sqrt(probs$Score * (1 - probs$Score) / probs$Obs)
probs$Naive <- sub(".*\\.", "", as.character(probs$`interaction(Sp, Naive)`))

colnames(probs) <- c("Sp", "Seen", "Obs", "Score", "BodyMass", "SE", "Naive")

plot(probs$Score ~ probs$BodyMass,
 col = as.factor(probs$Naive),
 pch = 15,
 xlab = "Body mass (kg)",
 ylab = "Detection score (red = uninformed)",
 log = "x")

probs$Taxa <- c("Mamm","Ois","Mamm","Ois","Ois","Ois","Mamm","Mamm","Ois","Mamm","Mamm","Ois","Mamm","Ois","Mamm")
probs$Naive<- as.factor(probs$Naive)
levels(probs$Naive) <- c("Naive", "Informed")

## Observers performance -----------------------
observers <- aggregate(Detect ~ ObsID + Naive, data, FUN = sum)
observers$tot <- aggregate(dummy ~ ObsID + Naive, data, FUN = sum)[, 3]
observers$Detect <- observers$Detect / observers$tot

mm_informed <- mean(observers$Detect[observers$Naive == "Informe"])
se_informed <- sqrt((mm_informed * (1 - mm_informed)) / 17)

mm_informed
mm_informed + se_informed
mm_informed - se_informed

mm_naive <- mean(observers$Detect[observers$Naive == "Naif"])
se_naive <- sqrt((mm_naive * (1 - mm_naive)) / 78)

mm_naive
mm_naive - se_naive
mm_naive + se_naive

ggplot(aes(x = Detect), data = observers[observers$Naive == "Informe", ]) +
 geom_histogram(bins = 30) +
 xlim(-0.1, 1) +
 theme_classic() +
 theme(axis.text = element_text(size = 14),
 axis.title = element_text(size = 14)) +
 xlab("Percentage of decoys successfully detected") +
 ylab("Number of observers") +
 geom_vline(xintercept = mm_informed, lty = 1) +
 geom_vline(xintercept = mm_informed + se_informed, lty = 3) +
 geom_vline(xintercept = mm_informed - se_informed, lty = 3)

ggplot(aes(x = Detect), data = observers[observers$Naive == "Naif", ]) +
 geom_histogram(bins = 10) +
 xlim(-0.1, 1) +
 theme_classic() +
 theme(axis.text = element_text(size = 14),
 axis.title = element_text(size = 14)) +
 xlab("Percentage of decoys successfully detected") +
 ylab("Number of observers") +
 geom_vline(xintercept = mm_naive, lty = 1) +
 geom_vline(xintercept = mm_naive + se_naive, lty = 3) +
 geom_vline(xintercept = mm_naive - se_naive, lty = 3)


## Are passengers better at detecting? ---------

# for naive experiments (in informed ones, all were drivers)

# mean detection by passengers
mp <- mean(data$Detect[data$Naive == "Naif" & data$Seat == "P"], na.rm = TRUE)
mp
sqrt(mp * (1 - mp) / length(data$Detect[data$Naive == "Naif" & data$Seat == "P"]))

# mean detection by drivers
md <- mean(data$Detect[data$Naive == "Naif" & data$Seat == "C"], na.rm = TRUE)
md
sqrt(md * (1 - md) / length(data$Detect[data$Naive == "Naif" & data$Seat == "C"]))

## GLM ----

## Detection ~ bodymass
model_lin <- glmer(as.factor(Detect) ~ Bmkg * Taxa + log(Bmkg) * Naive + Seat + (1 | ObsID),
 family = binomial,
 data)
summary(model_lin)
simulationOutput <- simulateResiduals(fittedModel = model_lin, plot = T)


## Detection ~ log(bodymass)
model <- glmer(as.factor(Detect) ~ log(Bmkg) * Taxa + log(Bmkg) * Naive + Seat + (1 | ObsID),
 family = binomial,
 data)
summary(model)
simulationOutput <- simulateResiduals(fittedModel = model, plot = T)
testOutliers(simulationOutput)

check_collinearity(model)


AIC(model) - AIC(model_lin)

## Final plot --------------

model_naive <- glmer(as.factor(Detect) ~ log(Bmkg) + Naive + Seat + (1 | ObsID),
 family = binomial, data = data) # plot over animal's taxa for clarity

bm_seq <- seq(0.1, 15, length.out = 200)
naive_levels <- levels(data$Naive)

grid <- expand.grid(Bmkg = bm_seq, Naive = naive_levels)
grid$ObsID <- levels(data$ObsID)[1]
grid$Seat <- levels(data$Seat)[1]
pred_link <- predict(model_naive, newdata = grid, type = "link", se.fit = TRUE, re.form = NA)
grid$pred <- plogis(pred_link$fit)
grid$pred_low <- plogis(pred_link$fit - 1.96 * pred_link$se.fit)
grid$pred_high <- plogis(pred_link$fit + 1.96 * pred_link$se.fit)
grid$Naive <- factor(grid$Naive, levels = c("Naif", "Informe"), labels = c("Uninformed", "Informed"))
grid <- arrange(grid, Naive, Bmkg)

probs$Naive <- factor(sub(".*\\.", "", as.character(probs$Sp)), levels = c("Naif", "Informe"), labels = c("Uninformed", "Informed"))
grid <- arrange(grid, Naive, Bmkg)

pdf("LR_plot.pdf")
ggplot(grid, aes(x = Bmkg, y = pred, linetype = Naive, shape = Naive, group = Naive)) +
 geom_line(linewidth = 1, color = "black") +
 geom_ribbon(aes(ymin = pred_low, ymax = pred_high), fill = "grey80", alpha = 0.2, color = NA) +
 geom_point(data = probs, aes(x = BodyMass, y = Score, shape = Naive), inherit.aes = FALSE, size = 2.5, color = "black") +
 geom_errorbar(data = probs, aes(x = BodyMass, ymin = Score - SE, ymax = Score + SE), inherit.aes = FALSE, width = 0.1, color = "black") +
  scale_x_log10(name = "Body mass (kg)", breaks = c(0.1, 0.3, 1, 3, 10), labels = c("0.1", "0.3", "1", "3", "10")) +
  scale_y_continuous(name = "Detection probability", limits = c(0, 1)) +
 scale_shape_manual(values = c(15, 17), name = "Observer status") +
 scale_linetype_manual(values = c("solid", "dotted"), name = "Observer status") +
 theme_classic(base_size = 14) +
 theme(
 legend.position = "right",
 axis.text.x = element_text(angle = 45, hjust = 1)
 )
dev.off()


pdf("graph_abs.pdf")
grid$Naive <- factor(grid$Naive, levels = c("Informed", "Uninformed"), labels = c("Actively searching", "Casually driving"))

ggplot(grid, aes(x = Bmkg, y = pred, linetype = Naive, shape = Naive, group = Naive)) +
  geom_line(linewidth = 1, color = "black") +
  geom_ribbon(aes(ymin = pred_low, ymax = pred_high), fill = "grey80", alpha = 0.2, color = NA) +
  scale_x_log10(name = "Body mass (kg)", breaks = c(0.1, 0.3, 1, 3, 10), labels = c("0.1", "0.3", "1", "3", "10")) +
  scale_y_continuous(name = "Detection probability", limits = c(0, 1)) +
  scale_shape_manual(values = c(15, 17), name = "Observer status") +
  scale_linetype_manual(values = c("solid", "dotted"), name = "Observer status") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()
## Collision rates estimates ---------

# Faune-AuRA reports for 2022 (not included for GDPR compliance)
reports <- c(607, 267, 63, 335)

# 1000 bootstrap replicates of T_R, in hours
load("TR_bootstrap.RData")
quantile(TR_boot, c(0.025, 0.5, 0.975))

species_data <- data.frame(
  Species = c("Red fox", "Buzzard", "Wood pigeon", "Red squirrel"),
  Bmkg = c(4.3, 1.012, 0.49, 0.321),
  Taxa = factor(c("Mamm", "Ois", "Ois", "Mamm"), levels = levels(data$Taxa)),
  Naive = factor("Naif", levels = levels(data$Naive)),
  Seat = factor("C", levels = levels(data$Seat))
)

model_formula <- as.formula("as.factor(Detect) ~ log(Bmkg) * Taxa + log(Bmkg) * Naive + Seat + (1 | ObsID)")
bootstrap_detect <- function() {
  boot_data <- NULL
  for (g in levels(data$Naive)) {
    ids <- sample(unique(as.character(data$ObsID[data$Naive == g])), replace = TRUE)
    for (i in seq_along(ids)) {
      d <- data[data$ObsID == ids[i], ]
      d$ObsID <- factor(paste0(ids[i], "_", i))
      boot_data <- rbind(boot_data, d)
    }
  }
  fit <- tryCatch(glmer(model_formula, family = binomial, data = boot_data), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, nrow(species_data)))
  predict(fit, newdata = species_data, type = "response", re.form = NA)
}

set.seed(2026)
n_draws <- 1000
p_boot <- t(replicate(n_draws, bootstrap_detect()))
colnames(p_boot) <- species_data$Species
sum(is.na(p_boot[, 1]))

TR_draws <- sample(TR_boot, n_draws, replace = length(TR_boot) < n_draws)
lambda_boot <- 24 * sweep(1 / p_boot, 2, reports, "*") / TR_draws

summary_table <- data.frame(
  Species = species_data$Species,
  Reports = reports,
  p = apply(p_boot, 2, median, na.rm = TRUE),
  p_low = apply(p_boot, 2, quantile, 0.025, na.rm = TRUE),
  p_high = apply(p_boot, 2, quantile, 0.975, na.rm = TRUE),
  lambda = apply(lambda_boot, 2, median, na.rm = TRUE),
  lambda_low = apply(lambda_boot, 2, quantile, 0.025, na.rm = TRUE),
  lambda_high = apply(lambda_boot, 2, quantile, 0.975, na.rm = TRUE)
)
summary_table$pct_reported <- 100 * summary_table$Reports / summary_table$lambda
summary_table[, -1] <- round(summary_table[, -1], 3)
summary_table


