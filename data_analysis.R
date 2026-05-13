# set the variable "participant_of_interest" to be the id of the participant of interest (1 or 2). For example,
# participant_of_interest <- 1

# load packages
library(betareg)
library(tidyverse)

# import data from the N-of-1 trials on acne severity described in Fu et al. (https://arxiv.org/abs/2302.07547)
dat_o <- read_csv("https://raw.githubusercontent.com/HIAlab/Acne_Multimodal_Nof1/refs/heads/main/Data/scores_unscaled_combined.csv")

# define number of iterations for g-computation algorithm and parametric bootstrap
G <- 500
B <- 500


#### Data Management

# rename variables and drop unnecessary columns
dat <- dat_o %>% rename(time="Timestamp (From Photo)(MMDD-YYYY-HHMMSS)",
                        id="Id",
                        image_id="Image Id(Id-Timestamp)",
                        temperature="Temperature (°F)",
                        activity="Activity Level (Categorical), i.e none, light, intense",
                        lotion="Applied Lotion/Makeup\r\n(Boolean)",
                        treatment="Intervention\r\n(Boolean)") %>% select(-"Unnamed: 7")

# create outcome as average of scores
dat <- mutate(dat, outcome = rowMeans(select(dat, starts_with("scores_"))))

# drop scores
dat <- dat %>% select(-starts_with("scores_"))

# fix wrong timestamp for 2 data points: if the timestamp does not match the image file name, use the time in the image file name as timestamp
dat$time <- ifelse(dat$time!=substr(dat$image_id,3,18), substr(dat$image_id,3,18), dat$time)

# convert time in R format
dat$time <- lubridate::mdy_hms(dat$time)

# add a variable indicating the measurement number for each participant
dat <- dat %>% arrange(id,time) %>% group_by(id) %>% mutate(time_discrete=1:n()) %>% ungroup()

# add a variable indicating the moment of the day
dat$day_moment <- ifelse(dat$time_discrete %in% seq(1,60,by=3),"wakeup", 
                         ifelse(dat$time_discrete %in% seq(2,60,by=3),"sec_meal","bedtime")) 

# filter only the two participants of interest
dat <- dat %>% filter(id==1 | id==2)

# visualize data for participants 1 and 2
p <- ggplot(dat, aes(x=time_discrete, y=outcome, col=treatment, group=id, shape=day_moment)) + geom_point(size=3) + geom_line() +
  geom_point(data=dat,aes(x=time_discrete, y=-0.05, group=id, fill=temperature), shape=22, size=5, inherit.aes = FALSE)  +
  scale_fill_gradient("Temperature (°F)",low="blue", high="white") + 
  scale_color_manual("Treatment", values=c("#F8766D","#619CFF"))  + 
  scale_shape_manual("Moment of the day", values = c(15,16,17), labels = c("Bedtime", "Second meal", "Wake up")) +
  ylim(-0.05,1) + ylab("Outcome") + xlab("Time points") + facet_grid(vars(id))

p

# save plot
pdf("Fig4.pdf",width=8,height=6)
p
dev.off()
rm(p)

# select participant of interest (1 or 2)
dat_i <- dat %>% filter(id==participant_of_interest)


#### Data analysis 1

# test for strong conditional stationarity
betareg(outcome ~ time_discrete, data = dat_i %>% filter(treatment==TRUE)) %>% summary()
betareg(outcome ~ time_discrete, data = dat_i %>% filter(treatment==FALSE)) %>% summary()

# estimate U-CATE under basic causal model
y_a1_i <- dat_i[dat_i$treatment==TRUE,"outcome"] %>% pull()
y_a0_i <- dat_i[dat_i$treatment==FALSE,"outcome"] %>% pull() 
nt1 <- length(y_a1_i)
nt0 <- length(y_a0_i)
se <- sqrt(var(y_a1_i)/nt1 + var(y_a0_i)/nt0)

# mean difference
tau_i <- mean(y_a1_i) - mean(y_a0_i)
tau_i

# 95%CI
c( tau_i - qnorm(1-0.05/2)*se, tau_i - qnorm(0.05/2)*se )

# remove created objects
rm(y_a1_i, y_a0_i, nt1, nt0, se, tau_i)


#### Data analysis 2

set.seed(12345*participant_of_interest)

# training data with relevant variables (also lagged)
dat_train <- dat_i %>% select(id, time_discrete, day_moment, temperature, treatment, outcome)
dat_train <- bind_cols(dat_train %>% add_row(.before=1) %>% filter(row_number() < n()) %>% setNames(paste0('lag1_', names(.))), dat_train)

# beta regression to predict the outcome
model_outcome_o <- betareg(outcome ~ treatment + temperature + day_moment + lag1_treatment + lag1_outcome, data=dat_train)

# ordinal logistic regression to predict the covariate temperature from one day to the next
model_temperature_o <- rms::lrm(temperature ~ lag1_temperature, data=dat_train %>% filter(day_moment=="wakeup"), penalty=0.5)

# create a function to initialize a dataset that will be filled with counterfactual outcomes (same structure as dat_train)
initialize_counterfactual_dataset <- function(treatment_strategy) {
  temp_d <- dat_train %>% mutate_at(.vars=vars(lag1_day_moment,lag1_temperature,lag1_treatment,lag1_outcome,day_moment,temperature,treatment,outcome), .funs=function(x) x=NA)
  temp_d$treatment <- treatment_strategy
  # fix first time point
  temp_d[2,c("lag1_day_moment", "lag1_temperature", "lag1_treatment", "lag1_outcome")] <- temp_d[1,c("day_moment", "temperature", "treatment", "outcome")] <- dat_train[1,c("day_moment", "temperature", "treatment", "outcome")]
  return(temp_d)
}

# create a function to fill the initialized dataset with counterfactual outcomes based on the beta and ordinal logistic models
fill_counterfactual_dataset <- function(temp_d, mod_out, mod_temp) {
  for (i in 2:nrow(temp_d)) {
    temp_d[i,"day_moment"] <- ifelse(pull(temp_d[i,"lag1_day_moment"])=="bedtime", "wakeup", ifelse(pull(temp_d[i,"lag1_day_moment"])=="wakeup", "sec_meal", "bedtime"))
    if (pull(temp_d[i,"lag1_day_moment"])=="bedtime") {
      probs <- as.numeric(predict(mod_temp, newdata=temp_d[i,], type="fitted.ind"))
      if (length(mod_temp$yunique)==2) {probs <- c(1-probs,probs)}
      temp_d[i, "temperature"] <- sample(mod_temp$yunique, size = 1, replace=T, prob = probs)
      rm(probs)
      }
    if (pull(temp_d[i,"lag1_day_moment"])!="bedtime") {temp_d[i, "temperature"] <- temp_d[i, "lag1_temperature"]}
    shape1 <- predict(mod_out, newdata=temp_d[i,], type="response")*as.numeric(mod_out$coefficients$precision)
    shape2 <- as.numeric(mod_out$coefficients$precision) - shape1
    temp_d[i, "outcome"] <- rbeta(1, shape1=shape1, shape2=shape2)
    
    if (i+1 <= nrow(temp_d)) {temp_d[i+1,c("lag1_treatment","lag1_day_moment", "lag1_temperature", "lag1_outcome")] <- temp_d[i,c("treatment","day_moment","temperature","outcome")]}
    rm(shape1, shape2)
  }
  return(temp_d)
}

# create a matrix to store results from the g-computation algorithm
point_effect <- matrix(NA_real_, nrow = G, ncol = nrow(dat_train))

# g-computation algorithm
for (g in 1:G) {
  # obtain outcome under "never treatment" strategy
  dat_counterfactual <- initialize_counterfactual_dataset(treatment_strategy = FALSE)
  dat_counterfactual <- fill_counterfactual_dataset(dat_counterfactual, model_outcome_o, model_temperature_o)
  count_outcomes_0 <- dat_counterfactual$outcome
  
  # obtain outcome under "always treatment" strategy
  dat_counterfactual <- initialize_counterfactual_dataset(treatment_strategy = TRUE)
  dat_counterfactual <- fill_counterfactual_dataset(dat_counterfactual, model_outcome_o, model_temperature_o)
  count_outcomes_1 <- dat_counterfactual$outcome
  
  # take difference of the outcomes between the two scenarios
  point_effect[g,] <- count_outcomes_1-count_outcomes_0
  rm(dat_counterfactual, count_outcomes_0, count_outcomes_1)
}

# take the mean across G repetitions, this is the point estimate
point_effect_o <- apply(point_effect, MARGIN=2, FUN=mean)
rm(point_effect)

# confidence interval via parametric bootstrap
effect <- matrix(NA_real_, nrow = B, ncol = nrow(dat_train))

for (b in 1:B) {
  # generate data from the models fitted in the original sample under the assigned treatment schedule
  dat_boot <- initialize_counterfactual_dataset(treatment_strategy = dat_train$treatment)
  dat_boot <- fill_counterfactual_dataset(dat_boot, model_outcome_o, model_temperature_o)
  
  # repeat the procedure to estimate the point estimate in the bootstrapped sample
  model_outcome <- betareg(outcome ~ treatment + temperature + day_moment + lag1_treatment + lag1_outcome, data=dat_boot)
  model_temperature <- rms::lrm(temperature ~ lag1_temperature, data=dat_boot %>% filter(day_moment=="wakeup"), penalty=0.5)
  
  point_effect <- matrix(NA_real_, nrow = G, ncol = nrow(dat_train))
  
  for (g in 1:G) {
    dat_counterfactual <- initialize_counterfactual_dataset(treatment_strategy = FALSE)
    dat_counterfactual <- fill_counterfactual_dataset(dat_counterfactual, model_outcome, model_temperature)
    count_outcomes_0 <- dat_counterfactual$outcome
    
    dat_counterfactual <- initialize_counterfactual_dataset(treatment_strategy = TRUE)
    dat_counterfactual <- fill_counterfactual_dataset(dat_counterfactual, model_outcome, model_temperature)
    count_outcomes_1 <- dat_counterfactual$outcome
    
    point_effect[g,] <- count_outcomes_1-count_outcomes_0
    rm(dat_counterfactual, count_outcomes_0, count_outcomes_1)
  }
  
  # save point estimate obtained in the bootstrapped sample
  effect[b,] <- apply(point_effect, MARGIN=2, FUN=mean)
  rm(model_outcome, model_temperature, point_effect)
}

# save point estimates and results from the parametric bootstrap
save(point_effect_o, file=paste0("pointestimates_participant",participant_of_interest,".RData"))
save(effect, file=paste0("bootstrapruns_participant",participant_of_interest,".RData"))

# build 95% confidence intervals
report <- data.frame(point=point_effect_o, std=apply(effect, MARGIN=2, FUN=sd))
report$lower <- report$point -qnorm(1-0.05/2)*report$std
report$upper <- report$point -qnorm(0.05/2)*report$std
report$time_discrete <- 1:nrow(report)

# visualize estimates
p <- ggplot(report %>% filter(time_discrete!=1), aes(x=time_discrete, y=point)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3) + xlab("Time point") + ylab("Estimated individual-specific causal effect")
p

# save plot
pdf(paste0("Fig5",ifelse(participant_of_interest==1,"a","b"),".pdf"),width=8,height=6)
p
dev.off()
rm(p)
