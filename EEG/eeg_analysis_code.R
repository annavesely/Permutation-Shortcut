#devtools::install_github("bnicenboim/eeguana")

# http://headit.ucsd.edu/studies/3316f70e-35ff-11e3-a2a9-0050563f2612
# Description: Subjects listened to voice recordings that suggest
# an emotional feeling and ask subjects to imagine an emotional scenario
# or to recall an experience in which they have felt that emotion before.


setwd("C:/Users/coshk/Desktop/Perm_project/EEG")
library(devtools)
library(magrittr)
library(eeguana)
library(plyr)
library(dplyr)
library(ggplot2)


load("data_seg_median_sample.Rdata")
names(data_seg) # .signal, .events, .segments
dim(data_seg$.signal) # 32832 x 258: .id, .sample, 256 channels
dim(data_seg$.events) # 0 x 4: .id, .initial, .final, .channel
dim(data_seg$.segments) # 64 x 4: .id, condition, .recording, subj

names(data_seg$.signal)
table(data_seg$.signal$.id)
data_seg$.signal$.sample

str(data_seg$.segments)
str(data_seg$.events)

table(data_seg$.segments$condition) # 32 neg, 32 pos
table(data_seg$.segments$subj) # 2 segments per subject

# str(data_seg)



#Plot of all the Event-Related Potential of one electrode
data_seg %>%
   select(D25) %>%
   ggplot(aes(x = .time, y = .value)) +
   geom_line()

#Plot ERP of each subject (average across condition):
data_seg %>%
   select(D25) %>%
   ggplot(aes(x = .time, y = .value)) +
   geom_line(aes(group = subj))  +
   stat_summary(
     fun = "mean", geom = "line", alpha = 1, size = 1.5,
     aes(color = "red"),show.legend = FALSE
   )

# data_seg$.signal$

# By-sampl means of the channels
mean_ch = chs_mean(data_seg,na.rm = TRUE)  
str(mean_ch)
mean_ch


#Plot ERP of condition (average across subject)
mean_ch %>%
  filter(between(as_time(.sample, unit = "s"), .0, .5)) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = .1, aes(group = .id, color = condition)) +
  stat_summary(
    fun = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = condition)
  ) +
  facet_wrap(~.key) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = .17, linetype = "dotted") +
  theme(legend.position = "bottom")


ERP_data <-  data_seg %>%
  group_by(.sample, condition) %>%
  summarize_at(channel_names(.), mean, na.rm = TRUE)

ERP_plot<- ERP_data %>%
  # select(.sample,.id,D2,D5,D6,D10,D11,D13) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(aes(color = condition)) +
  facet_wrap(~.key) +
  theme(legend.position = "bottom") +
  ggtitle("ERPs for disgust vs object") +
  theme_eeguana()
ERP_plot


ERP_plot %>% plot_in_layout()

data_seg %>%
  filter(between(as_time(.sample, unit = "s"), .05, .15)) %>%
  group_by(condition) %>%
  summarize_at(channel_names(.), mean, na.rm = TRUE) %>%
  plot_topo() +
  annotate_head() +
  geom_contour() +
  geom_text(colour = "black") +
  facet_grid(~condition)

############################### TO GET THE DATA TO TEST:

range_s=c(.05, .15)
file_name="data_seg_median_sample.Rdata"
load(file_name)

D=data_seg%>%filter(between(as_time(.sample, unit = "s"), range_s[1], range_s[2])) %>%
  group_by(subj,condition)%>%summarize_all(mean, na.rm = TRUE)

names(D)
dim(D$.signal) # 64 x 258
D$.signal$.id
D$.segments
table(D$.segments$subj,D$.segments$condition)






D %>%
  select(A1)


Dnew <- D %>%
  group_by(condition) %>%
  #group_by(channel_names(.)) %>%
  mutate(
    tA1 = D$.signal$A1[1] - D$.signal$A1[2]
  ) %>%
  select(tA1)


t_comp <- function(X, ch){
  
}




cn <- channel_names(D)

Dnew <- D %>%
  group_by(subj) %>%
  (function(x) (D$.signal$x[1] - D$.signal$x[2]))
  


ddply(D$.signal, channel_names(D), function(x) {
  mean.count <- mean($count)
  sd.count <- sd(x$count)
  cv <- sd.count/mean.count
  data.frame(cv.count = cv)
  })




names(D$.signal)

Dnew <- D %>%
  group_by(condition) %>%
  summarize_at(channel_names(.), difference, na.rm = TRUE)





