# devtools::install_github("bnicenboim/eeguana")

setwd("C:/Users/coshk/Desktop/Perm_project/EEG")
require(eeguana)
require(httr)
require(dplyr)
require(ggplot2)

# Pre-processed EEG data from BrainVision 2.0.
# The data belong to a simple experiment where a participant was presented 100 faces
# and 100 assorted images in random order. The task of the experiment was to
# mentally count the number of faces.

#GET("https://osf.io/q6b7x//?action=download",
#    write_disk("./faces.vhdr", overwrite = TRUE),
#    progress()
#)
#GET("https://osf.io/ft5ge//?action=download",
#    write_disk("./faces.vmrk", overwrite = TRUE),
#    progress()
#)
#GET("https://osf.io/85dgj//?action=download",
#    write_disk("./faces.dat", overwrite = TRUE),
#    progress()
#)

# The file faces.vhdr contains the metadata and links to the other two files,
# faces.vmrk contains the triggers and other events in the samples,
# and faces.dat contains the signals at every sample for every channel recorded.


# The function read_vhdr() creates a list with data frames for the signal, events,
# segments information, and incorporates in its attributes generic EEG information.
faces <- read_vhdr("faces.vhdr")
names(faces) # .signal, .events, .segments
# faces$.signal contains 525207 samples from 32 channels + EOGV and EOGH

# Some intervals were marked as "bad" by BrainVision, and so we'll remove them from
# the data. We'll also segment and baseline the data. In this experiment,
# the trigger "s70" was used for faces and "s71" for no faces. We'll segment the data
# using these two triggers (200 segments found).
faces_segs <- faces %>%
  eeg_segment(.description %in% c("s70", "s71"),
              lim = c(-.2, .25)
  ) %>%
  eeg_events_to_NA(.type == "Bad Interval") %>%
  eeg_baseline()

names(faces_segs)

# Once the eeg_lst is segmented, the segments table includes the relevant columns
# from the events table (but without the leading dots).
segments_tbl(faces_segs)

# We modify the entire object
faces_segs_some <- faces_segs %>%
  mutate(
    condition =
      if_else(description == "s70", "faces", "non-faces")
  ) %>%
  select(-type)


# ggplot() applied to an eeg_lst object will downsample the signals (when needed),
# and convert them to a long-format data frame that is feed into ggplot.
# This object can then be customized.
faces_segs_some %>%
  select(O1, O2, P7, P8) %>%
  ggplot(aes(x = .time, y = .value)) +
  geom_line(alpha = .1, aes(group = .id, color = condition)) +
  stat_summary(
    fun.y = "mean", geom = "line", alpha = 1, size = 1.5,
    aes(color = condition)
  ) +
  facet_wrap(~.key) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = .17, linetype = "dotted") +
  theme(legend.position = "bottom")

# Another possibility is to create a topographic plot of the two conditions,
# by first making segments that include only the interval .1-.2 s after the onset
# of the stimuli, creating a table with interpolated amplitudes
# and using the ggplot wrapper plot_topo.
faces_segs_some %>%
  filter(between(as_time(.sample, unit = "milliseconds"), 100, 200)) %>%
  group_by(condition) %>%
  summarize_at(channel_names(.), mean, na.rm = TRUE) %>%
  plot_topo() +
  annotate_head() +
  geom_contour() +
  geom_text(colour = "black") +
  facet_grid(~condition)
