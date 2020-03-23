# devtools::install_github("bnicenboim/eeguana")



require(eeguana)
require(httr)

# Pre-processed EEG data from BrainVision 2.0.
# The data belong to a simple experiment where a participant was presented 100 faces
# and 100 assorted images in random order. The task of the experiment was to
# mentally count the number of faces.

GET("https://osf.io/q6b7x//?action=download",
    write_disk("./faces.vhdr", overwrite = TRUE),
    progress()
)
GET("https://osf.io/ft5ge//?action=download",
    write_disk("./faces.vmrk", overwrite = TRUE),
    progress()
)
GET("https://osf.io/85dgj//?action=download",
    write_disk("./faces.dat", overwrite = TRUE),
    progress()
)

# The file faces.vhdr contains the metadata and links to the other two files,
# faces.vmrk contains the triggers and other events in the samples,
# and faces.dat contains the signals at every sample for every channel recorded.


# The function read_vhdr() creates a list with data frames for the signal, events,
# segments information, and incorporates in its attributes generic EEG information.
faces <- read_vhdr("faces.vhdr")

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

#34






