#require(devtools)
#devtools::install_github("angeella/ARIpermutation")

setwd("C:/Users/coshk/Desktop/Perm_project/fMRI_auditory")
require(ARIpermutation)

copes <- list()
sub_ids <- sapply(c(21:40),function(x) paste0(0,x))
for (sid in 1:length(sub_ids)) {  
  copes[[sid]] <- RNifti::readNifti(system.file("extdata/AuditoryData", paste0("/sub-", sub_ids[sid] , ".nii.gz"), package = "ARIpermutation"))
  
}

img_dims <- c(91, 109, 91)
img <- array(NA, c(img_dims, length(copes)))
for(sid in (1:length(copes))){
  img[,,,sid] <- copes[[sid]]
}

scores <- matrix(img, nrow=(91*109*91), ncol=length(copes))
mask <- system.file("extdata/AuditoryData", "mask.nii.gz", package = "ARIpermutation")



B <- 500
st <- signTest(X=scores, B, alternative="two.sided", mask=mask)
pvs <- c(st$pv,st$pv_H0)
pmat <- matrix(pvs, nrow=B+1) # B+1=501 rows, 168211 columns




