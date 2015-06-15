# Set path for all machines
if(.Platform$OS.type == 'windows') {
  path <- 'D:/SimonBesnard'
} else {
  info <- Sys.info()
  if (info['nodename'] == 'Chouchen') {
    path <- '/media/simonbesnard/External_SB/SimonBesnard/PhD_MPI/RO1/Inputs'
  } else if (info['nodename'] == 'papaya') {
    path <- '/media/DATA3/besna001'
  } else if (info['nodename'] == 'tanargue') {
    path <- 'media/whatever/'
  }
  
} # For some reasons the empty line at the bottom is important
