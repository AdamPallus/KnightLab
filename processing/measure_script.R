#' 2019 full pipeline
#' 1) Use matlab function to put .csv files in a folder for analysis
#' 2) execute the script below
source('processing/loadfunctions.R')
source('processing/knighthelperfunctions.R')

files_location = 'data/'
task = 'ST'
samplerate<- 304.7508/1000

h <- loadGazeFiles(path = files_location, task = task)

#'Next: remove noise, blinks etc
h%>%
  select(G,H,Targ,block,subject,task,blocknum)%>%
  group_by(block) %>%

  mutate(time=row_number(),
         Graw=G,
         G = replace(G, is.na(G), 0), #replace na with zero since that's how it was originally
         G=replace(smooth(G,"3R"),G==0,NA), #mark missing data as NA rather than 0
         Gnospline=G,
         G= applyspline(G,6),
         target.velocity=parabolicdiff(Targ,7)*samplerate,
         Gv=parabolicdiff(G,7)*samplerate, #calculate velocity
         Hv=parabolicdiff(H,7)*samplerate,
         # gazeshifts=markMovementsDouble(Gv,threshold1=100,threshold2=10),
         gazeshifts=markSaccadesDouble(Gv,threshold1=100,threshold2=20,
                                       driftcorrect = FALSE,markFixations = FALSE),
         headmovement=markMovementsDouble(Hv,threshold1=10,threshold2=4)) %>%
  do(markTagetMovements(t=.,buffer=200,threshold=200,trial.length=500))%>%
  filter(!is.na(trialnum))->
  h

h %>%
  group_by(block, trialnum) %>%
  mutate(Targ = replace(Targ, (counter< 200 & task=='AS'), 0)) ->
  h

#add function to save the files as CSVs for opening in Spike2
#write.csv(filter(h, block == 'aj28-ST-20150818'), 'spike2test.csv') ##



#'Next: measure trials

h %>%
  filter(!is.na(trialnum)) %>%
  group_by(task,subject,block,trialnum) %>%
  do(measureTrial(.))->
  hm


h <- left_join(h, hm)
saveRDS(h, paste0('dashboards/output/',task,'/newtrialtest', Sys.Date(),'.RDS'))

output_file = paste0('dashboards/output/',task,'/newdashboard', Sys.Date(),'.RDS')
saveRDS(hm, output_file)






