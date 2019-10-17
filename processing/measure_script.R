#' 2019 full pipeline
#' 1) Use matlab function to put .csv files in a folder for analysis
#' 2) execute the script below
source('processing/loadfunctions.R')
source('processing/knighthelperfunctions.R')

files_location = 'data/'
samplerate<- 304.7508/1000

h <- loadGazeFiles(path = files_location)

#CHECK TO SEE IF RAW DATA LOOK APPROPRIATE (HEAD END GAZE MOVING IN SAME DIRECTION AS TARGET?)

# manipulate(ggplot(h %>% filter(block == 'cg01ST1',
#                                time > timestart, 
#                                time<timestart+1000))+
#              geom_line(aes(time, E), color = 'hotpink')+
#              geom_line(aes(time, H), color = 'darkblue')+
#              geom_line(aes(time, na.spline(G)), color = 'darkgreen')+
#              geom_point(aes(time, Targ)),
#            timestart = slider(1, 200000, step = 900 ))
# 


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
                                       driftcorrect = TRUE,markFixations = FALSE),
         headmovement=markMovementsDouble(Hv,threshold1=10,threshold2=4)) %>%
  do(markTagetMovements(t=.,buffer=200,threshold=200,trial.length=500))%>%
  filter(!is.na(trialnum))->
  h

#CHECK AGAIN
# manipulate({
#   trials = unique(h$trialnum)
#   ggplot(h %>% filter(block == 'CP48ST1',
#                       trialnum == trials[chosenTrial],))+
#              geom_line(aes(time, G-H), color = 'hotpink')+
#              geom_line(aes(time, H), color = 'darkblue')+
#              geom_line(aes(time, G), color = 'darkgreen')+
#              geom_point(aes(time, Targ))},
#            chosenTrial = slider(1, length(unique(h$trialnum))))

#'Next: measure trials

h %>%
  filter(!is.na(trialnum)) %>%
  group_by(task,subject,block,trialnum) %>%
  do(measureTrial(.))->
  hm


# #CHECK
# hh <- left_join(h, hm)
# 
# manipulate({
#   trials = unique(hh$trialnum)
#   ggplot(hh %>% filter(block == 'cg01ST1',
#                       trialnum == trials[chosenTrial]) %>%
#            mutate(counter = counter*samplerate))+
#     geom_line(aes(counter, G-H), color = 'hotpink')+
#     geom_line(aes(counter, H), color = 'darkblue')+
#     geom_line(aes(counter, G), color = 'darkgreen')+
#     geom_point(aes(counter, Targ))+
#     geom_vline(aes(xintercept = gaze.onset*samplerate), color = 'darkgreen')+
#     geom_vline(aes(xintercept = gaze.offset*samplerate), linetype = 2, color = 'darkgreen')+
#     geom_vline(aes(xintercept = total.gaze.offset*samplerate),linetype = 2, color = 'darkgreen')+
#     geom_vline(aes(xintercept = head.onset*samplerate),color = 'darkblue')+
#     geom_vline(aes(xintercept = head.offset*samplerate), linetype = 2,color = 'darkblue')+
#     xlab('Time (ms)')
#     
#   },
#   chosenTrial = slider(1, length(unique(hh$trialnum))))

## append to dashboard data
# hm %>%
#   ungroup() %>%
#   mutate(subject = paste0(subject, 'n'))->
#   hm


# hm_old <- readRDS('dashboard2018-02-10.RDS')
# 
# hm_old %>% 
#   filter(subject != 'CP48n') %>%
#   bind_rows(hm)->
#   hmfull
output_file = paste0('dashboards/output/dashboard',Sys.Date(),'.RDS')
saveRDS(hm, output_file)






