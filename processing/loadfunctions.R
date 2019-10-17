
#1/9/2019
#how to use:
#Use the .m file to combine the two mat files and put them in a directory
#The files should be named XXXXST1 where XXXX is the subject ID, ST is the task type
#and the final number is the block number

#1) use h<-loadGazeFiles on the folder where the .csv files are
#2) use hm<- measureTrials(h) to measure the trials
#3) the restulting data frame can be used for plotting or sent to the dashboard

get_pixel_conversion<- function(){
  pixels <- c(640,704,768,832,896,960,1024,1088,1152,1216,
             576,512,448,384,320,256,192,128,64)
  
  deg <- c(0,10.81,20.58,29.3,37,43.6,49.2,53.7,57.2,59.6,
           -10.5,-20,-28.5,-36,-42.6,-48.1,-52.7,-56.2,-58.8)
  
  d <- data.frame(Targ =pixels,
                  deg = deg)
  
 mod <- lm(deg~poly(Targ,3),data=d)

}


fixH<- function(t){
  #fixes head drift when converting velocity to position
  mod<- lm(H~poly(time,3),data=t)
  t %>%
    mutate(H=H-predict(mod,newdata=t))->
    t
  return(t)
}

loadfromMAT<- function(filepath = '~/knight/ConcussionGaze/fromtom/testoutput/',
                       filename = 'AP17ST1.csv' ){
  
  samplerate<- 304.7508/1000
  
  h = read.csv(paste0(filepath, filename),na.strings = "NaN")
  
  names<-str_match(filename,"(^[a-zA-Z0-9]{4})([a-zA-Z0-9]{3})")
  
  #I don't know why there are NA values here in the time column, but na.rm works
  target_start = min(h$time[!is.na(h$raw_targ)],na.rm=T)
  
  h <- filter(h, sampletime>= target_start)

  h %>%
    fill(raw_targ) %>%
    select(sampletime,rep,raw_targ,hhv) %>%
    rename(time = sampletime,
           E = rep,
           HV = hhv,
           Targ = raw_targ) %>%
    mutate(HV = replace(HV, abs(HV) > 400, 0),
           HV = replace(HV, is.na(HV), 0),
           HV = HV,
           HV = HV*-1, #flip head left/right
           E = E*-1, #flip eye left/right
           E=atan(E)*180/pi*2,
           H=cumsum(HV)/samplerate/1000,
           block = names[1],
           blocknum = as.numeric(str_sub(names[3],3)),
           task = str_sub(names[3],1,2),
           subject = names[2])->
    h
  


   h %>%
     group_by(block) %>%
     #mutate(H_old = H) %>%
     do(fixH(.)) %>%
     mutate(G=H+E)->
     h
  # h
  

}


loadGazeFiles<- function(referencefile=NULL,path='~/knight/ConcussionGaze/fromtom/testoutput/'){
  require(stringr)
  require(dplyr)
  require(data.table)
  require(tidyr)

  #get names of all files in path
  files <- list.files(path=path,pattern='*.csv')
  nfiles<-length(files)
  
  if (nfiles>0){
    message(c('New Files: ',files))

    loadedfiles<-lapply(files, loadfromMAT,filepath = path)

    t <-rbindlist(loadedfiles,fill = TRUE)
    
    #load model for converting pixels to degrees
    mod <- get_pixel_conversion()
    t$Targ = predict(mod,newdata=t)
    return(t)
    
  }else{
    message('********NO NEW DATA********')
    t<-NULL
  }

}


measureTrials<- function(h){
  require(dplyr)
  source('~/knight/ConcussionGaze/Scripts/knighthelperfunctions.R')
  
  h%>%
    select(G,H,Targ,block,subject,task,blocknum)%>%
    group_by(block) %>%
    filter(!is.na(G)) %>%
    mutate(time=row_number(),
           Graw=G,
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
  
  #Force antisaccade gaze to be on target...
  h %>%
    ungroup() %>%
    mutate(Targ=replace(Targ,task=='AS' &  counter<200,0))->
    h
  #measure trials join to previous table if needed
  h %>%
    filter(!is.na(trialnum)) %>%
    group_by(task,subject,block,trialnum) %>%
    do(measureTrial(.))->
    hm
  
  hm %>%
    mutate(recenter=(task=='AS' & abs(IGP)>10)) %>%
    filter(!recenter) %>%
    select(-recenter)->
    hm
  return(hm)
}
