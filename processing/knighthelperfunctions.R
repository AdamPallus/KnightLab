parabolicdiff <- function(pos,n=7){
  require(cladoRcpp)
  #Differentiates by convolving the position trace with an n-point parabola
  q <- sum(2*((1:n)^2))
  convoutput<- rcpp_convolve(pos,c(-n:-1, 1:n))
  convoutput<- convoutput[(n*2):(length(pos)-((n*2)+1))]
  vels<- c(array(convoutput[1],dim=n*2),convoutput,array(convoutput[length(convoutput)],dim=n*2))
  vels <- vels/q*-1000
}

maxabs<- function(x){
  m1<-max(x,na.rm=T)
  m2<-min(x,na.rm=T)
  if (abs(m1)>abs(m2)) {
    return(m1)
  } else{
    return(m2)
  }
}

markTagetMovements<-function(t,buffer=200,threshold=1000,trial.length=1000){
  #This function breaks the data files into chunks based on movements of the target
  #The function allows for overlap such that if you group_by(trialnum) each group will be the same size 
  #The code here is reused from saccade marking functions which is why it references saccades in variable names
  findSaccades<-function(ev,threshold=40){
    mindur<- 1
    i<-which(abs(ev)>threshold) #find all the times when speed > threshold
    sacoff<-which(diff(i)>mindur) #minimum duration of an accepted saccade
    sacon<-c(1,sacoff+1) #first saccade
    sacoff<-c(sacoff,length(i)) #end of last saccade
    saccade.onset<-i[sacon] #get actual times
    saccade.offset<-i[sacoff] 
    
    return(data.frame(saccade.onset,saccade.offset))
  }
  
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 2: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$trialnum<- stimes[[4]]
    # df$saccade.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes<-findSaccades(t$target.velocity,threshold)
  stimes %>%
    mutate(dur=saccade.offset-saccade.onset,
           s=row_number(),
           saccade.onset=saccade.onset-buffer,
           saccade.offset=saccade.onset+trial.length+2*buffer)->
    stimes
  
  x<- as.data.frame(rbindlist(apply(stimes,1,jsac)))
  x %>%
    group_by(trialnum) %>%
    mutate(counter=time-first(time)) ->
    x
  t<-left_join(t,x,by='time')
  select(t,-target.velocity)
}

loadnewheadfree<- function(referencefile=NULL,path="~/GitHub/ConcussionGaze/kdata/"){
  require(stringr)
  require(dplyr)
  require(data.table)
  #This function loads .csv files in a particular folder. They must have the same columns for rbind
  #Saves time by only reading the csv when necessary
  
  #get names of all files in path
  files <- list.files(path=path,pattern='*.txt')
  #extract neuron name eg. Bee-01
  # names<-sapply(files, str_match,"^[a-zA-Z]+-[0-9]+",USE.NAMES=FALSE)
  names<-sapply(files, str_match,"^[a-zA-Z0-9]+",USE.NAMES=FALSE)
  # check for new cells
  if (!is.null(referencefile)){
    files<-files[!names %in% referencefile$neuron] #comparison
  }
  
  nfiles<-length(files)
  
  if (nfiles>0){
    message(c('New Files: ',files))
    loadedfiles <- lapply(paste(path,files,sep=''),read.csv,sep='\t',header=FALSE)
    t<-data.frame()
    temp<- NULL
    # t.old<-NULL
    for (i in 1:nfiles) {
      f<- files[i]
      message(paste('Loading:',f))
      # temp[[i]]=loadedfiles[[i]]
      # names<-str_match(f,"(^[a-zA-Z]+)-([0-9]+)")
      # loadedfiles[[i]]$block<-names[1]
      # loadedfiles[[i]]$subject<-names[2]
      # loadedfiles[[i]]$blocknum<-as.numeric(names[3])
      #updated for 2-25-2018 
      #remove extra rows
      # if (ncol(loadedfiles[[i]])>8){
      #   loadedfiles[[i]]<- select(loadedfiles[[i]],1:8)
      # }
      #files are named with Subject ID bU07 and block ST1 with ST standing for Saccade Task
      names<-str_match(f,"(^[a-zA-Z0-9]{4})([a-zA-Z0-9]{3})")
      loadedfiles[[i]]$block<-names[1]
      loadedfiles[[i]]$subject<-names[2]
      loadedfiles[[i]]$blocknum<-as.numeric(str_sub(names[3],3))
      task<-str_sub(names[3],1,2)
      loadedfiles[[i]]$task<- task
      if (task=='AS' & ncol(loadedfiles[[i]])>8){ #antisaccade dask needs target data fixed
        #'There are three targets in the anti saccade files
        #'but only the second is useful. The central fixation is always at zero
        #'If this changes, the value in V7 should be used instead of zero
        targ<- loadedfiles[[i]]$V8
        targ<- replace(targ,targ< -80,0)
        loadedfiles[[i]]$V8=targ
      }
      
    }
    t <-rbindlist(loadedfiles,fill = TRUE)
    # t<- dplyr::select(t, -thp,-tvp,-time)
  }else{
    message('********NO NEW DATA********')
    t<-NULL
  }
  if (ncol(t) == 12){
    names(t)<- c('G','GV','H','HV','E','EV','sampletime','Targ','block','subject','blocknum','task')
  }
  
  if (ncol(t)==14){
    t$V9<-NULL
  }
  
  
  if (ncol(t) == 13){
    names(t)<- c('G','GV','H','HV','E','EV','sampletime','Targ','block','subject','blocknum','task','time')
  }
  
  
  if (ncol(t) == 7){
    names(t)<-c('HV','E','Targ','block','subject','blocknum','task')
  }
  t %>%
    mutate(time=row_number())->
    t
  return(t)
}

removeBlinks<-function(x,buffer=5){
  #This function finds periods of missing data
  #and removes additional points nearby, defined by buffer
  
  #This is necessary beacuse there is a period where data are not missing, 
  #but they are not accurate. 
  
  i<-which(is.na(x)) #find all the times when the data are missing
  sacoff<-which(diff(i)>1)
  sacon<-c(1,sacoff+1) #first saccade
  sacoff<-c(sacoff,length(i)) #end of last saccade
  missing.onset<-i[sacon] #get actual times
  missing.offset<-i[sacoff] 
  
  blinks<- data.frame(onset=missing.onset,offset=missing.offset)
  
  xnew<- x
  
  blinks %>%
    mutate(expanded.onset=onset-buffer,
           expanded.offset=offset+buffer,
           expanded.onset=replace(expanded.onset,expanded.onset<0,0),
           expanded.offset=replace(expanded.offset,expanded.offset>length(x),length(x)))->
    blinks
  
  for (k in 1:nrow(blinks)){
    xnew[blinks$expanded.onset[k]:blinks$expanded.offset[k]]<- NA
  }  
  
  xnew
}

applyspline<-function(G,buffer=5){
  #This function is for interpolating the periods of missing data
  
  xout<- which(is.na(G))
  # xout<-which(is.na(hh$G))
  if (length(xout)>0){
    interp<- removeBlinks(G,buffer)
    interp<-splinefun(1:length(G),interp,method="monoH.FC")(1:length(G))
    # interp<-splinefun(hh$time,interp,method="natural")(hh$time)
    interp
    # hh<-mutate(hh,interp=replace(G,xout,spline(time,G,xout=xout,method='fmm')$y))
  }else{
    G
  }
}

markMovementsDouble<- function(v, threshold1=60,threshold2=20,min.dur=3,mark.fixations=FALSE){
  #This function is an R implementation of a two-threshold event marker
  #The algorithm works like this: Find all the times when velocity is above the high threshold
  #Extend this out until velocity is below the lower threshold
  
  #in practice, I'm identifying all the times that the saccades cross the low threshold and then rejecting 
  #any that don't meet the higher threshold
  #I'm also rejecting events that are below a certain duration
  
  #this algorithm also assigns positive ID numbers to the saccades and 
  #negative ID numbers to the non-saccades (fixations?)
  #after running this function, you can group_by(event) and measure the fixations or saccades as you wish
  
  require(dplyr)
  require(data.table) #for rbindlist - a fast version of do.call('rbind') that uses data.table
  datalength<-length(v)
  i<-which(abs(v)>threshold2) #find all the times when speed is above the lower threshold
  
  #For continuous numbers, diff=1. If there is a larger jump, it means there's a gap. 
  #That indicates another saccade
  sacoff<-which(diff(i)>1) #sacoff now contains the indices of the ends of all the saccades
  #sacon has all the indices of when the saccades start
  #After an offset, the next index will be the onset of the next saccade.
  #find the onsets by looking at the next index after the offset
  sacon<-c(1,sacoff+1) #first saccade always starts at first index
  sacoff<-c(sacoff,length(i)) #end of last saccade is always at the end
  event.onset<-i[sacon] #Convert from the indices to actual times
  event.offset<-i[sacoff] 
  
  #event.onset now has the time (in sample time) of all saccade onsets
  #set up stimes as a data.frame with two columns. Onset and offset. 
  stimes<- data.frame(event.onset,event.offset) 
  
  #this is a little function that works with the weirdness of R's "apply" family of functions
  #it just takes the onset and offset and returns the whole range. 
  #if you give it [10, 20], it will return [10 11 12 13 14 15 16 17 18 19 20]
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 4: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$event<- stimes[[4]]
    df$event.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes %>%
    mutate(dur=event.offset-event.onset, #calculate duration of each saccade
           s=row_number())-> #assign an ID number to each saccade
    stimes
  
  #Use "apply" to run the function "jsac" (above) on each line of the "stimes" data.frame
  #the result is the times of when all the saccades happen
  x<-rbindlist(apply(stimes,1,jsac)) 
  
  v<- data.frame(v=v) #Make the original velocity trace into a data.frame
  v<- mutate(v, time=row_number()) #add time to keep track
  
  #join the marked saccades and the velocity
  #the result is the velocity trace plus a row that just indicates whether there's a saccade
  #each saccade is identified by it's unique marker "event" that comes from df$event<- stimes[[4]] above 
  xx<- left_join(v,x,by='time')
  
  xx %>%
    group_by(event) %>% #This means we analyze each saccade individually
    summarize(max.vel=max(abs(v)), #calculate max velocity
              dur=n()) %>% #calculate duration
    filter(max.vel>threshold1, #reject all saccades that fail to exceed the large threshold
           dur>min.dur)-> #reject all saccades that fail to exceed the minimum duration
    xm #xm is a summary which means it just lists the saccades and their measured values
  
  xx %>% #go back to the full data set and now reject all the saccades that were rejected above
    filter(event %in% unique(xm$event)) %>% 
    dplyr::select(time,event) -> #All we need is the time and the eventID
    g
  
  
  if (nrow(g)==0){
    return(1:datalength*NA)
  } else{
    
    g<- left_join(v,g,by='time')
    return(g$event)
  }
}

measureTrial<- function(tt, buffer=200){
  #This function will receive data from one trial and return a one-row data frame
  
  trial.length<- nrow(tt)
  
  tt$E = tt$G-tt$H
  # task<-first(tt$task)
  
  #identify the id number of the first saccade that satisfies the criteria 
  firstshift=min(tt$gazeshifts[buffer:trial.length],na.rm=T)
  
  if (is.infinite(firstshift)){
    #Don't bother measuring this trial if there isn't a primary gaze shift
    # message(paste0('No primary gaze shift in: ',tt$trialnum[1]))
    return(data.frame())
  }
  
  #In this block, we will reject any saccades with a delay of more than 150ms since the last saccade
  #we will also reject saccades where the peak velocity is less than 200 ms as we assume they are corrective
  
  if (length(unique(tt$gazeshifts)) > 2){
    tt %>%
      filter(!is.na(gazeshifts)) %>%
      group_by(gazeshifts) %>%
      summarize(start.time=first(time),
                stop.time=last(time),
                peak.gaze.velocity=max(Gv)) %>%
      mutate(interval= start.time-lag(stop.time)) %>%
      filter(interval<47, peak.gaze.velocity>200) -> #reject saccades that are separated by too much time (47 samples ~ 150 ms)
      isi
    good.saccades<- c(firstshift, isi$gazeshifts)
    
    #delete them from the marking
    good.gazeshifts=replace(tt$gazeshifts,!tt$gazeshifts %in% good.saccades,NA)
  }else{
    good.saccades<- firstshift
    good.gazeshifts<- tt$gazeshifts
  }
  
  #measure gaze
  fixation.velocity=mean(tt$Gv[buffer-100:buffer])
  
  gaze.onset=max(first(tt$counter[good.gazeshifts==firstshift&!is.na(good.gazeshifts)]),1)
  
  gaze.offset=last(tt$counter[good.gazeshifts==firstshift&!is.na(good.gazeshifts)])
  gaze.dur=gaze.offset-gaze.onset
  gaze.amp=tt$G[gaze.offset]-tt$G[gaze.onset]
  
  total.gaze.offset=max(tt$counter[!is.na(good.gazeshifts)])
  total.gaze.amp=tt$G[tt$counter==total.gaze.offset]-tt$G[tt$counter==gaze.onset]
  total.gaze.dur=tt$G[total.gaze.offset]-tt$G[gaze.onset]
  peak.gaze.velocity=maxabs(tt$Gv[good.gazeshifts==firstshift])
  IHP=tt$H[tt$counter==gaze.onset]
  IGP=tt$G[tt$counter==gaze.onset]
  IEP=IGP-IHP
  IEPr = IEP * sign(gaze.amp) #relative initial eye position
  final.EP=tt$E[tt$counter==total.gaze.offset]
  head.contribution=tt$H[tt$counter==gaze.offset]-IHP
  total.head.contribution=tt$H[tt$counter==total.gaze.offset]-IHP #how much head moved during gaze shifts
  
  
  #measure head
  firsthead=min(tt$headmovement[buffer:trial.length],na.rm=T)
  #set up a flag to determine whether we will be measuring head movements
  head.onset<-NA
  head.offset<-NA
  peak.head.velocity<-NA
  head.amp<-NA

  if (!is.infinite(firsthead)){
    head.onset=first(tt$counter[tt$headmovement==firsthead&!is.na(tt$headmovement)])
    head.offset=last(tt$counter[tt$headmovement==firsthead&!is.na(tt$headmovement)])

    if (head.onset<total.gaze.offset){ #ignore head movements that start after the gaze shifts are done
      peak.head.velocity=maxabs(tt$Hv[tt$headmovement==firsthead])
      head.amp=tt$H[tt$counter==head.offset]-tt$H[tt$counter==head.onset]
    }
    if (head.onset>(gaze.onset+50)){
      head.onset=NA
      head.offset=NA
      peak.head.velocity=NA
      head.amp=NA
    }
  }
  if (is.na(head.offset)){
    post.VOR.EP = NA
  }else{
    post.VOR.EP =  tt$E[tt$counter == head.offset]
  }
  #measure target
  target.amp=tt$Targ[300]-tt$Targ[1]
  
  #assess missing
  Graw=replace(tt$Graw,tt$Graw==0,NA)
  total.missing=sum(is.na(Graw))
  missing.critical=sum(is.na(Graw[(buffer-100):(buffer+100)]))
  missing.gs=sum(is.na(Graw[gaze.onset:total.gaze.offset]))
  
  #gaze steps
  gaze.steps=length(good.saccades)
  
  #VOR
  post.saccade.drift=mean(tt$Gv[(gaze.offset+10):(gaze.offset+20)])
  post.saccade.head=mean(tt$Hv[(gaze.offset+10):(gaze.offset+20)])
  Ev=tt$Gv-tt$Hv
  post.saccade.eye=mean(Ev[(gaze.offset+10):(gaze.offset+20)])
  post.saccade.VOR=post.saccade.eye/post.saccade.head
  
  
  
  tt %>%
    ungroup() %>%
    summarize_at(c('block','trialnum','blocknum'),funs(first))->
    d

  # message(paste0('Attempting: ',tt$block[1],'-',tt$trialnum[1],'\n'))
  
  
  df<- data.frame(#firstshift,
            fixation.velocity,
            gaze.onset,
            gaze.offset,
            gaze.dur,
            gaze.amp,
            total.gaze.offset,
            total.gaze.amp,
            total.gaze.dur,
            peak.gaze.velocity,
            gaze.steps,
            IHP,
            IGP,
            IEP,
            IEPr,
            head.contribution,
            total.head.contribution,
            post.saccade.drift,
            post.saccade.head,
            post.saccade.eye,
            post.saccade.VOR,
            post.VOR.EP,
            #firsthead,
            head.onset,
            head.offset,
            peak.head.velocity,
            head.amp,
            target.amp,
            total.missing,
            missing.critical,
            missing.gs#,
            # task
  )
  
  # message(paste0('Success: ',tt$trialnum[1],'\n'))
  if (nrow(df) ==1){
    d<- bind_cols(d,df)
  }else{
    return(data.frame())
  }
  
}

AdjustCalibration<-function(h,eyegain=1,headgain=1,eyeoffset=0,headoffset=0,
                            applyfilter=TRUE,filterfreq=0.00007,skipsamples=500,
                            samplerate=0.3047508){
  require(dplyr)

  

  h %>%
    mutate(time=row_number()) %>%
    dplyr::filter(time>skipsamples) %>%
    mutate(H=cumsum(HV*headgain)/samplerate/1000+headoffset)->
    h

  
  if (applyfilter){
    # require(signal) #for butterworth
    # filter<-dplyr::filter
    # 
    # 
    # filterButter<- function(y,freqs=0.00007,type='high'){
    #   require(signal)
    #   
    #   bf <- butter(2, freqs,type=type)
    #   # bf <- butter(2, c(.0005,.006))
    #   return(filtfilt(bf, y))
    # }
    # h <-mutate(h,H=filterButter(H,freqs=filterfreq))
  }
  h %>%
    mutate(E=E*eyegain+eyeoffset,
           G=H+E)->
    h
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function
  #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

markSaccadesDouble<- function(v, threshold1=60,threshold2=20,min.dur=5,maxreject=1000,
                              driftcorrect=FALSE,markFixations=TRUE){
  #This function is an R implementation of a two-threshold event marker
  #The algorithm works like this: Find all the times when velocity is above the high threshold
  #Extend this out until velocity is below the lower threshold
  
  #in practice, I'm identifying all the times that the saccades cross the low threshold and then rejecting 
  #any that don't meet the higher threshold
  #I'm also rejecting events that are below a certain duration
  
  #this algorithm also assigns positive ID numbers to the saccades and 
  #negative ID numbers to the non-saccades (fixations?)
  #after running this function, you can group_by(event) and measure the fixations or saccades as you wish
  require(dplyr)
  require(data.table) #for rbindlist - a fast version of do.call('rbind') that uses data.table
  datalength<-length(v)
  i<-which(abs(v)>threshold2) #find all the times when speed is above the lower threshold
  
  #For continuous numbers, diff=1. If there is a larger jump, it means there's a gap. 
  #That indicates another saccade
  sacoff<-which(diff(i)>1) #sacoff now contains the indices of the ends of all the saccades
  #sacon has all the indices of when the saccades start
  #After an offset, the next index will be the onset of the next saccade.
  #find the onsets by looking at the next index after the offset
  sacon<-c(1,sacoff+1) #first saccade always starts at first index
  sacoff<-c(sacoff,length(i)) #end of last saccade is always at the end
  event.onset<-i[sacon] #Convert from the indices to actual times
  event.offset<-i[sacoff] 
  
  #event.onset now has the time (in sample time) of all saccade onsets
  #set up stimes as a data.frame with two columns. Onset and offset. 
  stimes<- data.frame(event.onset,event.offset) 
  
  #this is a little function that works with the weirdness of R's "apply" family of functions
  #it just takes the onset and offset and returns the whole range. 
  #if you give it [10, 20], it will return [10 11 12 13 14 15 16 17 18 19 20]
  jsac<- function(stimes){
    summary(stimes)
    #input should be an array of length 4: c(onsettime,offsettime, saccade.number,saccade.dur)
    df<- data.frame(time=stimes[[1]]:stimes[[2]])
    df$event<- stimes[[4]]
    df$event.dur<- stimes[[3]]
    return(df)
    # return(stimes[[1]]:stimes[[2]])
  }
  
  stimes %>%
    mutate(dur=event.offset-event.onset, #calculate duration of each saccade
           s=row_number())-> #assign an ID number to each saccade
    stimes
  
  #Use "apply" to run the function "jsac" (above) on each line of the "stimes" data.frame
  #the result is the times of when all the saccades happen
  x<-rbindlist(apply(stimes,1,jsac)) 
  
  v<- data.frame(v=v) #Make the original velocity trace into a data.frame
  v<- mutate(v, time=row_number()) #add time to keep track
  
  #join the marked saccades and the velocity
  #the result is the velocity trace plus a row that just indicates whether there's a saccade
  #each saccade is identified by it's unique marker "event" that comes from df$event<- stimes[[4]] above 
  xx<- left_join(v,x,by='time')
  
  if (driftcorrect){
    #'The way that drift correct works is that it finds when the velocity 
    #'exceeds the higher threshold and then look for the first time after that
    #'when the velocity is below 100 deg/s AND the acceleration is below 10,000
    #'This comes from Mark's published work. He determined these values by comparing
    #'the firing of MLBs and matching the duration of the burst to the duration of the saccades
    #'The reason for this correction is that monkeys with strabismus often show a post-saccadic
    #'drift where the eyes seem to re-accelerate after the end of the saccade which can cause
    #'velocity-only algorithms to overestimate the duration of the saccade. 
    xx%>%
      mutate(event2=event, #you're not allowed to change grouping variables, so make a copy first
             acc=parabolicdiff(v,7)) %>%
      group_by(event2) %>%
      mutate(realstart=first(time[v>threshold1]),
             total.end=last(time),
             driftcut=first(time[time>(realstart+10)&abs(v)<100&abs(acc)<10000]),
             rejectdrift=time>driftcut,
             event=replace(event,time>driftcut,NA)) %>%
      ungroup()->
      xx
    
    # xx %>%
    #   ungroup() %>%
    #   mutate(acc=parabolicdiff(v,7)) %>%
    #   mutate(event=replace(event,abs(v)<100 & abs(acc)<10000,NA)) ->
    #   xx
  }
  
  xx %>%
    group_by(event) %>% #This means we analyze each saccade individually
    summarize(max.vel=max(abs(v)), #calculate max velocity
              dur=n()) %>% #calculate duration
    filter(max.vel>threshold1, #reject all saccades that fail to exceed the large threshold
           max.vel<maxreject, #reject all saccades over the max threshold
           dur>min.dur)-> #reject all saccades that fail to exceed the minimum duration
    xm #xm is a summary which means it just lists the saccades and their measured values
  
  xx %>% #go back to the full data set and now reject all the saccades that were rejected above
    filter(event %in% unique(xm$event)) %>% 
    dplyr::select(time,event) -> #All we need is the time and the eventID
    g
  
  if (markFixations){
    #this next part goes through and assigns an ID to all the non-saccade portions of the data
    stimes2<- filter(stimes,s %in% unique(xm$event))
    
    ftimes<-data.frame(fix.onset=c(1,stimes2$event.offset+1),
                       fix.offset=c(stimes2$event.onset-1,datalength))
    
    ftimes %>%
      filter(fix.onset>0,fix.offset>0)%>%
      mutate(dur=fix.offset-fix.onset,
             s=row_number()) %>%
      filter(fix.onset<datalength)->
      ftimes
    
    f<-rbindlist(apply(ftimes,1,jsac))
    
    f<- select(f,-event.dur)
    
    #the code below isn't very elegant, but it just combines the fixations and saccades
    #and assigns negative IDs to the fixations and positive IDs to the saccades
    f$issaccade=FALSE 
    g$issaccade=TRUE
    fg<-rbind(f,g)
    fg<- arrange(fg,time)
    
    fg$event[!fg$issaccade]=fg$event[!fg$issaccade]*-1
    
    #This is a debugging message in case the result isn't the correct length
    #we have to return a vector of the same length as the input
    if (length(fg$event)!=datalength){
      message('FAIL')
      message(length(fg$event))
      message(datalength)
      # message(paste('FAILED: ',length(fg$event,datalength,sep='-')))
    }
    else{
      # message('SUCCESS')
    }
    
    fg$event #return just an array of the IDs of saccades and fixations
  }else{
    g<- left_join(select(xx,time),g,by='time')
    return(g$event)
  }
  
}
