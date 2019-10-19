function auto_combine(taskfile, escfile, savepath, savefile,task)

HHV_COLUMN = 71;
VHV_COLUMN = 70;
REP_COLUMN = 37;
LEP_COLUMN = 8;
REPV_COLUMN = 36;
LEPV_COLUMN = 7;

st = load(taskfile);
esc = load(escfile);

hhv = esc.Data(:,HHV_COLUMN);
vhv = esc.Data(:,VHV_COLUMN);
rep = esc.Data(:,REP_COLUMN);
lep = esc.Data(:,LEP_COLUMN);
repv = esc.Data(:,REPV_COLUMN);
lepv = esc.Data(:,LEPV_COLUMN);

time= esc.Data(:,1);

%experimental data
t=table(time,hhv,vhv, rep,lep,repv, lepv, 'variablenames',...
                         {'time','hhv','vhv','rep','lep','repv','lepv'});


%add column for sample number
sample_rate = t.time(2);
t.sampletime = ceil(t.time/sample_rate);

%convert esc date format to number of seconds since midnight
datatime=split(esc.Date,':');
get_hour = split(datatime(1)," ");
datatime=str2num(datatime{2})*60+str2num(datatime{3})+360*str2num(get_hour{2}); %just take seconds
stimtime=360*st.systemtime(4)+st.systemtime(5)*60+st.systemtime(6);

file_diff = datatime- stimtime;

sttable = table(st.xposition, ...
                st.yposition, ...
                st.approxtime(1:end-1), ...
                'variablenames', ...
                {'raw_targ', ...
                 'raw_vtarg',...
                 'approxtime'});
        
             
if task == 'AS'
    sttable.raw_targ = sttable.raw_targ+640;
end
        
sttable.adjustedtime = sttable.approxtime - file_diff;
sttable.sampletime = ceil(sttable.adjustedtime/sample_rate);


tt= outerjoin(t,sttable,'key','sampletime','MergeKeys',true);

writetable(tt, [savepath savefile])