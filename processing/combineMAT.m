function combineMAT()

files = {};
while length(files) ~= 2
    [files path] = uigetfile('*.mat','Choose Experiment File',...
        'multiselect', 'on');
end

if strcmp(files{1}(1:2), 'ST')
    stfile = files{1};
    escfile = flies{2};
else
    stfile = files{2};
    escfile = files{1};
end

disp(path)

% [stfile, stpath] = uigetfile('*.mat','Choose Experiment File');
% if stfile == 0
%     return
% end
% [escfile, escpath] = uigetfile('*.mat','Choose Data File');
% if escfile == 0
%     return
% end
[savefile, savepath] = uiputfile('*.csv', 'Save Combined File');
if savepath == 0
    %disp('ABORTING...')
    %return
%     ts = split(path,'\');
%     sprintf('Auto saving as: %s', ts{10})
%     sprintf('%s.csv',ts{10})
%     savepath = path;
end
% xx = split(path, '\');
% savefile = xx{10};
% savepath = path;

st = load([path stfile]);
esc = load([path escfile]);

hhv = esc.Data(:,71);
rep = esc.Data(:,37);
lep = esc.Data(:,8);
time= esc.Data(:,1);

%experimental data
t=table(time,hhv,rep,lep,'variablenames',{'time','hhv','rep','lep'});


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
                st.approxtime(1:end-1), ...
                'variablenames', ...
            {'raw_targ', 'approxtime'});
        
sttable.adjustedtime = sttable.approxtime - file_diff;
sttable.sampletime = ceil(sttable.adjustedtime/sample_rate);


tt= outerjoin(t,sttable,'key','sampletime','MergeKeys',true);

writetable(tt, [savepath savefile])