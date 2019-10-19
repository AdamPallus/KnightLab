function batch_combiner()
    % Read Spreadsheet
     disp('CHOOSE MASTER SHEET')
    [file, path] = uigetfile({'*.xlsx;*.csv', 'data sheets'},'Load Master Spreadsheet');

    %
    try
        disp('Master Sheet Loaded!')
        df = readtable([path, file]);

    catch
        disp('Failed to load sheet')
        return
    end
    %
    disp('CHOOSE DATA FOLDER')
    path = uigetdir('Identify Data Folder');
    display(['LOOKING IN' ,path])
    %
    disp('CHOOSE OUTPUT FOLDER')
    outputpath = uigetdir('Choose Save Folder');

    for row = 1:height(df)
        behaviorfilename = df.BehaviorFile{row};
        stimfile = df.StimulusFile{row};
        task = df.Task{row};
        id = df.ID{row};

        date_append = behaviorfilename(1:8); % should be yyyymmdd
        savefilename = [id,'-',task,'-',date_append,'.csv'];

        if isfile([outputpath,'/',savefilename])
            disp(['Skipping ',savefilename, ': already combined'])
            continue
        end

        escfilefull = [path, ['/',behaviorfilename, '.mat']];
        if ~isfile(escfilefull)
            disp('Checking alternative name')
            escfilefull = [path, ['/',stimfile]];
            if ~isfile(escfilefull)
                disp(['Cannot find ',escfilefull])
                disp(['Aborting ', savefilename])
                continue
            else
               disp('********************* Used alternative name*********') 
               disp(escfilefull)
            end
        end

        taskfilefull = [path,'/',task,stimfile];
        if ~isfile(taskfilefull)
            %check with '-'
            disp('Checking alternative name')
            taskfilefull = [path,'/',task,behaviorfilename,'.mat'];
            if ~isfile(taskfilefull)
                disp(['Cannot find ',taskfilefull])
                disp(['Aborting ', savefilename])
                continue
            else
                disp('********************* Used alternative name*********')
                disp(taskfilefull)
            end
        end
        %try
            auto_combine(taskfilefull, escfilefull, [outputpath,'/'], savefilename, task)
        %    disp(['****SUCCESS*** Saved ', savefilename,' Successfully'])
        %catch
        %    disp(['Failed to combine ', savefilename])
        end
    end