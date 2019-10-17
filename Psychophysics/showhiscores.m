function showhiscores(windowPtr,subjid,condnr)

try

    % read scores
    if condnr>=9
        files = dir(['output_color/*_' num2str(condnr) '_*.mat']);
    else
        files = dir(['output/*_' num2str(condnr) '_*.mat']);
    end
    subjlist.perf = [];
    for ii=1:length(files)
        allsubjid{ii} = files(ii).name(1:strfind(files(ii).name,'_')-1);
        if condnr>=9
            load(['output_color/' files(ii).name]);
        else
            load(['output/' files(ii).name]);
        end        
        allscores(ii) = 100*sum(data(:,5))/length(data(:,5));
        perf(ii) = 100*sum(data(:,5))/length(data(:,5));
        if strcmp(subjid,allsubjid{ii})
            subjlist.perf(end+1) = perf(ii);
        end
    end
    subjlist.lastperf = subjlist.perf(end);
    subjlist.perf = sort(subjlist.perf,'descend');
    subjlist.lastidx = min(find(subjlist.perf==subjlist.lastperf));
    
    % build group hi score list
    uID = unique(allsubjid);
    for ii=1:length(uID)
        groupscores(ii)=0;
        for jj=1:length(allsubjid)
            if strcmp(allsubjid{jj},uID{ii}) && allscores(jj)>groupscores(ii)
                groupscores(ii)=allscores(jj);
            end
        end
    end
    [groupscores I] = sort(groupscores,'descend');
    grouplist.perf = groupscores;
    grouplist.subjid = uID(I);
    
    % screen info
    screenNumber=max(Screen('Screens'));       % use external screen if exists
    [w h]=Screen('WindowSize', screenNumber);  % screen resolution
    screen_resolution = [w h];                 % screen resolution
    windowOpened = (windowPtr < 0);
    
    % read header
    headerim = imread('hof.png');    
    headerimsize = size(headerim(:,:,1));
    headerimsize = headerimsize([2 1]);
    headerimrect = [0 0 headerimsize];  
    headerimrect([1 3]) = headerimrect([1 3]) + round((screen_resolution(1) - headerimsize(2))/2) - 50;
        
    if windowPtr < 0
        windowPtr = screen('OpenWindow',screenNumber,[128 128 128],[],32,2);
    end
    
    % show hiscore lists
    screen('fillRect',windowPtr,[128 128 128]);
    screen('PutImage',windowPtr,headerim,headerimrect);

    screen('TextSize',windowPtr,20);
    screen('DrawText',windowPtr,['CONDITION #' num2str(condnr)],530,340,[80 200 80]); 
    
    % Personal scores
    currx = 660;
    curry = 400;
    dy = 20;
    screen('DrawText',windowPtr,'PERSONAL TOP 20',currx,curry,[80 80 180]); curry=curry+2*dy;
    screen('TextSize',windowPtr,16);
    for ii=1:min(length(subjlist.perf),20)
        perfstr = num2str(subjlist.perf(ii),3);
        if isempty(strfind(perfstr,'.'))
            perfstr = [perfstr '.0'];
        end
        screen('DrawText',windowPtr,[num2str(ii) '. '],currx,curry,[255 255 255]);
        screen('DrawText',windowPtr,[perfstr '%'],currx+40,curry,[255 255 255]);
        if (ii==subjlist.lastidx)
            screen('DrawText',windowPtr,'<-- most recent',currx+120,curry,[180 180 80]);
        end
        curry=curry+dy;
    end    
    if windowOpened
        screen('DrawText',windowPtr,['avg = ' num2str(mean(subjlist.perf),3) '%'],currx,curry+dy,[180 180 80]);            
    end
    
    % Group scores
    currx = 260;
    curry = 400;
    dy = 20;
    screen('TextSize',windowPtr,20);
    screen('DrawText',windowPtr,'GROUP HISCORE LIST',currx,curry,[80 80 180]); curry=curry+2*dy;
    screen('TextSize',windowPtr,16);
    for ii=1:length(grouplist.perf)      
        perfstr = num2str(grouplist.perf(ii),3);
        if isempty(strfind(perfstr,'.'))
            perfstr = [perfstr '.0'];
        end
        screen('DrawText',windowPtr,[num2str(ii) '. '],currx,curry,[255 255 255]);
        if strcmp(grouplist.subjid{ii},subjid)
            screen('DrawText',windowPtr,grouplist.subjid{ii},currx+40,curry,[255 255 0]);
        else
            screen('DrawText',windowPtr,grouplist.subjid{ii},currx+40,curry,[255 255 255]);
        end
        screen('DrawText',windowPtr,[' ' perfstr '%'],currx+80,curry,[255 255 255]);
        curry=curry+dy;
    end    

    
    % show "press any key"
    screen('DrawText',windowPtr,'Press any key to continue',20,screen_resolution(2)*.9,[255 255 255]);
    screen('Flip',windowPtr);
    
    waitForKey;
  
    if (windowOpened)
        screen('closeall');
    end
    
catch    
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end


function keyCode = waitForKey
keyCode = ones(1,256);
while sum(keyCode(1:254))>0
    [keyIsDown,secs,keyCode] = KbCheck;
end
while sum(keyCode(1:254))==0
    [keyIsDown,secs,keyCode] = KbCheck;
end
keyCode = min(find(keyCode==1));
