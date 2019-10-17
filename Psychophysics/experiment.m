% function experiment(subject, conditionNr)
% 
% Run the experiment described in "Optimal inference of sameness" by 
% Van den Berg, Vogel, Josic, and Ma, PNAS 2012.
%
% CONDITION 0: PRACTICE
% CONDITION 1: SIGMA_NZ THRESHOLDS
% CONDITION 2: ELLIPSE ECCENTRICITY THRESHOLDS
% CONDITION 3: ACTUAL EXPERIMENT, LOW EPSILON        (EXPERIMENT 2 in paper)
% CONDITION 4: ACTUAL EXPERIMENT, MIXED EPSILON      (EXPERIMENT 2 in paper)
% CONDITION 5: ACTUAL EXPERIMENT, HIGH EPSILON       (EXPERIMENT 2 in paper)
% CONDITION 6: ACTUAL EXPERIMENT, HIGH EPSILON, N=6  (EXPERIMENT 1 in paper)
% CONDITION 7: ACTUAL EXPERIMENT, HIGH EPSILON, N=4  (EXPERIMENT 1 in paper)
% CONDITION 8: ACTUAL EXPERIMENT, HIGH EPSILON, N=2  (EXPERIMENT 1 in paper)
% CONDITION 9: COLOR EXPERIMENT                      (EXPERIMENT 3 in paper)
% CONDITION 10: COLOR EXPERIMENT PRACTICE
% CONDITION 11: DEMO

% Written by RvdB, 2010

function experiment(subject, conditionNr)

% init randomizer
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(s);

%-%-%-%-%-
%- INIT %-
%-%-%-%-%-
subject = upper(subject);
makeScreenShot  = 0;     % if 1, then screenshots of stimuli will be made
screen_width    = 40;    % in cm (Dell@T115A: ~48cm; Dell@T101C: ~40 cm)

% load settings
settings = LoadExperimentSettings(subject,conditionNr);

% file names etc
lutfile    = 'BCM_BOOTHC.txt';  % file with 8-bit lookup table
if strcmp(settings.exptype,'exp3_color')
    outputfile = ['output_color/' upper(subject) '_condition_' num2str(conditionNr) '_' datestr(clock,30) '.mat'];
else
    outputfile = ['output/' upper(subject) '_condition_' num2str(conditionNr) '_' datestr(clock,30) '.mat'];
end

% screen info
screenNumber=max(Screen('Screens'));       % use external screen if exists
[w h]=Screen('WindowSize', screenNumber);  % screen resolution
screen_resolution = [w h];                 % screen resolution
screen_distance = 60;                      % distance between observer and screen (in cm)
screen_angle = 2*(180/pi)*(atan((screen_width/2) / screen_distance)) ; % total visual angle of screen
screen_ppd = screen_resolution(1) / screen_angle;  % pixels per degree
screen_fixposxy = screen_resolution .* [.5 .5]; % fixation position

% read calibration file
[DAC R G B W] = textread(lutfile,'%d%f%f%f%f','headerlines',1);
lut = [R G B W];

% get colors
colMat = getColMat;

% open screen
gray=GrayIndex(screenNumber);
windowPtr = screen('OpenWindow',screenNumber,gray,[],32,2);
HideCursor;

% show start screen
screen('TextSize',windowPtr,30);
screen('DrawText',windowPtr,'Press any key to start',20,screen_resolution(2)/2,[255 255 255]);
screen('Flip', windowPtr);
waitForKey;
screen('FillRect', windowPtr, gray);
screen('Flip', windowPtr);
screen('TextSize', windowPtr, 15);

sameKey = kbName('s');
diffKey = kbName('d');
escKey  = kbName('esc');

%-%-%-%-%-%-%-%-%-%-%-%-%
%- LOOP THROUGH TRIALS %-
%-%-%-%-%-%-%-%-%-%-%-%-%
sessionStartTime = getSecs;
nConsecutiveCorrect = 0;
fixsize = 4;
fixcol = round(interp1(W,DAC,settings.bglum))*1.2;
nextFlipTime = 0; % just to initialize...
aborted = 0;
trialnr = 0;
infoTrials = round(settings.noTrials .* [.5]);
if settings.noTrials<50
    infoTrials = 1000;
end
nCorrect = 0;

while (trialnr < settings.noTrials) && ~aborted
    trialnr = trialnr + 1;
    
    % Randomly choose an N
    N = settings.N(ceil(rand*length(settings.N)));
    
    % Randomly choose a sigma_NZ
    sigmaNz = settings.sigmaNZ(ceil(rand*length(settings.sigmaNZ)));
    
    % Randomly assign epsilons
    if strcmp(settings.exptype,'actual_mixed')
        nLow = ceil(rand*(settings.N-1));
        setEpsilons = [ones(1,nLow)*settings.epsilons(1) ones(1,settings.N-nLow)*settings.epsilons(2)];
        setEpsilons = setEpsilons(randperm(settings.N));
    end
    
    for ii=1:N
        setEpsilons(ii) = settings.epsilons(ceil(rand*length(settings.epsilons)));
    end
    if ~(settings.multiEps)
        setEpsilons(2:end) = setEpsilons(1);
    end
    
    % Assign orientations (or colors)
    setMu = rand*180;
    if (rand>settings.pSame)
        setSigma = sigmaNz;
        setOrts = normrnd(setMu,setSigma,1,N);
    else
        setSigma = 0;
        setOrts = ones(1,N)*setMu;
    end    
    
    % Set some remaining variables
    BGDAC = round(interp1(W,DAC,settings.bglum));
    settings.bglum;
    
    % create stimulus patches
    clear stimtex;
    for ii=1:N
        fgdac = round(interp1(lut(:,4),0:255,settings.fglum));
        if strcmp(settings.exptype,'exp3_color')
            setOrts(setOrts<0.5)=setOrts(setOrts<0.5)+180;
            setOrts(setOrts>180.5)=setOrts(setOrts>180.5)-180;
            radius = round(sqrt(settings.stimsize/pi)*screen_ppd);
            template = fspecial('disk',radius);
            for jj=1:3                
                template(template~=0)=colMat(round(setOrts(ii)),jj);
                im(:,:,jj)=template;               
            end
            im(im==0)=BGDAC;
        else
            im = drawEllipse((settings.stimsize * screen_ppd)^2,setEpsilons(ii),setOrts(ii),fgdac,BGDAC);
        end
        patchsize(ii,1) = size(im,2);
        patchsize(ii,2) = size(im,1);
        stimtex(ii)=screen('MakeTexture', windowPtr, im);
    end
    
    % Compute stimulus positions
    if (conditionNr>5)
        dAngle = pi/4;
    else
        dAngle = pi/3;
    end
    setPosThetas(1) = rand*2*pi;
    setPosThetas(2:N) = mod(setPosThetas(1) + dAngle .* ((2:N)-1),2*pi);
    setPosRho(1:N) = screen_ppd * settings.ecc;
    for ii=1:N
        [x y] = pol2cart(setPosThetas(ii),setPosRho(ii));
        stimPos(ii,1) = round(.5*screen_resolution(1)+x);
        stimPos(ii,2) = round(.5*screen_resolution(2)+y);
    end
    
    % jitter
    stimPos = normrnd(stimPos,settings.posJitterSigma*screen_ppd);
    
    % SCREEN 1: FIXATION
    screen('fillRect',windowPtr,BGDAC);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize,1);
    %         screen('DrawText',windowPtr,['Trial # = ' num2str(trialnr)],20,0,round(BGDAC*1.2));
    screen('flip',windowPtr,nextFlipTime);
    nextFlipTime = getSecs + .5;
    
    % SCREEN 2: STIMULUS
    screen('fillRect',windowPtr,BGDAC);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize,1);
    for ii=1:N
        srcrect = [0 0 patchsize(ii,:)];
        destrect = centerRectOnPoint(srcrect,stimPos(ii,1),stimPos(ii,2));
        screen('drawtexture',windowPtr,stimtex(ii),srcrect,destrect);
    end
    screen('flip',windowPtr,nextFlipTime);
    stimStartTime = getSecs;
    nextFlipTime = getSecs + settings.stimtime/1000;
    
    if (makeScreenShot)
        grabSize = 2.5 * screen_ppd * settings.ecc;
        grabrect = centerRectOnPoint([0 0 grabSize grabSize],screen_fixposxy(1),screen_fixposxy(2));
        ss = screen('getimage',windowPtr,grabrect);
        imwrite(ss,['screenshots/stim_' num2str(trialnr) '.png'],'png');
    end
    
    % SCREEN 3: RESPONSE
    screen('fillRect',windowPtr,BGDAC);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize,1);
    screen('flip',windowPtr,nextFlipTime);
    REALSTIMTIME = round(1000*(getSecs - stimStartTime));
    
    responseStartTime = getSecs;
    done=0;
    while ~done
        keyCode = waitForKey;
        if (keyCode==sameKey)
            done=1;
            RESP='same';
            RESPIDX=0;
        elseif (keyCode==diffKey)
            done=1;
            RESP='diff';
            RESPIDX=1;
        elseif (keyCode == escKey)
            aborted=1;
            break;
        end
        
    end
    
    CORRECT = (RESPIDX==0 && setSigma==0) || (RESPIDX==1 && setSigma>0);
    nCorrect = nCorrect + CORRECT;
    
    % in case of training, adjust stimulus time based on performance
    if strcmp(settings.exptype,'train')
        if ~CORRECT
            settings.stimtime = min(settings.stimtime*1.3, settings.maxstimtime);
            fprintf('ERROR   --> stimulus time set to %d ms\n',round(settings.stimtime));
        else
            nConsecutiveCorrect = nConsecutiveCorrect + 1;
            if (nConsecutiveCorrect==3)
                settings.stimtime = max(settings.stimtime/1.3, settings.minstimtime);
                nConsecutiveCorrect = 0;
                fprintf('CORRECT --> stimulus time set to %d ms\n',round(settings.stimtime));
            else
                fprintf('CORRECT\n');
            end
        end
    end
    
    RT = round(1000*(getSecs - responseStartTime));
    
    % store data
    realSetSigma = std(setOrts);
    setEpsilons(N+1:8)=0;
    setPosThetas(N+1:8)=0;
    setOrts(N+1:8)=999;
    data(trialnr,:) = [setMu setSigma realSetSigma RESPIDX CORRECT RT REALSTIMTIME setEpsilons setOrts setPosThetas];
    save(outputfile,'settings','data');
    
    % SCREEN 4: INTER TRIAL DISPLAY
    screen('fillRect',windowPtr,BGDAC);
    if settings.feedback
        if (CORRECT)
            drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),[0 200 0],fixsize*1,0);
        else
            drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),[200 0 0],fixsize*1,0);
        end
    else
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize,0);
    end
    screen('flip',windowPtr);
    nextFlipTime = getSecs + (settings.ITT/1000);
    
    % SCREEN 5: Some intermediate information
    if any(infoTrials==trialnr)
        % show start screen
        screen('TextSize',windowPtr,14);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize,1);
        textx = 450;
        texty = screen_resolution(2)/2 + 30;
        dy = 25;
        pFinished = round(100*trialnr/settings.noTrials);
        screen('DrawText',windowPtr,['You completed ' num2str(pFinished) '% of this block'],textx,texty,[180 180 180]); texty = texty+dy;
        if abs(pFinished-50)<5 && settings.noTrials>50
            screen('DrawText',windowPtr,['Performance so far: ' num2str(100*nCorrect/trialnr,2) '%'],textx+75,texty,[255 255 255]); texty = texty+dy;
            if (nCorrect/trialnr) > .5
                if (nCorrect/trialnr)<.7
                    screen('DrawText',windowPtr,['Good job!'],textx+130,texty,[80 180 80]); texty = texty+dy;
                elseif (nCorrect/trialnr)<.80
                    screen('DrawText',windowPtr,['Excellent!'],textx+130,texty,[80 180 80]); texty = texty+dy;
                elseif (nCorrect/trialnr)<.90
                    screen('DrawText',windowPtr,['Excellent^2!'],textx+130,texty,[80 180 80]); texty = texty+dy;
                else
                    screen('DrawText',windowPtr,['Excellent^3'],textx+130,texty,[80 180 80]); texty = texty+dy;
                end
            end
        end
        screen('Flip', windowPtr);
        waitForKey;
        screen('TextSize',windowPtr,14);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize,1);
        screen('Flip', windowPtr);
        WaitSecs(.5);
    end
    
    
end

delete('output_color/DEMO*.mat')
delete('output_color/PRACTICE*.mat')
if settings.noTrials>50
    showhiscores(windowPtr,subject,conditionNr)
end

%-%-%-%-%-%-%-
%- FINALIZE %-
%-%-%-%-%-%-%-
noSecs = round(getSecs - sessionStartTime);
noMin = floor(noSecs / 60);
noSec = mod(noSecs,60);
if (aborted)
    fprintf('\n\nSESSION ABORTED (total time = %dm%ds)\n\n',noMin,noSec);
    delete(outputfile);
else
    fprintf('\n\nSESSION COMPLETED (total time = %dm%ds)\n\n',noMin,noSec);
end
screen('closeall');
ShowCursor;


%-%-%-%-%-%-%-%-%-%-%-%-%- HELPER FUNCTIONS %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-

function keyCode = waitForKey
keyCode = ones(1,256);
while sum(keyCode(1:254))>0
    [keyIsDown,secs,keyCode] = KbCheck;
end
while sum(keyCode(1:254))==0
    [keyIsDown,secs,keyCode] = KbCheck;
end
keyCode = min(find(keyCode==1));

function drawfixation(windowPtr,x,y,color,size,vert)
Screen('DrawLine',windowPtr,color,x-size,y,x+size,y,2);
if (vert)
    Screen('DrawLine',windowPtr,color,x,y-size,x,y+size,2);
end

function colMat = getColMat()

% This conversion is based on Bruce Lindbloom's excellent website about color spaces
% http://brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html

% define stuff
ref_white  = [0.95047 1 1.08883];   % D65 reference white (this is a point in XYZ-space)
CIE_center = [58 12 13];   % (L,a,b) center point in CIE Lab space
CIE_radius = 50; % radius of the circle along which we'll pick our colors
epsilon = 216/24389; % some junction point in some function; see CIE documentation 
kappa = 24389/27; % another variable in the same function; see CIE documentation

% compute the (L,a,b)-values of the colorwheel 
L = CIE_center(1) * ones(1,180);
a = CIE_center(2) + CIE_radius*cos((1:180)*2*pi/180);
b = CIE_center(3) + CIE_radius*sin((1:180)*2*pi/180);

%-%-%-%-%-%-%-%-%-%-%-
% Step 1: Lab to XYZ %
%-%-%-%-%-%-%-%-%-%-%-
f_y = (L+16)/116;
f_x = a/500 + f_y;
f_z = f_y - b/200;

x_r                     = f_x.^3;
x_r(x_r<=epsilon)       = (116*f_x(x_r<=epsilon)-16)/kappa;
y_r(L>kappa*epsilon)    = ((L(L>kappa*epsilon)+16)/116).^3;
y_r(L<=kappa*epsilon)   = L(L<=kappa*epsilon)/kappa;
z_r                     = f_z.^3;
z_r(z_r<=epsilon)       = (116*f_z(z_r<=epsilon)-16)/kappa;

X = x_r * ref_white(1);
Y = y_r * ref_white(2);
Z = z_r * ref_white(3);

%-%-%-%-%-%-%-%-%-%-%-
% Step 2: XYZ to RGB %
%-%-%-%-%-%-%-%-%-%-%-

% From Bruce Lindbloom's website: "Given an XYZ color whose components are 
% in the nominal range [0.0, 1.0] and whose reference white is the same as 
% that of the RGB system, the% conversion to companded RGB is done in two steps."
%
% --> we need to check if our monitor's white point is D65 (that's the whitepoint that Wen-Chuang was using and that i'm using heret oo)

% STEP A: XYZ to Linear RGB
M = [ 3.2404542 -1.5371385 -0.4985314; ...   % conversion matrix for sRGB with D65 reference white 
     -0.9692660  1.8760108  0.0415560; ...   % source: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
      0.0556434 -0.2040259  1.0572252];
RGB = M*[X; Y; Z];
  
% STEP B: Companding
r = RGB(1,:); g = RGB(2,:); b = RGB(3,:);
R(r<=0.0031308) = 12.92*r(r<=0.0031308);
R(r >0.0031308) = 1.055*r(r >0.0031308).^(1/2.4);
G(g<=0.0031308) = 12.92*g(g<=0.0031308);
G(g >0.0031308) = 1.055*g(g >0.0031308).^(1/2.4);
B(b<=0.0031308) = 12.92*b(b<=0.0031308);
B(b >0.0031308) = 1.055*b(b >0.0031308).^(1/2.4);

colMat = round(255*[R; G; B]');
