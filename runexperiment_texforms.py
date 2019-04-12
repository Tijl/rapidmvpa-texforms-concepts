#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:08:35 2018

@author: tgro5258
"""

from psychopy import core, event, visual, parallel, gui
import random,sys,json,requests,os,itertools
from glob import glob
import pandas as pd
import numpy as np

# debug things
debug_testsubject = 0
debug_usedummytriggers = 0
debug_windowedmode = 0

objectstimuli = sorted(glob('Stimuli/Orig*/*/*.png'))
texformstimuli = sorted(glob('Stimuli/Tex*/*/*.png'))
print(len(objectstimuli))
print(len(texformstimuli))
assert(len(objectstimuli)==len(texformstimuli))

if debug_testsubject:
    subjectnr = 0
else:
    # Get subject info
    subject_info = {'Subject number (update participants.tsv)':''}
    if not gui.DlgFromDict(subject_info,title='Enter subject info (update participants.tsv):').OK:
        print('User hit cancel at subject information')
        exit()
    try:
        subjectnr = int(subject_info['Subject number (update participants.tsv)'])
    except:
        raise

outfn = 'sub-%02i_task-rsvp_events.csv'%subjectnr
if not debug_testsubject and os.path.exists(outfn):
    raise Exception('%s exists'%outfn)

nstimuli = len(objectstimuli)
nsequence = 30 #number of sequences per condition (texform vs object)
nrepeats = 1 #repeats of all stimuli per sequence

random.seed(subjectnr)

refreshrate = 60
feedbackduration = .5 - .5/refreshrate
fixationduration = 1 - .5/refreshrate
stimdurations,isidurations = [],[]
for d in [.016, .03, .05, .2]:
    stimdurations.append(d - .5/refreshrate)
    isidurations.append(d - .5/refreshrate)

trigger_stimon = 1
trigger_stimoff = 2
trigger_sequencestart = 3
trigger_duration = 0.005
trigger_port = 0xcff8

webhook_url='https://hooks.slack.com/services/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'

stimnum = list(range(nstimuli))

trialstructures,sequencenumber,presentationnumber,istarget = [],[],[],[];

for i in range(nsequence):
    random.shuffle(stimnum)
    stream = []
    for j in range(nrepeats):
        x=[]
        while not x or stream and len(set(stream[-4:]).intersection(set(x[:4]))):
            x = random.sample(range(nstimuli),nstimuli)
        stream+=x
    #insert 1backs at two random pos in stream
    ntargets = random.randint(1,4)
    targetpos=[1, 1]
    t = [0 for x in range(len(stream))]
    while len(targetpos)>1 and any([y-x<20 for x,y in zip(targetpos,targetpos[1:])]):
        targetpos = sorted(random.sample(range(10,len(stream)-10),ntargets))
    for p in targetpos:
        t[p]=1
    trialstructures+=stream
    istarget+=t
    sequencenumber+=[i for x in range(len(stream))]
    presentationnumber+=list(range(len(stream)))

df = pd.DataFrame({'ustreamnumber':sequencenumber,
                   'stimnumber':trialstructures,
                   'presentationnumber':presentationnumber,
                   'istarget':istarget})

conditiontuples = list(itertools.product(range(len(stimdurations)),range(nsequence)))
random.shuffle(conditiontuples)

eventlist=pd.DataFrame()
for n,(c,s) in enumerate(conditiontuples):
    x = pd.DataFrame(df[df.ustreamnumber==s])
    x['streamnumber'] = n
    x['durationcondition'] = c
    eventlist=eventlist.append(x,1)
eventlist['duration']=[stimdurations[i]+.5/refreshrate for i in eventlist['durationcondition']]

#now double the streams for our two conditions
eventlist['stimcondition']=0
eventlist['stimpath']=[texformstimuli[i] for i in eventlist['stimnumber']]
x = eventlist.copy()
x['stimcondition']=1
x['streamnumber'] = x['streamnumber']+np.max(x['streamnumber'])+1
x['stimpath'] = [objectstimuli[i] for i in x['stimnumber']]
eventlist=eventlist.append(x,1)
eventnr=0
    
print(eventlist)
    
def writeout(eventlist):
    with open(outfn,'w') as out:
        eventlist.to_csv(out,index_label='eventnumber')

#writeout(eventlist)
#raise Exception

# =============================================================================
# %% START
# =============================================================================
try:
    if debug_windowedmode:
        win=visual.Window([700,700],units='pix',color=(-.1, -.1, -.1))
    else:
        win=visual.Window(units='pix',fullscr=True)
    mouse = event.Mouse(visible=False)

    fixation = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=10,
        name='fixation', autoLog=False)
    feedback = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=fixation.size,
        name='feedback', autoLog=False)
    progresstext = visual.TextStim(win,text='',pos=(0,100),name='progresstext')
    sequencestarttext = visual.TextStim(win,text='Press button 1 or 4 to start the sequence\nPress button 2 or 3 when the fixation dot changes colour',pos=(0,50),name='sequencestarttext')

    filesep='/'
    if sys.platform == 'win32':
        filesep='\\'
        
    objectstimtex=[]
    for (i,y) in enumerate(objectstimuli):
        objectstimtex.append(visual.ImageStim(win,y,size=256,name=y.split(filesep)[1]))
    
    texformstimtex=[]
    for (i,y) in enumerate(texformstimuli):
        texformstimtex.append(visual.ImageStim(win,y,size=256,name=y.split(filesep)[1]))
    
    targetstimtex = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=fixation.size,
        name='target', autoLog=False, color='red')
    
    def send_dummy_trigger(trigger_value):
        core.wait(trigger_duration)
            
    def send_real_trigger(trigger_value):
        trigger_port.setData(trigger_value)
        core.wait(trigger_duration)
        trigger_port.setData(0)
    
    if debug_usedummytriggers:
        sendtrigger = send_dummy_trigger
    else:
        trigger_port = parallel.ParallelPort(address=trigger_port)
        trigger_port.setData(0)
        sendtrigger = send_real_trigger

    nevents = len(eventlist)
    nsequences = eventlist['streamnumber'].iloc[-1]+1
    sequencenumber = -1
    for eventnr in range(len(eventlist)):
        first = eventlist['streamnumber'].iloc[eventnr]>sequencenumber
        if first:
            fixation.draw()
            time_stimoff=win.flip()
            if eventnr:
                eventlist.at[eventnr-1, 'time_stimoff'] = time_stimoff;
                eventlist.at[eventnr-1, 'stimdur'] = time_stimoff-eventlist.at[eventnr-1, 'time_stimon'];
            writeout(eventlist)
            sequencenumber = eventlist['streamnumber'].iloc[eventnr]
            stimcondition = eventlist['stimcondition'].iloc[eventnr]
            durationcondition = eventlist['durationcondition'].iloc[eventnr]
            stimduration = stimdurations[durationcondition]
            isiduration = isidurations[durationcondition]
            last_target = -99
            correct=0
            
            if not debug_testsubject:
                try:
                    slack_data={'text':'pp%i seq %i/%i <@tijlgrootswagers> <@amanda> <@U9C24ECQ7>'%(subjectnr,sequencenumber,nsequences),'channel':'#eeglab','username':'python'}
                    response = requests.post(webhook_url, data=json.dumps(slack_data),headers={'Content-Type': 'application/json'})
                except:
                    pass

            progresstext.text = '%i / %i'%(1+sequencenumber,1+nsequences)
            progresstext.draw()
            sequencestarttext.draw()
            fixation.draw()
            win.flip()
            k=event.waitKeys(keyList='afq', modifiers=False, timeStamped=True)
            if k[0][0]=='q':
                raise Exception('User pressed q')
            fixation.draw()
            time_fixon = win.flip()
            sendtrigger(trigger_sequencestart)
            while core.getTime() < time_fixon + fixationduration:pass
        
        response=0
        rt=0
        
        if stimcondition:
            stim = objectstimtex[eventlist['stimnumber'].iloc[eventnr]]
        else:
            stim = texformstimtex[eventlist['stimnumber'].iloc[eventnr]]
            
        stim.draw()
        fixation.draw()
        if eventlist['istarget'].iloc[eventnr]:
            targetstimtex.draw()
            
        time_stimon=win.flip()
        sendtrigger(trigger_stimon)
        
        if eventlist['istarget'].iloc[eventnr]:
            last_target=time_stimon
        
        correct=0
        k=event.getKeys(keyList='sdq', modifiers=False, timeStamped=True)
        if k:
            response=k[0][0]
            rt=k[0][1]
            if response=='q':
                raise Exception('User pressed q')
            else:
                response=1
            #correct = any(eventlist['istarget'].iloc[[eventnr-x for x in range(6)]])
            correct = rt-last_target < 1
        
        eventlist.at[eventnr, 'response'] = int(response);
        eventlist.at[eventnr, 'rt'] = rt-last_target if correct else 0;
        eventlist.at[eventnr, 'correct'] = int(correct);
        eventlist.at[eventnr, 'time_stimon'] = time_stimon;
        if not first:
            eventlist.at[eventnr-1, 'time_stimoff'] = time_stimon;    
            eventlist.at[eventnr-1, 'stimdur'] = time_stimon-eventlist.at[eventnr-1, 'time_stimon'];
        #+= [response, rt, int(correct), time_stimon, time_stimoff, time_stimoff-time_stimon]   
        while core.getTime() < time_stimon + isiduration:pass

finally:
    fixation.draw()
    time_stimoff=win.flip()
    if eventnr:
        eventlist.at[eventnr-1, 'time_stimoff'] = time_stimoff;
        eventlist.at[eventnr-1, 'stimdur'] = time_stimoff-eventlist.at[eventnr-1, 'time_stimon'];
    writeout(eventlist)
    progresstext.text = 'Experiment finished!'
    progresstext.draw()
    win.flip()
    core.wait(1)
    win.close()
    exit()


