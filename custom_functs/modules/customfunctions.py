import numpy as np
import pandas as pd
import calim


def concat_hdf_files(list_files):
    list_files = list_files
    data_hdf = calim.Project()
    print(f"Concatenate {len(list_files)} files:")
    for filename in list_files:
        print(filename)
        data_file_hdf = calim.Project()
        data_file_hdf.from_hdf(filename) 
        for recording in data_file_hdf.recordings:
            data_hdf.append(data_file_hdf.recordings[recording])
    print("Done.")
    return data_hdf


def getActivity(start, end, c, dt):
    num_events  = len(list(c.get_event(range(start, end)))) # Number of events
    if num_events  == 0:
        activity = 0
    if num_events != 0:
        activity = (num_events/(end-start))/dt # Events/s
    return activity, num_events


def getBurstParams(list_bursts, dt=0.1, cutoffBurst=5, start=0, end=0):
    if len(list_bursts) == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    
    if len(list_bursts) != 0:        
        intraburst_freq = list_bursts.burst_freq.mean()
        #if segment begins or end with a burst, drop that respective burst
        if list_bursts.iloc[0].first_frame < start+(cutoffBurst/dt):
            list_bursts = list_bursts.drop(0)
        if len(list_bursts) != 0:
            if list_bursts.iloc[-1].last_frame > end-(cutoffBurst/dt):
                list_bursts = list_bursts[:-1]
        
        #print(list_bursts)

        burst_num_spikes = list_bursts.burst_num_spikes.mean()
        burst_length = list_bursts.burst_length.mean()
        burst_n = list_bursts.burst_n.mean()
        burst_per_cell_per_sec = list_bursts.burst_n.mean()/((end-start)*dt)
        
        return intraburst_freq, burst_num_spikes, burst_length, burst_n, burst_per_cell_per_sec
        
def getBursts(eventList, dt=0.1, maxEventLength=5, minEventInCluster=10):
        
    maxEventLength = maxEventLength/dt # Convert the maxEventLength into number of frames
    burstFrequency = pd.DataFrame([])
    list_ISI = pd.DataFrame([])
    
    if len(eventList) <= 0:
        list_ISI = pd.concat([list_ISI, pd.DataFrame(data={"frame":[np.nan], "pos": "---"})])
        burstFrequency = pd.DataFrame(data={"burst_freq": [np.nan], \
                                        "burst_num_spikes":[np.nan],\
                                        "burst_length":[np.nan],\
                                        "burst_n":[np.nan],\
                                       "first_frame":[np.nan],\
                                      "last_frame":[np.nan]})
        
        #print(burstFrequency)
        #print(list_ISI)
        return burstFrequency, list_ISI   

    eventList = np.array([e.frame for e in eventList if e.use])
    
    last = eventList[0]
    cluster = np.array([])
    n_bursts = 0

    for i in eventList[1:]:
        if i - last <= maxEventLength:
            if len(cluster) == 0:
                cluster = np.append(cluster, last)
            cluster = np.append(cluster,i)
            #print(i)
        else:
            if len(cluster)>= minEventInCluster:
                #print("cluster ende", cluster)
                meanFreq = 1 / (np.mean(np.diff(cluster))*dt)
                lenEvents = len(cluster)
                length = (cluster[-1] - cluster[0])*dt
                n_bursts += 1
                d = pd.DataFrame(data={"burst_freq": [meanFreq], \
                                       "burst_num_spikes":[lenEvents], "burst_length":[length], "burst_n":[n_bursts],\
                                       "first_frame":[cluster[0]],\
                                      "last_frame":[cluster[-1]]})
                list_ISI = pd.concat([list_ISI, pd.DataFrame(data={"frame":[cluster[0]], "pos": "first"})])
                list_ISI = pd.concat([list_ISI, pd.DataFrame(data={"frame":[cluster[-1]], "pos": "last"})])
                
                burstFrequency = pd.concat([burstFrequency, d])

                #print(meanFreq, lenEvents, length, n_bursts)

            cluster = np.array([])
        last = i
        
    if len(cluster) >= minEventInCluster:
        meanFreq = 1 / (np.mean(np.diff(cluster))*dt)
        lenEvents = len(cluster)
        length = (cluster[-1] - cluster[0])*dt
        n_bursts += 1
        d = pd.DataFrame(data={"burst_freq": [meanFreq], \
                               "burst_num_spikes":[lenEvents], "burst_length":[length], "burst_n":[n_bursts],\
                               "first_frame":[cluster[0]],\
                              "last_frame":[cluster[-1]]})
        
        list_ISI = pd.concat([list_ISI, pd.DataFrame(data={"frame":[cluster[0]], "pos": "first"})])
        list_ISI = pd.concat([list_ISI, pd.DataFrame(data={"frame":[cluster[-1]], "pos": "last"})])
        burstFrequency = pd.concat([burstFrequency, d])
        
    if len(list_ISI) != 0:
        list_ISI = list_ISI.reset_index(drop=True)
        #print(list_ISI)
        list_ISI = list_ISI.drop(0)
        list_ISI = list_ISI[:-1]
    elif len(list_ISI) == 0:
        list_ISI = pd.concat([list_ISI, pd.DataFrame(data={"frame":[np.nan], "pos": "---"})])
    burstFrequency["burst_n"] = n_bursts
    burstFrequency = burstFrequency.reset_index(drop=True)

    
    return burstFrequency, list_ISI