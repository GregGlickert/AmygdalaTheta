import argparse
import os
import sys
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fooof import FOOOF
from fooof.sim.gen import gen_aperiodic
from scipy.signal import decimate, welch
from scipy.signal.windows import hann as hanning

scale = 1

def raster(spikes_df,node_set,skip_ms=0,ax=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms] 
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        ax.scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])
        
        #depth_of_mod = 2
        #freq = 8
        #t = np.arange(0,15000,0.001)
        #ax.plot(100*depth_of_mod*np.sin(2 * np.pi * freq * t)+500)
    
    vpsi_input = h5py.File('./vpsi_inh_spikes.h5')


    vpsi_spikes = pd.DataFrame({'node_ids':vpsi_input['spikes']['vpsi_inh']['node_ids'],
                            'timestamps':vpsi_input['spikes']['vpsi_inh']['timestamps']})


    ax.scatter(vpsi_spikes['timestamps'],vpsi_spikes['node_ids']+1000,s=0.25,c='black',label='VPSI spikes')
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels))
    ax.grid(True)
    #ax.set_xlim(6500, 9000)

path = 'outputECP'
dt = 0.1
steps_per_ms = 1/dt
skip_seconds = 5
skip_ms = skip_seconds*1000
skip_n = int(skip_ms * steps_per_ms)
end_ms = 15000

spikes_location = os.path.join(path,'spikes.h5')

print("loading " + spikes_location)
f = h5py.File(spikes_location)
spikes_df = pd.DataFrame({'node_ids':f['spikes']['BLA']['node_ids'],'timestamps':f['spikes']['BLA']['timestamps']}) 
print("done")

node_set = [
    {"name":"PN","start":0*scale,"end":799*scale,"color":"blue"},
    {"name":"PV","start":800*scale,"end":892*scale,"color":"red"},
    {"name":"SOM","start":893*scale,"end":943*scale,"color":"green"},
    {"name":"CR","start":944*scale,"end":999*scale,"color":"purple"}
]
fig, (ax1) = plt.subplots(1,1,figsize=(15,4.8))#6.4,4.8 default
fig.suptitle('Amygdala Theta Analysis')
raster(spikes_df,node_set,skip_ms=skip_ms,ax=ax1)
plt.show()