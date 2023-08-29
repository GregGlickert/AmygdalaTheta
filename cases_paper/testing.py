import h5py
import matplotlib.pyplot as plt
import pandas as pd

vpsi_input = h5py.File('./vpsi_inh_spikes.h5')
spikes = vpsi_input['spikes/vpsi_inh/node_ids'][:]
time = vpsi_input['spikes/vpsi_inh/timestamps'][:]


vpsi_spikes = pd.DataFrame({'node_ids':vpsi_input['spikes']['vpsi_inh']['node_ids'],
                            'timestamps':vpsi_input['spikes']['vpsi_inh']['timestamps']})


plt.scatter(vpsi_spikes['timestamps'],vpsi_spikes['node_ids'],s=0.25)
plt.xlim((5000,9000))
plt.title('VPSI input for 4 seconds')
plt.show()