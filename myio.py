import pickle
import numpy as np
import scipy
import scipy.io
import bioread
import mne

def save_pkl(data, out_name):
    fout = open(out_name, 'wb')
    pickle.dump(data, fout)
    fout.close()

    return

def load_pkl(in_name):
    fin = open(in_name, 'rb')
    data = pickle.load(fin)
    fin.close()

    return data

def load_biopac(fname, target_channel):

    if '.acq' in str(fname):
        data = my_load_acq(fname)

    elif '.mat' in str(fname):
        data = scipy.io.loadmat(fname)
        data['sampling_rate'] = np.squeeze(data['sampling_rate'])
        data['labels'] = [i.rstrip() for i in data['labels']]
        data['events_tind'] = np.squeeze(data['events_tind'])
        data['events_tsec'] = np.squeeze(data['events_tsec'])

    target_channel_ind = data['labels'].index(target_channel)
    data['data'] = data['data'][target_channel_ind, :]
    data['labels'] = data['labels'][target_channel_ind]
    data['units'] = data['units'][target_channel_ind]
    data['sampling_rate'] = data['sampling_rate'][target_channel_ind]

    return data

def load_brainvision(file, target_channel):
    raw = mne.io.read_raw_brainvision(file)
    raw.set_channel_types({'ECG': 'ecg', 'HEOG': 'eog', 'VEOG': 'eog','PPG':'ecg'}, verbose=None)
    print(raw.info)

    raw.pick_channels(target_channel)
    data = dict()
    data['data'] = raw.get_data()
    data['sampling_rate'] = raw.info['sfreq']
    #data['labels'] = [i.rstrip() for i in data['labels']]
    #data['events_tind'] = np.squeeze(data['events_tind'])
    #data['events_tsec'] = np.squeeze(data['events_tsec'])

    return data

def my_load_acq(fname):

    data = bioread.read_file(fname)
    out_dict = dict()
    names, units, sf = [], [], []

    signals = np.stack([data.channels[i].data for i in range(len(data.channels))])
    for i in range(len(data.channels)):
        names.append(data.channels[i].name.rstrip())
        units.append(data.channels[i].units)
        sf.append(data.channels[i].samples_per_second)

    out_dict['data'] = signals
    out_dict['labels'] = names
    out_dict['units'] = units
    out_dict['sampling_rate'] = sf
    return out_dict

