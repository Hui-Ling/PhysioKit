#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import scipy.io
from pathlib import Path
import biosppy
from myio import *
from myplot import *
import argparse

def peakdet(v, delta):
    maxtab = []
    mintab = []

    #if delta <= 0:
    #    print('Input argument DELTA must be positive')

    mn, mx = np.inf, -np.inf
    mnpos, mxpos = np.nan, np.nan

    lookformax = 1

    for i in range(len(v)):
        this = v[i]
        if this > mx:
            mx, mxpos = this, i
        if this < mn:
            mn, mnpos = this, i

        if lookformax:
            if this < mx - delta:
                maxtab.append([mxpos, mx])
                mn = this
                mnpos = i
                lookformax = 0
        else:
            if this > mn + delta:
                mintab.append([mnpos, mn])
                mx, mxpos = this, i
                lookformax = 1
    return np.array(maxtab), np.array(mintab)

def event_identify(x, sampling_rate):
    peak_rise = 0.5
    delta = 60 / np.array([45, 120]) * sampling_rate
    est_hb_n_all_half = int(len(x) / (60 / 100 * sampling_rate) * 0.5) # 20210312: half of expected number of heartbeat
    # Take a first pass at peak detection:
    mxtb, mntb = peakdet(x, 1e-14)
    # set the threshold based on the 20th highest rather than the highest:
    sorted_peaks = np.sort(mxtb[:, 1])
    if len(sorted_peaks) > est_hb_n_all_half:  # 20210312: replace 20 with est_hb_n_all_half
        peak_resp = sorted_peaks[-est_hb_n_all_half]
    else:
        peak_resp = sorted_peaks[-1]


    # And now do a second pass, more robustly filtered, to get the actual peaks:
    mxtb, mntb = peakdet(x, peak_rise * peak_resp)
    events = np.zeros(x.shape)

    peaks = np.squeeze(mxtb[:, 0])
    dpeaks = np.diff(peaks)

    for ind, v in enumerate(dpeaks):
        if v > delta[0]:
            seg = x[int(peaks[ind]):int(peaks[ind + 1])]
            mxtb, mntb = peakdet(seg, 1e-15)
            # set the threshold:
            sorted_peaks = np.sort(mxtb[:, 1])

            est_hb_n = int(len(seg) / (60 / 100 * sampling_rate))
            try:
                peak_resp = sorted_peaks[-est_hb_n]
            except:
                peak_resp = sorted_peaks[-1]

            mxtb, mntb = peakdet(seg, peak_resp)

            if mxtb.shape[0] > 0:
                addpeaks = peaks[ind] + np.squeeze(mxtb[:, 0])
                try:
                    peaks = np.concatenate([peaks, addpeaks], axis=0)
                except:
                    peaks = np.concatenate((peaks, [addpeaks]))

    peaks = np.sort(peaks)
    #print(peaks.shape)
    # dpeaks = np.diff(peaks)

    # kppeaks = [i for i,v in enumerate(dpeaks) if v > delta[1]]
    ind = 0
    while ind < len(peaks) - 1:
        if peaks[ind + 1] - peaks[ind] < delta[1]:
            a = ind + np.argmin([x[int(peaks[ind])], x[int(peaks[ind + 1])]])
            peaks = np.delete(peaks, a)
        else:
            ind += 1

    # 20210312: added to remove the peaks with amplitude relatively lower than the previous and the next peaks
    peaks = np.array(list(map(lambda x: int(x), peaks)))
    ind = 1
    while ind < len(peaks) - 2:
        if (x[peaks[ind]] < x[peaks[ind + 1]] * 0.5) and (x[peaks[ind]] < x[peaks[ind - 1]] * 0.5):
            local_avg_diff = np.mean(np.diff(peaks[max(0, (ind - 10)):min(len(peaks) - 1, (ind + 10))]))
            local_std_diff = np.std(np.diff(peaks[max(0, (ind - 10)):min(len(peaks) - 1, (ind + 10))]))
            thres = local_avg_diff - 2*local_std_diff
            # check if the interval between peaks is really too short
            if (peaks[ind + 1] - peaks[ind]) < thres or (peaks[ind] - peaks[ind - 1]) < thres:
                peaks = np.delete(peaks, ind)
            else:
                ind += 1
        else:
            ind += 1

    #print(peaks.shape)
    newpeaks = peaks[[0].extend(peaks)]
    newpeaks = newpeaks.astype(int)
    newpeaks = np.squeeze(newpeaks[0, :])
    events[newpeaks] = 1
    #print(newpeaks.shape)
    # DEBUG CATCH:
    if len(newpeaks) < 5:
        print("WARNING: NUMBER OF PEAKS LOWER THAN 5!")

    return events, newpeaks

def ppg_preprocess(ppg, sampling_rate=1000, filter_type="butter", filter_band="bandpass", filter_frequency=[0.75, 3.5],
                   filter_order=1):
    # Transform to array
    ppg = np.array(ppg)

    # corr_data = outlier_correction(ppg)

    # Filter signal
    # ppg_hr_filter(physio.biopac_df['PPG'].values,[0.75, 3.5],physio.biopac_sampling_rate)
    if filter_type in ["FIR", "butter", "cheby1", "cheby2", "ellip", "bessel"]:
        # order = int(filter_order * sampling_rate)
        filtered, _, _ = biosppy.tools.filter_signal(signal=ppg,
                                                     ftype=filter_type,
                                                     band=filter_band,
                                                     order=filter_order,
                                                     frequency=filter_frequency,
                                                     sampling_rate=sampling_rate)
    else:
        filtered = ppg  # filtered is not-filtered

    # filtered = ppg_hr_filter(corr_data,[0.75, 3.5],physio.biopac_sampling_rate)
    rpeaks_signal, rpeaks = event_identify(filtered, sampling_rate)

    heart_rate_idx, heart_rate = biosppy.tools.get_heart_rate(beats=rpeaks,
                                                              sampling_rate=sampling_rate,
                                                              smooth=True,
                                                              size=7)

    # Prepare Output Dataframe
    # ==========================
    ppg_df = pd.DataFrame({"PPG_Raw": np.array(ppg)})  # Create a dataframe
    ppg_df["PPG_Filtered"] = filtered  # Add filtered signal

    # Add R peaks
    # rpeaks_signal = np.array([np.nan]*len(ppg))
    # rpeaks_signal[rpeaks_ind] = 1
    ppg_df["PPG_R_Peaks"] = rpeaks_signal

    # Get time indices
    length = len(ppg)
    T = (length - 1) / float(sampling_rate)
    ts = np.linspace(0, T, length, endpoint=False)
    heart_rate_times = ts[heart_rate_idx]
    heart_rate_times = np.round(heart_rate_times * sampling_rate).astype(int)  # Convert heart rate times to timepoints

    rpeaks_timep = ppg_df.loc[ppg_df["PPG_R_Peaks"]==1, :].index.tolist()
    rri_ms = (ts[rpeaks_timep[1:]] - ts[rpeaks_timep[:-1]]) * 1000
    rri_times = ts[rpeaks_timep[1:]]

    # Heart Rate
    try:
        heart_rate = interpolate(heart_rate, heart_rate_times, sampling_rate)  # Interpolation using 3rd order spline
        ppg_df["Heart_Rate"] = heart_rate
    except TypeError:
        print("NeuroKit Warning: ecg_process(): Sequence too short to compute heart rate.")
        ppg_df["Heart_Rate"] = np.nan

    # Store Additional Feature
    # ========================
    processed_ppg = {"df": ppg_df,
                     "PPG":{
                         "R_Peaks": rpeaks,
                         "rri_ms": rri_ms,
                         "rri_times": rri_times,
                         "times": ts
                     }
                     }
    return processed_ppg

def interpolate(values, value_times, sampling_rate=1000):
    """
    3rd order spline interpolation.
    Parameters
    ----------
    values : dataframe
        Values.
    value_times : list
        Time indices of values.
    sampling_rate : int
        Sampling rate (samples/second).
    Returns
    ----------
    signal : pd.Series
        An array containing the values indexed by time.
    Example
    ----------
    #>>> import neurokit as nk
    #>>> signal = interpolate([800, 900, 700, 500], [1000, 2000, 3000, 4000], sampling_rate=1000)
    #>>> pd.Series(signal).plot()
    Notes
    ----------
    *Authors*
    - `Dominique Makowski <https://dominiquemakowski.github.io/>`_
    *Dependencies*
    - scipy
    - pandas
    """
#    values=RRis.copy()
#    value_times=beats_times.copy()
    # Preprocessing
    initial_index = value_times[0]
    value_times = np.array(value_times) - initial_index

    # fit a 3rd degree spline on the data.
    spline = scipy.interpolate.splrep(x=value_times, y=values, k=3, s=0)  # s=0 guarantees that it will pass through ALL the given points
    x = np.arange(0, value_times[-1], 1)
    # Get the values indexed per time
    signal = scipy.interpolate.splev(x=x, tck=spline, der=0)
    # Transform to series
    signal = pd.Series(signal)
    signal.index = np.array(np.arange(initial_index, initial_index+len(signal), 1))

    return(signal)

def run(input_file, target_channel, output_dir, isplot):

    input_file = Path(input_file)
    output_dir = Path(output_dir)

    print("Rpeak processing: Loading PPG data...")
    file_name = ''.join(input_file.name.split('.')[:-1])
    if input_file.suffix == '.vhdr':
        physio = load_brainvision(input_file, [target_channel])
    elif input_file.suffix == '.acq' or input_file.suffix == '.mat':
        physio = load_biopac(input_file, target_channel)

    ylim = [50, 110]
    plot_range = 90

    sampling_rate = physio['sampling_rate']
    target_data = dict()
    target_data['signals'] = np.squeeze(physio['data'])
    target_data['times'] = np.arange(0, max(physio['data'].shape) / sampling_rate, 1 / sampling_rate)

    preproc_dir = output_dir / target_channel / "signal"
    preproc_fname_pkl = preproc_dir / (target_channel + "_preproc_" + file_name + ".pkl")
    preproc_fname_csv = preproc_dir / (target_channel + "_preproc_" + file_name + ".csv")
    preproc_fname_txt = preproc_dir / (target_channel + "_preproc_" + file_name + ".txt")

    if not preproc_fname_pkl.exists():

        preproc_dir.mkdir(parents=True, exist_ok=True)
        print("Rpeak processing: Preprocessing PPG data...")
        df = ppg_preprocess(target_data['signals'], sampling_rate)
        # df['PPG']['rri_times']
        print("Rpeak processing: Saving preprocessed PPG data...")
        save_pkl(df, preproc_fname_pkl)
    else:
        print("Rpeak processing: Loading preprocessed PPG data...")
        df = load_pkl(preproc_fname_pkl)
        # df = dict()
        # df['df'] = pd.read_csv(ppg_preproc_file_name,index_col=0)

    df['df'].to_csv(preproc_fname_csv)
    with open(preproc_fname_txt, 'w') as fp:
        fp.writelines('%d\t%f\n' % (r_sp, r_t) for r_sp, r_t in zip(df['PPG']['R_Peaks'], df['PPG']['rri_times']))

    if isplot:
        print("Rpeak processing: Plotting time series of heart rate...")

        fig_dir = output_dir / target_channel / "fig"

        (fig_dir / "hr").mkdir(parents=True, exist_ok=True)
        (fig_dir / "rpeak").mkdir(parents=True, exist_ok=True)

        fig = plt.figure(figsize=(20, 5))
        plot_signal(df['df']['Heart_Rate'], df['df'].index / sampling_rate / 60, ylim, 'time (min)', 'HR')
        fig.savefig(fig_dir / "hr" / (file_name + "_hr.png"))
        # plt.show()

        # plot r peaks, ppg signals, and HR for each minute of data
        fig_dir = fig_dir / "rpeak" / file_name
        if not fig_dir.exists():
            fig_dir.mkdir()

        print("Rpeak processing: Plotting time series of rpeaks...")
        stp = 0
        while stp < len(df['PPG']['times']):
            plt.figure(figsize=(20, 3))
            plt.subplot(2, 1, 1)
            edp = min(int(stp + plot_range * sampling_rate), len(df['PPG']['times']))
            plt.plot(df['PPG']['times'][stp:edp], df['df'].loc[stp:(edp - 1), 'PPG_Filtered'], label='filtered PPG')
            # plt.plot(df['PPG']['times'][stp:edp], df['df'].loc[stp:(edp-1),'PPG_Raw'], label='raw', c='grey')

            seg = df['df'].loc[stp:edp, ['PPG_R_Peaks']]
            peaks = seg.loc[seg.PPG_R_Peaks == 1, :].index
            plt.scatter(df['PPG']['times'][peaks], df['df'].loc[peaks, 'PPG_Filtered'], c='r')

            plt.autoscale(tight=True)
            plt.legend(loc=1)
            plt.subplot(2, 1, 2)
            plt.plot(df['PPG']['times'][stp:edp], df['df'].loc[stp:(edp - 1), 'Heart_Rate'], label='HR (bpm)')
            plt.xlim(df['PPG']['times'][stp], df['PPG']['times'][edp - 1])
            plt.ylim(ylim)

            # plt.plot(df['PPG']['times'][stp:edp], df['df'].loc[stp:(edp-1),'Heart_Rate'], label='HR (bpm)', c='r')
            # plt.autoscale(tight=True)
            plt.legend(loc=1)
            plt.savefig(fig_dir / ("rpeak_and_hr_sec_%d.png" % (stp / sampling_rate)))
            # plt.show()
            stp += int(plot_range * sampling_rate)
            plt.close()

        print("Rpeak processing: Finish")

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Process PPG file to get R peaks.')
    parser.add_argument('-i', '--input_file', help='path to PPG file.')
    parser.add_argument('-o', '--output_dir', help='directory to store output files. default: ./output', default='./output')
    parser.add_argument('-c', '--channel_name', help='name of PPG channel. default: PPG', default='PPG')
    parser.add_argument('-p', '--plot', help='option to plot time series of heart rate and detected R peaks', action="store_true")

    args = parser.parse_args()

    signal_file = args.input_file #Path("/mnt/hdd8T/Projects/COI_IAPS2/rawdata/COI_IAPS2_001/COI_IAPS2_EEG_001.vhdr")
    output_dir = args.output_dir #Path("/mnt/hdd8T/Projects/COI_IAPS2/")
    ppg_channel_name = args.channel_name #'PPG'
    isplot = args.plot #True

    run(signal_file, ppg_channel_name, output_dir, isplot)



