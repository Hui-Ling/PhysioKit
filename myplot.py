import matplotlib.pyplot as plt

def plot_signal(signal, times, ylim, xlabel, ylabel, isShow=False):
    #fig = plt.figure(figsize=(20,5))
    plot_signal_overlaid_rectangles(signal, times, [], [], ylim, xlabel, ylabel, isShow)

def plot_signal_overlaid_rectangles(signal, times, x_st, x_ed, ylim, xlabel, ylabel, isShow=False):
    #fig = plt.figure(figsize=(20,5))

    plt.plot(times,signal)
    if len(x_st)>0 or len(x_ed)>0:
        for i in range(len(x_ed)):
            rectangle = plt.Rectangle((x_st[i],ylim[0]), x_ed[i]-x_st[i], ylim[1]-ylim[0], fc='grey',ec=None,alpha=0.5)
            plt.gca().add_patch(rectangle)

    plt.xlim(times[0],times[-1])
    plt.ylim(ylim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if isShow:
        plt.show()
