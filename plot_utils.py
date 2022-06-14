import matplotlib.pyplot as plt
import numpy as np

def bar_plot(title, xlabel, ylabel, xdata, ydata, colors, filename, show = False):
    plt.bar(xlabel, ylabel, colors)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45)
    plt.rc('xtick', labelsize=6)
    fname = f'./plot_print/{filename}_plot.png'
    plt.savefig(fname)
    if show:
        plt.show()
    plt.close()
