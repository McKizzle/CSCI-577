import matplotlib.pyplot as plt
import numpy as np

# Plot two systems together for comparing.
# http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
def plt2dcmpr(X1, Y1, X2, Y2, styles, legend_labels, xlab, ylab, title):
    plt1, = plt.plot(X1, Y1, styles[0])
    plt2, = plt.plot(X2, Y2, styles[1])
    plt.legend([plt1, plt2], legend_labels)
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(title)




