import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filename = '../../build/test/LIFneuronPool.csv'
data = pd.read_csv(filename, header=None, na_values='nan')
for i in range(1, len(data.columns)):
    x = data.iloc[:, i].values
    x = x[~np.isnan(x)]
    plt.eventplot(x, lineoffsets=-0.6*i, linelengths=0.5)
plt.show()

