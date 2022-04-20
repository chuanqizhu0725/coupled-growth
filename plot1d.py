import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ns = 200000
step_arr = np.arange(0, ns*201, ns)
for step in step_arr:
    df = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
    arr = df[0].values
    dft = pd.read_csv(f"data/temp/1d{step}.csv", header=None)
    arrt = dft[0].values
    fig, ax1 = plt.subplots()
    ax1.plot(arr, color="blue")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot(arrt, color="red")
    plt.savefig(f"figures/phi/1d{step}")
    plt.close()
