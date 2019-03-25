#!/bin/python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('TD.csv',index_col=0)
D = df.index.to_list()
E = df.columns.to_list()
T = df.values.tolist()

fig, ax = plt.subplots()
im = ax.pcolormesh(T)
ax.set_title("Transmisión")
ax.set_xticks(np.linspace(0,float(len(E)),5))
ax.set_xticklabels(np.linspace(float(E[0]),float(E[-1]),5))
ax.set_xlabel('E')
ax.set_yticks(np.linspace(0,float(len(D)),6))
ax.set_yticklabels(np.linspace(1+float(D[0]),1+float(D[-1]),6))
ax.set_ylabel('D')
fig.colorbar(im, ax=ax)
plt.show()
