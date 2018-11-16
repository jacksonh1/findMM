
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
np.set_printoptions(suppress=True)
from optimal_colors import op_colors
plt.style.use('custom_standard')
# print(plt.style.available)
# with plt.style.context(('dark_background')):
# plt.style.use(['custom_standard', 'custom_standard_bar'])
# plt.rcParams.keys()
# plt.rcParams.update()
# %%
# ==============================================================================
# // TITLE
# ==============================================================================
test = 'out_test/example_transcriptome-30-mers.fa'

os.path.splitext(test)
