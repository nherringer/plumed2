import numpy as np

data = np.loadtxt("clean.out")
print(data.shape)
data = data.reshape(-1, 852)
print(data.shape)
np.save("clean_data.npy", data)
