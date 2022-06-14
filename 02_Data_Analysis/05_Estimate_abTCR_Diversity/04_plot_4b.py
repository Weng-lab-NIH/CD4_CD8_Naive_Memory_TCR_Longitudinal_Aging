import numpy as np 
import matplotlib.pyplot as plt 
import pickle


with open("new_outs/cd4_all_alpha_ratios.p", "rb") as f:
    cd4_dat = 100*np.array(pickle.load(f))

with open("new_outs/cd8_all_alpha_ratios.p", "rb") as f:
    cd8_dat = 100*np.array(pickle.load(f))


cd4_x = np.random.normal(loc=0, scale=0.125, size=len(cd4_dat))
cd8_x = np.random.normal(loc=1, scale=0.125, size=len(cd8_dat))

plt.plot(cd4_x, cd4_dat, color="blue", marker="o", linestyle="none")
plt.plot(cd8_x, cd8_dat, color="red", marker="o", linestyle="none")
plt.xticks(ticks=[0,1], labels=["CD4", "CD8"])
plt.ylim(75,100)
plt.xlim(-0.5,1.5)

plt.savefig("new_outs/both/Fig_4b.png")

## write as text
with open("new_outs/cd4_alphas.txt", "w") as f:
    for i in cd4_dat:
        f.write(str(i)+"\n")

with open("new_outs/cd8_alphas.txt", "w") as f:
    for i in cd8_dat:
        f.write(str(i)+"\n")

