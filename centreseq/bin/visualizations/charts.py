import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

if __name__ == "__main__":
    sns.set(color_codes=True)

    x = np.random.normal(size=100)
    sns.distplot(x)
    plt.show()
