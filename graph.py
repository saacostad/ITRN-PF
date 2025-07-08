import matplotlib.pyplot as plt
import pandas as pd

# Read the tab-separated file
df = pd.read_csv("data.tsv", sep="\t")

# Check column names (optional)
print(df.head())

# Plot: y-axis in logarithmic scale
plt.figure(figsize=(8, 6))
plt.plot(df.iloc[:, 0], df.iloc[:, 1], marker="o", linestyle="-")

plt.xscale("linear")  # or 'log' if you also want x in log scale
plt.yscale("log")

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Plot from data.tsv (logarithmic scale on Y)")
plt.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()
