import matplotlib.pyplot as plt
import pandas as pd

energy = "25"

# Read the tab-separated file
df = pd.read_csv("PaperData" + energy + "MeV.tsv", sep="\t", header=None)
extracted = pd.read_csv("ExtractedData" + energy + "MeV.csv", sep=",", header=None)
dfFitted = pd.read_csv("FittedData" + energy + "MeV.tsv", sep="\t", header=None)

# Plot: y-axis in logarithmic scale
plt.figure(figsize=(10, 6))
plt.plot(
    df.iloc[:, 0],
    df.iloc[:, 1],
    linestyle="-",
    color="black",
    linewidth=3.0,
    label="Paper parameters",
)
plt.plot(
    dfFitted.iloc[:, 0],
    dfFitted.iloc[:, 1],
    linestyle="-",
    color="blue",
    linewidth=3.0,
    label="Fitted parameters",
)
plt.scatter(
    extracted.iloc[:, 0],
    extracted.iloc[:, 1],
    marker="o",
    linestyle="-",
    color="red",
    label="Exp. data",
)

plt.xscale("linear")  # or 'log' if you also want x in log scale
plt.yscale("log")

plt.xlabel(r"$\theta$ [$^\circ$]", fontsize=16)
plt.ylabel(r"$d\sigma / d\sigma_R$", fontsize=16)
plt.title(r"Cross-Section 12C-12C with $E_{lab}/A =$" + energy + "  [MeV]", fontsize=16)
plt.grid(True, which="both", ls="--")

plt.tick_params(axis="both", which="major", labelsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig("Images/" + energy + "MeV.png")  # Saves as PNG in the current directory
plt.close()  # Close the figure to free memory
