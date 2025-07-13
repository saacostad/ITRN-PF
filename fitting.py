import subprocess

import numpy as np
from scipy.optimize import minimize

expData = np.loadtxt("ExtractedData200MeV.csv", delimiter=",")

initialParams = [-102.9, 332.7, np.sqrt(0.91), np.sqrt(2.03)]

x_exp = expData[:, 0]
y_exp = np.log10(expData[:, 1])


def objective(params):
    V1, V2, a1, a2 = params

    subprocess.run(
        ["./Executable", str(V1), str(V2), str(a1), str(a2)],
        stdout=open("predictions.tsv", "w"),
        check=True,
    )

    predictions = np.loadtxt("predictions.tsv", delimiter="\t")
    y_pred = np.log10(predictions[:, 1])

    residuals = y_exp - y_pred
    return np.sum(residuals**2)


def Rsquared(y_exp, y_pred):
    ymean = np.mean(y_exp)
    RSS = np.sum((y_exp - y_pred) ** 2)
    TSS = np.sum((y_exp - ymean) ** 2)

    return 1 - (RSS / TSS)


result = minimize(objective, initialParams, method="Nelder-Mead")
fittedParams = result.x


predictions = np.loadtxt("predictions.tsv", delimiter="\t")[:, 1]
papers = np.loadtxt("data200MeV.tsv", delimiter="\t")[:, 1]


print(f"Fitted parameters: {fittedParams}")
print(f"R squared for fitting predicitons: {Rsquared(y_exp, np.log10(predictions))}")
# print(f"R squared for paper's predictions: {Rsquared(y_exp, np.log10(papers))}")


with open("200MeVData.txt", "w") as f:
    f.write(f"Best parameters = {fittedParams}\n")
    f.write(f"R2 Our Fit = {Rsquared(y_exp, np.log10(predictions))}\n")
    f.write(f"R2 papre's fit = {Rsquared(y_exp, np.log10(papers))}\n")

# print("Best parameters: ", result.x)
