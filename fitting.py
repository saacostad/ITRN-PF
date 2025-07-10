import subprocess

import numpy as np
from scipy.optimize import minimize

expData = np.loadtxt("ExtractedData85MeV.csv", delimiter=",")

initialParams = [-67.5, 689.6, np.sqrt(0.72), np.sqrt(2.55)]

x_exp = expData[:, 0]
y_exp = expData[:, 1]


def objective(params):
    V1, V2, a1, a2 = params

    subprocess.run(
        ["./Executable", str(V1), str(V2), str(a1), str(a2)],
        stdout=open("predictions.tsv", "w"),
        check=True,
    )

    predictions = np.loadtxt("predictions.tsv", delimiter="\t")
    y_pred = predictions[:, 1]

    residuals = y_exp - y_pred
    return np.sum(residuals**2)


result = minimize(objective, initialParams, method="Nelder-Mead")


print("Best parameters: ", result.x)
