from dataclasses import dataclass
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import numpy as np
import arviz as az
from pymc.gp.util import plot_gp_dist
import pymc as pm
from .isoform import Isoform


@dataclass
class GaussianProcess(Isoform):

    def _gp_model(self, data, col):
        with pm.Model() as gp_model:
            p = pm.HalfCauchy("p", 3)
            n = pm.HalfCauchy("n", 3)
            M = pm.gp.mean.Constant(0)
            K = (n**2) * pm.gp.cov.ExpQuad(1, p)
            theta = pm.HalfNormal("theta", 50)
            t1_marg = pm.gp.Marginal(mean_func=M, cov_func=K)
            t1_marg.marginal_likelihood(
                "t1_likelihood",
                X=data[col].values.reshape(-1, 1),
                y=data.distances.values,
                sigma=theta,
            )

        with gp_model:
            idata = pm.sample(1000)

        az.plot_trace(idata)

        X_new = np.linspace(data[col].min(), data[col].max(), 100)

        with gp_model:
            f_pred = t1_marg.conditional("f_pred", X_new[:, None])

        with gp_model:
            gp_model_samples = pm.sample_posterior_predictive(
                idata.sel(draw=slice(0, 9)), var_names=["f_pred"]
            )

        f_pred_samples = az.extract(
            gp_model_samples, group="posterior_predictive", var_names=["f_pred"]
        )
        return X_new, f_pred_samples

    def _plot_gp_model(self, data, col, X_new, f_pred_samples):
        # plot the results
        fig = plt.figure(figsize=(12, 5))
        ax = fig.gca()

        # plot the samples from the gp posterior with samples and shading

        # f_pred_samples = az.extract(gp_model_samples, group="posterior_predictive", var_names=["f_pred"])
        plot_gp_dist(ax, samples=f_pred_samples.T, x=X_new)
        # plot the data and the true latent function
        # plt.plot(X, f_true, "dodgerblue", lw=3, label="True f")
        plt.plot(
            data[col], data.distances, "ok", ms=3, alpha=0.5, label="Observed data"
        )

        # axis labels and title
        plt.xlabel("X")
        plt.title("Posterior distribution over $f(x)$ at the observed values")
        plt.legend()

    def _sample_gp_fit(self, X_new, f_pred_samples):
        sample = []
        sample_str = 100
        for i in range(len(f_pred_samples[0].values)):
            v = []
            for j in range(len(f_pred_samples)):
                v.append(1 / f_pred_samples[j][i].values)
            df = pd.DataFrame({"x": X_new, "inv_v": v})
            df["cdf"] = df.inv_v.cumsum(axis=0) / df.inv_v.sum()
            for i in range(sample_str):
                r = np.random.rand()
                x_r = df[df.cdf > r].x.iloc[0]
                sample.append(x_r)

        plt.hist(sample)

        with pm.Model() as model:
            n = pm.Uniform("n", lower=0, upper=1)
            s = pm.Uniform("s", lower=-5, upper=5)
            x = pm.TruncatedNormal("n2", mu=n, sigma=s, lower=0, upper=1)
            l = pm.TruncatedNormal(
                "norm", mu=n, sigma=s, lower=0, upper=1, observed=sample
            )

        with model:
            sample_data = pm.sample()

        az.plot_trace(sample_data)

        hdi = az.hdi(sample_data, hdi_prob=0.95)
        return sample_data, hdi
