import argparse

from utils import visualization

def argparser():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=str, default='10XPBMCs')
    parser.add_argument('--pkl', type=str, default='./')
    parser.add_argument('--annotation', type=str, default="10X_PBMC.label.csv")
    return parser.parse_args()

def make_plots(args):
    
    res = visualization.load_pickle(args.pkl)

    log_densities = res["log_densities"]
    logp = res["-logp(x|z)"]
    kldiv = res["kldiv"]
    mu = res["mu"]
    params = res["params"]

    model = args.model

    visualization.hist_densities(log_densities, model)
    print(log_densities)

    visualization.hist_param(params.reshape(-1), model)
    visualization.plot_logp(logp, model)
    visualization.plot_kldiv(kldiv, model)
    visualization.plot_latent(mu, model)
    visualization.plot_latent_with_annotation(mu, model, annotation="10X_PBMC.label.csv")
    visualization.plot_latent_with_rare(mu, model, log_densities, annotation="10X_PBMC.label.csv", prefix="./data", percentage=0.1)


if __name__ == '__main__':
    make_plots(argparser())
