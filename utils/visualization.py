#!/usr/bin/python3
"""
Utility functions for visualization
"""
import pickle
import numpy as np
import pandas as pd
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import operator
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
#plt.style.use('ggplot')
try:
    from MulticoreTSNE import MulticoreTSNE as TSNE
except BaseException:
    print("Missing MulticoreTSNE package.. Only important if evaluating other manifold learners.")


def plot_logp(logp, model):
    """
    Plot the results of the model, i.e. the precision of each method
    """
    plt.plot(logp)
    plt.legend(model)
    plt.xlabel('Iterations')
    plt.ylabel('-logp(x|z)')
    plt.savefig('./figures/logp_{}.png'.format(model))
    plt.close()


def plot_kldiv(kldiv, model):
    """
    Plot the results of the model, i.e. the precision of each method
    """
    plt.plot(kldiv)
    plt.legend(model)
    plt.xlabel('Iterations')
    plt.ylabel('KL-divergence')
    plt.savefig('./figures/kldiv_{}.png'.format(model))
    plt.close()

def hist_densities(log_densities, model):
    plt.hist(log_densities)
    plt.savefig('./figures/log_densities_{}.png'.format(model))
    plt.close()

def hist_param(params, model):
    plt.hist(params)
    plt.savefig('./figures/hist_param_{}.png'.format(model))
    plt.close()

def plot_latent(mu, model):
    fig=plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')
    ax3D.scatter(mu[:, 0], mu[:, 1], mu[:, 2])
    ax3D.set_xlabel('Latent dim 1')
    ax3D.set_ylabel('Latent dim 2') 
    ax3D.set_zlabel('Latent dim 3')            
    plt.savefig('./figures/projection_latent_{}.png'.format(model))

def plot_latent_with_annotation(mu, model, annotation):
    csv = pd.read_csv("./data/{}".format(annotation), delimiter=',')
    
    fig=plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')
    ax3D.scatter(mu[:, 0], mu[:, 1], mu[:, 2], c=csv["x"], cmap=cm.coolwarm, s=5)
    ax3D.set_xlabel('Latent dim 1')
    ax3D.set_ylabel('Latent dim 2') 
    ax3D.set_zlabel('Latent dim 3')            
    plt.legend()
    plt.savefig('./figures/projection_latent_with_annotation_{}.png'.format(model))

def plot_raw_with_rare(expr, model, log_densities, annotation="10X_PBMC.label.csv", percentage=0.1):
    x = pd.read_csv(expr, delimiter = ',')

    t0 = time.time()
    tsne = TSNE(n_jobs=16, random_state=1)
    hle = tsne.fit_transform(x)
    t1 = time.time()
    durationTsne = round(t1 - t0, ndigits=4)
    print("Total running tsne time is :" + str(durationTsne) + " s")
    
    fig=plt.figure()
    threshold = np.percentile(log_densities, percentage*100)
    idx = [i for i, x in enumerate(log_densities) if x < threshold]
    idx_opposite = [i for i, x in enumerate(log_densities) if x >= threshold]

    print("# of the rare cells is :{}, with threshold={}".format(len(idx), threshold))

    #background scatter
    plt.scatter(hle[idx, 0], hle[idx, 1], s=5, c='red', edgecolors=None, marker='.', linewidths=0, alpha=1)
    plt.scatter(hle[idx_opposite, 0], hle[idx_opposite, 1], s=5, c='lightgray', edgecolors=None, marker='.', linewidths=0, label='rare')
    
    plt.savefig('./figures/projection_raw_with_rare_{}.png'.format(model))

def plot_latent_with_rare(mu, model, log_densities, annotation="10X_PBMC.label.csv", percentage=0.1):

    csv = pd.read_csv("./data/{}".format(annotation), delimiter=',')
    
    fig=plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')

    ax3D.set_xlabel('Latent dim 1')
    ax3D.set_ylabel('Latent dim 2') 
    ax3D.set_zlabel('Latent dim 3')

    threshold = np.percentile(log_densities, percentage*100)

    idx = [i for i, x in enumerate(log_densities) if x < threshold]
    idx_opposite = [i for i, x in enumerate(log_densities) if x >= threshold]

    print("# of the rare cells is :{}, with threshold={}".format(len(idx), threshold))

    ax3D.scatter(mu[idx, 0], mu[idx, 1], mu[idx, 2],s=5, c='red', edgecolors=None, marker='.')
    ax3D.scatter(mu[idx_opposite, 0], mu[idx_opposite, 1], mu[idx_opposite, 2], s=5, c='lightgray', edgecolors=None, marker='.')

    plt.savefig('./figures/projection_latent_with_rare_{}.png'.format(model))

def load_pickle(path):
    """
    Load a dictionary containing the precision of each models
    """
    with open(path, 'rb') as pkl:
        results = pickle.load(pkl)
    return results


if __name__ == '__main__':
    pkl = "./results/10XPBMCs.pkl"
    
    res = load_pickle(pkl)
    log_densities = res["log_densities"]

    plot_raw_with_rare(expr='./data/10X_PBMC.csv', model="10XPBMCs", log_densities=log_densities)