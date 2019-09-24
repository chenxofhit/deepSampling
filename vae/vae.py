#!/usr/bin/python3
"""
Pytorch Variational Autoendoder Network Implementation
"""
from itertools import chain
import time
import json
import pickle
import numpy as np
import torch
from torch.autograd import Variable
from torch import nn
from torch import optim
from torch.nn import functional as F

class Encoder(nn.Module):
    """
    Probabilistic Encoder

    Return the mean and the variance of z ~ q(z|x). The prior
    of x is assume to be normal(0, I).

    Arguments:
        input_dim {int} -- number of features

    Returns:
        (tensor, tensor) -- mean and variance of the latent variable
            output from the forward propagation
    """
    def __init__(self, input_dim, config):
        super(Encoder, self).__init__()

        config_encoder = json.loads(config.get("encoder"))
        config_read_mu = json.loads(config.get("read_mu"))
        config_read_logvar = json.loads(config.get("read_sigma"))

        config_encoder[0]['in_features'] = input_dim

        encoder_network = []
        for layer in config_encoder:
            if layer['type'] == 'linear':
                encoder_network.append(nn.Linear(layer['in_features'], layer['out_features']))
            elif layer['type'] == 'relu':
                encoder_network.append(nn.ReLU())
            elif layer['type'] == 'tanh':
                encoder_network.append(nn.Tanh())
            elif layer['type'] == 'dropout':
                encoder_network.append(nn.Dropout(layer['rate']))
            elif layer['type'] == 'batch_norm':
                encoder_network.append(nn.BatchNorm1d(layer['num_features']))

        self.encoder_network = nn.Sequential(*encoder_network)
        self.read_mu = nn.Linear(config_read_mu['in_features'], config.getint('latent_dim'))
        self.read_logvar = nn.Linear(config_read_logvar['in_features'], config.getint('latent_dim'))
        self.initialize_parameters()

    def initialize_parameters(self):
        """
        Xavier initialization
        """
        for layer in self.modules():
            if isinstance(layer, nn.Linear):
                bound = 1 / np.sqrt(layer.in_features)
                layer.weight.data.uniform_(-bound, bound)
                layer.bias.data.zero_()

    def forward(self, inputs):
        """
        Forward propagation
        """
        hidden_state = self.encoder_network(inputs)
        mean = self.read_mu(hidden_state)
        logvar = self.read_logvar(hidden_state)
        return mean, logvar


class Decoder(nn.Module):
    """
    Decoder
    """
    def __init__(self, input_dim, config):
        super(Decoder, self).__init__()
        config_decoder = json.loads(config.get("decoder"))
        self._distr = config['distribution']

        decoder_network = []
        for layer in config_decoder:
            if layer['type'] == 'linear':
                decoder_network.append(nn.Linear(layer['in_features'], layer['out_features']))
            elif layer['type'] == 'relu':
                decoder_network.append(nn.ReLU())
            elif layer['type'] == 'relu6':
                decoder_network.append(nn.ReLU6())
            elif layer['type'] == 'tanh':
                decoder_network.append(nn.Tanh())
            elif layer['type'] == 'sigmoid':
                decoder_network.append(nn.Sigmoid())
            elif layer['type'] == 'dropout':
                decoder_network.append(nn.Dropout(layer['rate']))
            elif layer['type'] == 'batch_norm':
                decoder_network.append(nn.BatchNorm1d(layer['num_features']))
            elif layer['type'] == 'read_x':
                decoder_network.append(nn.Linear(layer['in_features'], input_dim))
        self.decoder = nn.Sequential(*decoder_network)
        if self._distr == 'poisson':
            self.read_alpha = nn.Sequential(
                nn.Linear(config.getint('latent_dim'), input_dim),
                nn.ReLU6()
            )
        self.initialize_parameters()

    def initialize_parameters(self):
        for layer in self.modules():
            if isinstance(layer, nn.Linear):
                bound = 1 / np.sqrt(layer.in_features)
                layer.weight.data.uniform_(-bound, bound)
                layer.bias.data.zero_()

    def forward(self, z):
        if self._distr == 'poisson':
            alpha = 0.5 * self.read_alpha(z)
            return alpha * self.decoder(z)
        else:
            return self.decoder(z)

class VAE(nn.Module):
    """
    VAE, x --> mu, log_sigma_sq --> N(mu, log_sigma_sq) --> z --> x
    """
    def __init__(self, input_dim, config, checkpoint_directory):
        super(VAE, self).__init__()
        self.config = config
        self.model_name = '{}'.format(config['model']['name'])
        self.checkpoint_directory = checkpoint_directory
        self._distr = config['model']['distribution']
        self._device = config['model']['device']
        self._encoder = Encoder(input_dim, config['model'])
        self._decoder = Decoder(input_dim, config['model'])

        self.num_epochs = config.getint('training', 'n_epochs')

        self._optim = optim.Adam(
            self.parameters(),
            lr=config.getfloat('training', 'lr')
        )

        self.mu = None
        self.logvar = None

        self.cur_epoch = 0
        self._save_every = config.getint('model', 'save_every')

    def parameters(self):
        return chain(self._encoder.parameters(), self._decoder.parameters())

    def _sample_z(self, mu, logvar):
        epsilon = torch.randn(mu.size())
        epsilon = Variable(epsilon, requires_grad=False).type(torch.FloatTensor).to(self._device)
        sigma = torch.exp(logvar / 2)
        return mu + sigma * epsilon

    def forward(self, inputs):
        """
        Forward propagation
        """
        self.mu, self.logvar = self._encoder(inputs)
        latent = self._sample_z(self.mu, self.logvar)
        theta = self._decoder(latent)
        return theta, self.mu, self.logvar

    def _to_numpy(self, tensor):
        return tensor.data.cpu().numpy()

    def poisson_cross_entropy(self, logtheta, inputs):
        return - inputs * logtheta + torch.exp(logtheta)

    def loglikelihood(self, reduction):
        """
        Return the log-likelihood
        """
        if self._distr == 'poisson':
            if reduction == 'none':
                return self.poisson_cross_entropy
            return nn.PoissonNLLLoss(reduction=reduction)
        elif self._distr == 'bernoulli':
            return nn.BCELoss(reduction=reduction)
        else:
            raise ValueError('{} is not a valid distribution'.format(self._distr))

    def fit(self, trainloader, print_every=1):
        """
        Train the neural network
        """

        start_time = time.time()

        storage = {
            'loss': [], 'kldiv': [], '-logp(x|z)': [],
            'log_densities': None, 'params': None
        }

        for epoch in range(self.cur_epoch, self.cur_epoch + self.num_epochs):

            self.cur_epoch += 1

            # temporary storage
            losses, kldivs, neglogliks = [], [], []

            for inputs in trainloader:
#            for inputs, _ in trainloader:
                self.train()
                inputs = inputs.to(self._device).float()
                logtheta, _, _= self.forward(inputs)
                loglikelihood = -self.loglikelihood(reduction='sum')(logtheta, inputs) / inputs.shape[0]
                kl_div = -0.5 * torch.sum(1 + self.logvar - self.mu.pow(2) - self.logvar.exp()) / inputs.shape[0]
                loss = - loglikelihood + kl_div
                loss.backward()
                self._optim.step()
                self._optim.zero_grad()
                losses.append(self._to_numpy(loss))
                kldivs.append(self._to_numpy(kl_div))
                neglogliks.append(self._to_numpy(-loglikelihood))

            storage['loss'].append(np.mean(losses))
            storage['kldiv'].append(np.mean(kldivs))
            storage['-logp(x|z)'].append(np.mean(neglogliks))

            if (epoch + 1) % print_every == 0:
                epoch_time = self._get_time(start_time, time.time())

                print('epoch: {} | loss: {:.3f} | -logp(x|z): {:.3f} | kldiv: {:.3f} | time: {}'.format(
                    epoch + 1,
                    storage['loss'][-1],
                    storage['-logp(x|z)'][-1],
                    storage['kldiv'][-1],
                    epoch_time))

        storage['log_densities'] = self._get_log_densities(trainloader)
        storage['mu'] = self._get_mu(trainloader)
        storage['params'] = self._get_parameters(trainloader)
        with open('./results/{}.pkl'.format(self.model_name), 'wb') as _f:
            pickle.dump(storage, _f, pickle.HIGHEST_PROTOCOL)

    def _get_time(self, starting_time, current_time):
        total_time = current_time - starting_time
        minutes = round(total_time // 60)
        seconds = round(total_time % 60)
        return '{} min., {} sec.'.format(minutes, seconds)

    def _get_parameters(self, dataloader):
        self.eval()
        parameters = []
        for inputs in dataloader:
            inputs = inputs.to(self._device).float()
            logtheta, _, _= self.forward(inputs)
            logtheta = self._to_numpy(logtheta)
            parameters.extend(logtheta)
        if self._distr == 'poisson':
            parameters = np.exp(np.array(parameters))
        else:
            parameters = np.array(parameters)
        return parameters

    def _get_mu(self, dataloader):
        self.eval()
        all_mu = []

        for inputs in dataloader:
            with torch.no_grad():
                inputs = inputs.to(self._device).float()
                _, mu ,_ = self.forward(inputs)
                all_mu.extend(self._to_numpy(mu))
        all_mu_array = np.array(all_mu)
        return all_mu_array


    def _get_log_densities(self, dataloader):
        all_log_densities = []
        #for inputs, _ in dataloader:
        for inputs in dataloader:
            mini_batch_log_densities = self._evaluate_probability(inputs)
            all_log_densities.extend(mini_batch_log_densities)
        all_log_densities = np.array(all_log_densities)
        return all_log_densities

    def _evaluate_probability(self, inputs):
        self.eval()
        with torch.no_grad():
            inputs = inputs.to(self._device).float()
            logtheta, _ ,_ = self.forward(inputs)
            log_likelihood = -self.loglikelihood(reduction='none')(logtheta, inputs)
            log_likelihood = torch.sum(log_likelihood, 1)
            assert inputs.shape[0] == log_likelihood.shape[0]
            return self._to_numpy(log_likelihood)


    def restore_model(self, epoch):
        """
        Retore the model parameters
        """
        model_path = '{}{}_{}.pt'.format(
            self.config['paths']['checkpoints_directory'],
            self.model_name,
            epoch)
        checkpoint = torch.load(model_path)
        self.load_state_dict(checkpoint['model_state_dict'])
        self._optim.load_state_dict(checkpoint['optimizer_state_dict'])
        self.cur_epoch = epoch