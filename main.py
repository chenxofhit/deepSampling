"""
Main experiment
"""
import json
import os
import argparse
import torch
from torch.utils.data import DataLoader
from configparser import ConfigParser
from datetime import datetime
import pandas as pd
import torch


from vae.vae import VAE
from utils.data import scRNASeqDataset

def argparser():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser(description='VAE spam detector')
    parser.add_argument('--model', type=str, required=True)
    parser.add_argument(
        '--n_epochs', type=int, default=500, 
        help="epoch numbers in training"
    )
    parser.add_argument(
        '--globals', type=str, default='./configs/globals.ini', 
        help="Path to the configuration file containing the global variables "
             "e.g. the paths to the data etc. See configs/globals.ini for an "
             "example."
    )
    parser.add_argument(
        '--config', type=str, default=None,
        help="Id of the model configuration file. If this argument is not null, "
             "the system will look for the configuration file "
             "./configs/{args.model}/{args.model}{args.config}.ini"
    )
    parser.add_argument(
        '--restore', type=str, default=None, 
        help="Path to a model checkpoint containing trained parameters. " 
             "If provided, the model will load the trained parameters before "
             "resuming training or making a prediction. By default, models are "
             "saved in ./checkpoints/<args.model><args.config>/<date>/"
    )
    return parser.parse_args()


def load_config(args):
    """
    Load .INI configuration files
    """
    config = ConfigParser()

    # Load global variable (e.g. paths)
    config.read(args.globals)

    # Path to the directory containing the model configurations
    model_config_dir = os.path.join(config['paths']['configs_directory'], '{}/'.format(args.model))

    # Load default model configuration
    default_model_config_filename = '{}.ini'.format(args.model)
    default_model_config_path = os.path.join(model_config_dir, default_model_config_filename)
    config.read(default_model_config_path)

    if args.config:
        model_config_filename = '{}{}.ini'.format(args.model, args.config)
        model_config_path = os.path.join(model_config_dir, model_config_filename)
        config.read(model_config_path)

    config.set('model', 'device', 'cuda' if torch.cuda.is_available() else 'cpu')
    if args.n_epochs is not None:
        config.set('training', 'n_epochs', str(args.n_epochs))
    return config


def train(config, trainloader, devloader=None):
    current_time = datetime.now().strftime('%Y_%m_%d_%H_%M')
    checkpoint_directory = os.path.join(
        config['paths']['checkpoints_directory'],
        '{}'.format(config['model']['name']),
        current_time)
    os.makedirs(checkpoint_directory, exist_ok=True)

    input_dim = trainloader.dataset.input_dim_ # input_dim_ refers to the input features number
    vae = VAE(input_dim, config, checkpoint_directory)
    vae.to(config['model']['device'])
    vae.fit(trainloader)

if __name__ == '__main__':
    args = argparser()
    config = load_config(args)

    # Get data path
    data_dir = config.get("paths", "data_directory")
    train_data_file_name = config.get("paths", "train_data_file_name")
    train_csv_path = os.path.join(data_dir, train_data_file_name)

    #load scRNA-seq dataset
    train_data = scRNASeqDataset(train_csv_path)

    trainloader = DataLoader(
        train_data,
        batch_size=config.getint("training", "batch_size"),
        shuffle=True,
        num_workers=0,
        pin_memory=False)

    train(config, trainloader)
