#!/usr/bin/python3
import pandas as pd
import torch
from torch.utils.data import Dataset
from torch.autograd import Variable


class scRNASeqDataset(Dataset):
    """
    Abstract class for the smsSpamCollection

    Args
        path: (string) path to the dataset
        split: (list) list of float [train_pct, valid_pct, test_pct]
    """
    def __init__(self, csv_path,  transform=None):
        self._data = pd.read_csv(csv_path, delimiter = ',')

    def __len__(self):
        return len(self._data)

    def __getitem__(self, idx):
        # Input processing
        input = self._data.values[idx]
        return input

    @property
    def input_dim_(self):
        return len(self[0])