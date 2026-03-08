"""
Deep Learning Models for Methylation Classification
"""

from .cnn_1d import MethylationCNN1D, ResNetMethylation, create_cnn_model
from .transformer import MethylationTransformer, TabTransformer, create_transformer_model

__all__ = [
    'MethylationCNN1D',
    'ResNetMethylation',
    'create_cnn_model',
    'MethylationTransformer',
    'TabTransformer',
    'create_transformer_model',
]
