#!/usr/bin/env python3
"""
1D Convolutional Neural Network for Methylation Classification
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class ConvBlock1D(nn.Module):
    """Convolutional block with batch norm and dropout."""
    
    def __init__(self, in_channels, out_channels, kernel_size=3, dropout=0.2):
        super().__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size, padding=kernel_size//2)
        self.bn = nn.BatchNorm1d(out_channels)
        self.dropout = nn.Dropout(dropout)
    
    def forward(self, x):
        x = self.conv(x)
        x = self.bn(x)
        x = F.relu(x)
        x = self.dropout(x)
        return x


class MethylationCNN1D(nn.Module):
    """
    1D CNN for methylation-based classification.
    
    Treats methylation features as a 1D signal and applies
    convolutional layers to capture local patterns.
    """
    
    def __init__(
        self,
        n_features: int,
        n_classes: int,
        n_filters: list = [32, 64, 128],
        kernel_sizes: list = [7, 5, 3],
        fc_units: list = [256, 128],
        dropout: float = 0.3,
        pool_size: int = 2
    ):
        super().__init__()
        
        self.n_features = n_features
        self.n_classes = n_classes
        
        # Convolutional layers
        conv_layers = []
        in_channels = 1
        
        for filters, kernel_size in zip(n_filters, kernel_sizes):
            conv_layers.append(ConvBlock1D(in_channels, filters, kernel_size, dropout))
            conv_layers.append(nn.MaxPool1d(pool_size))
            in_channels = filters
        
        self.conv_layers = nn.Sequential(*conv_layers)
        
        # Calculate size after convolutions
        with torch.no_grad():
            dummy = torch.zeros(1, 1, n_features)
            conv_out = self.conv_layers(dummy)
            self.flat_size = conv_out.view(1, -1).size(1)
        
        # Fully connected layers
        fc_layers = []
        in_units = self.flat_size
        
        for units in fc_units:
            fc_layers.append(nn.Linear(in_units, units))
            fc_layers.append(nn.ReLU())
            fc_layers.append(nn.Dropout(dropout))
            in_units = units
        
        self.fc_layers = nn.Sequential(*fc_layers)
        
        # Output layer
        self.output = nn.Linear(in_units, n_classes)
    
    def forward(self, x):
        # Add channel dimension: (batch, features) -> (batch, 1, features)
        x = x.unsqueeze(1)
        
        # Convolutional layers
        x = self.conv_layers(x)
        
        # Flatten
        x = x.view(x.size(0), -1)
        
        # Fully connected layers
        x = self.fc_layers(x)
        
        # Output
        x = self.output(x)
        
        return x
    
    def get_features(self, x):
        """Extract features before the final classification layer."""
        x = x.unsqueeze(1)
        x = self.conv_layers(x)
        x = x.view(x.size(0), -1)
        x = self.fc_layers(x)
        return x


class ResidualBlock1D(nn.Module):
    """Residual block for 1D CNN."""
    
    def __init__(self, channels, kernel_size=3, dropout=0.2):
        super().__init__()
        self.conv1 = nn.Conv1d(channels, channels, kernel_size, padding=kernel_size//2)
        self.bn1 = nn.BatchNorm1d(channels)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size, padding=kernel_size//2)
        self.bn2 = nn.BatchNorm1d(channels)
        self.dropout = nn.Dropout(dropout)
    
    def forward(self, x):
        residual = x
        x = F.relu(self.bn1(self.conv1(x)))
        x = self.dropout(x)
        x = self.bn2(self.conv2(x))
        x = x + residual
        x = F.relu(x)
        return x


class ResNetMethylation(nn.Module):
    """
    ResNet-style 1D CNN for methylation classification.
    Uses residual connections for better gradient flow.
    """
    
    def __init__(
        self,
        n_features: int,
        n_classes: int,
        n_blocks: int = 3,
        base_filters: int = 64,
        dropout: float = 0.3
    ):
        super().__init__()
        
        self.n_features = n_features
        self.n_classes = n_classes
        
        # Initial convolution
        self.init_conv = nn.Sequential(
            nn.Conv1d(1, base_filters, 7, padding=3),
            nn.BatchNorm1d(base_filters),
            nn.ReLU(),
            nn.MaxPool1d(2)
        )
        
        # Residual blocks
        res_blocks = []
        for i in range(n_blocks):
            res_blocks.append(ResidualBlock1D(base_filters, dropout=dropout))
            if i < n_blocks - 1:
                res_blocks.append(nn.MaxPool1d(2))
        
        self.res_blocks = nn.Sequential(*res_blocks)
        
        # Global average pooling
        self.global_pool = nn.AdaptiveAvgPool1d(1)
        
        # Classifier
        self.classifier = nn.Sequential(
            nn.Linear(base_filters, 128),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(128, n_classes)
        )
    
    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.init_conv(x)
        x = self.res_blocks(x)
        x = self.global_pool(x)
        x = x.view(x.size(0), -1)
        x = self.classifier(x)
        return x


def create_cnn_model(n_features: int, n_classes: int, model_type: str = 'standard', **kwargs):
    """Factory function to create CNN models."""
    if model_type == 'standard':
        return MethylationCNN1D(n_features, n_classes, **kwargs)
    elif model_type == 'resnet':
        return ResNetMethylation(n_features, n_classes, **kwargs)
    else:
        raise ValueError(f"Unknown model type: {model_type}")
