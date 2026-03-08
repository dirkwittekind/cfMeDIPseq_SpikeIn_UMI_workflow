#!/usr/bin/env python3
"""
PyTorch Datasets for Methylation Classification
"""

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler


class MethylationDataset(Dataset):
    """PyTorch Dataset for methylation features."""
    
    def __init__(self, features: np.ndarray, labels: np.ndarray = None, transform=None):
        """
        Args:
            features: Feature matrix (n_samples x n_features)
            labels: Label array (n_samples,)
            transform: Optional transform to apply
        """
        self.features = torch.FloatTensor(features)
        self.labels = torch.LongTensor(labels) if labels is not None else None
        self.transform = transform
    
    def __len__(self):
        return len(self.features)
    
    def __getitem__(self, idx):
        x = self.features[idx]
        
        if self.transform:
            x = self.transform(x)
        
        if self.labels is not None:
            return x, self.labels[idx]
        return x


class MethylationDataModule:
    """Data module for handling methylation data loading and splitting."""
    
    def __init__(
        self,
        features_file: str,
        labels_file: str,
        batch_size: int = 32,
        test_size: float = 0.2,
        val_size: float = 0.1,
        random_state: int = 42,
        num_workers: int = 4
    ):
        self.batch_size = batch_size
        self.test_size = test_size
        self.val_size = val_size
        self.random_state = random_state
        self.num_workers = num_workers
        
        # Load data
        self.X = pd.read_csv(features_file, sep='\t', index_col=0).values
        labels_df = pd.read_csv(labels_file, sep='\t', index_col=0)
        
        # Encode labels
        self.label_encoder = LabelEncoder()
        self.y = self.label_encoder.fit_transform(labels_df.iloc[:, 0])
        self.n_classes = len(self.label_encoder.classes_)
        self.class_names = list(self.label_encoder.classes_)
        
        # Feature info
        self.n_features = self.X.shape[1]
        self.n_samples = self.X.shape[0]
        
        # Scale features
        self.scaler = StandardScaler()
        self.X = self.scaler.fit_transform(self.X)
        
        # Split data
        self._split_data()
    
    def _split_data(self):
        """Split data into train/val/test sets."""
        # First split: train+val / test
        X_trainval, X_test, y_trainval, y_test = train_test_split(
            self.X, self.y,
            test_size=self.test_size,
            stratify=self.y,
            random_state=self.random_state
        )
        
        # Second split: train / val
        val_ratio = self.val_size / (1 - self.test_size)
        X_train, X_val, y_train, y_val = train_test_split(
            X_trainval, y_trainval,
            test_size=val_ratio,
            stratify=y_trainval,
            random_state=self.random_state
        )
        
        self.X_train, self.y_train = X_train, y_train
        self.X_val, self.y_val = X_val, y_val
        self.X_test, self.y_test = X_test, y_test
        
        print(f"Train: {len(self.X_train)}, Val: {len(self.X_val)}, Test: {len(self.X_test)}")
    
    def get_train_loader(self) -> DataLoader:
        """Get training data loader."""
        dataset = MethylationDataset(self.X_train, self.y_train)
        return DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
            pin_memory=True
        )
    
    def get_val_loader(self) -> DataLoader:
        """Get validation data loader."""
        dataset = MethylationDataset(self.X_val, self.y_val)
        return DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=True
        )
    
    def get_test_loader(self) -> DataLoader:
        """Get test data loader."""
        dataset = MethylationDataset(self.X_test, self.y_test)
        return DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=True
        )
    
    def get_full_loader(self) -> DataLoader:
        """Get loader with all data."""
        dataset = MethylationDataset(self.X, self.y)
        return DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=True
        )


class CrossValidationDataModule:
    """Data module with K-fold cross-validation support."""
    
    def __init__(
        self,
        features_file: str,
        labels_file: str,
        n_folds: int = 5,
        batch_size: int = 32,
        random_state: int = 42,
        num_workers: int = 4
    ):
        self.n_folds = n_folds
        self.batch_size = batch_size
        self.random_state = random_state
        self.num_workers = num_workers
        
        # Load data
        self.X = pd.read_csv(features_file, sep='\t', index_col=0).values
        labels_df = pd.read_csv(labels_file, sep='\t', index_col=0)
        
        # Encode labels
        self.label_encoder = LabelEncoder()
        self.y = self.label_encoder.fit_transform(labels_df.iloc[:, 0])
        self.n_classes = len(self.label_encoder.classes_)
        self.class_names = list(self.label_encoder.classes_)
        
        # Feature info
        self.n_features = self.X.shape[1]
        self.n_samples = self.X.shape[0]
        
        # Scale features
        self.scaler = StandardScaler()
        self.X = self.scaler.fit_transform(self.X)
        
        # Create CV splitter
        self.cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_state)
        self.folds = list(self.cv.split(self.X, self.y))
    
    def get_fold_loaders(self, fold_idx: int) -> tuple:
        """Get train and val loaders for a specific fold."""
        train_idx, val_idx = self.folds[fold_idx]
        
        X_train, y_train = self.X[train_idx], self.y[train_idx]
        X_val, y_val = self.X[val_idx], self.y[val_idx]
        
        train_dataset = MethylationDataset(X_train, y_train)
        val_dataset = MethylationDataset(X_val, y_val)
        
        train_loader = DataLoader(
            train_dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
            pin_memory=True
        )
        
        val_loader = DataLoader(
            val_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=True
        )
        
        return train_loader, val_loader
