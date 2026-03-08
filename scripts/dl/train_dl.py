#!/usr/bin/env python3
"""
Deep Learning Training Script for Methylation Classification
Supports CNN and Transformer models with cross-validation.
"""

import argparse
import json
import pickle
import numpy as np
from pathlib import Path
from datetime import datetime

import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim.lr_scheduler import ReduceLROnPlateau, CosineAnnealingLR
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, classification_report

# Local imports
import sys
sys.path.insert(0, str(Path(__file__).parent))

from datasets import MethylationDataModule, CrossValidationDataModule
from models.cnn_1d import create_cnn_model
from models.transformer import create_transformer_model


class EarlyStopping:
    """Early stopping handler."""
    
    def __init__(self, patience: int = 10, min_delta: float = 0.001, mode: str = 'min'):
        self.patience = patience
        self.min_delta = min_delta
        self.mode = mode
        self.counter = 0
        self.best_score = None
        self.early_stop = False
    
    def __call__(self, score):
        if self.best_score is None:
            self.best_score = score
        elif self._is_improvement(score):
            self.best_score = score
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
    
    def _is_improvement(self, score):
        if self.mode == 'min':
            return score < self.best_score - self.min_delta
        else:
            return score > self.best_score + self.min_delta


def train_epoch(model, loader, criterion, optimizer, device):
    """Train for one epoch."""
    model.train()
    total_loss = 0
    all_preds = []
    all_labels = []
    
    for batch_x, batch_y in loader:
        batch_x = batch_x.to(device)
        batch_y = batch_y.to(device)
        
        optimizer.zero_grad()
        outputs = model(batch_x)
        loss = criterion(outputs, batch_y)
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
        preds = torch.argmax(outputs, dim=1)
        all_preds.extend(preds.cpu().numpy())
        all_labels.extend(batch_y.cpu().numpy())
    
    avg_loss = total_loss / len(loader)
    accuracy = accuracy_score(all_labels, all_preds)
    
    return avg_loss, accuracy


def evaluate(model, loader, criterion, device):
    """Evaluate model."""
    model.eval()
    total_loss = 0
    all_preds = []
    all_probs = []
    all_labels = []
    
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)
            
            outputs = model(batch_x)
            loss = criterion(outputs, batch_y)
            
            total_loss += loss.item()
            probs = torch.softmax(outputs, dim=1)
            preds = torch.argmax(outputs, dim=1)
            
            all_preds.extend(preds.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())
            all_labels.extend(batch_y.cpu().numpy())
    
    avg_loss = total_loss / len(loader)
    all_probs = np.array(all_probs)
    
    metrics = {
        'loss': avg_loss,
        'accuracy': accuracy_score(all_labels, all_preds),
        'f1': f1_score(all_labels, all_preds, average='weighted'),
    }
    
    # ROC-AUC
    if len(np.unique(all_labels)) == 2:
        metrics['roc_auc'] = roc_auc_score(all_labels, all_probs[:, 1])
    else:
        metrics['roc_auc'] = roc_auc_score(all_labels, all_probs, multi_class='ovr', average='weighted')
    
    return metrics, all_preds, all_probs, all_labels


def train_model(
    model,
    train_loader,
    val_loader,
    n_epochs: int,
    learning_rate: float,
    device,
    patience: int = 10,
    scheduler_type: str = 'plateau'
):
    """Full training loop."""
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=0.01)
    
    if scheduler_type == 'plateau':
        scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5)
    else:
        scheduler = CosineAnnealingLR(optimizer, T_max=n_epochs)
    
    early_stopping = EarlyStopping(patience=patience, mode='min')
    
    best_val_loss = float('inf')
    best_model_state = None
    history = {'train_loss': [], 'train_acc': [], 'val_loss': [], 'val_acc': []}
    
    for epoch in range(n_epochs):
        # Train
        train_loss, train_acc = train_epoch(model, train_loader, criterion, optimizer, device)
        
        # Validate
        val_metrics, _, _, _ = evaluate(model, val_loader, criterion, device)
        val_loss = val_metrics['loss']
        val_acc = val_metrics['accuracy']
        
        # Update scheduler
        if scheduler_type == 'plateau':
            scheduler.step(val_loss)
        else:
            scheduler.step()
        
        # Save history
        history['train_loss'].append(train_loss)
        history['train_acc'].append(train_acc)
        history['val_loss'].append(val_loss)
        history['val_acc'].append(val_acc)
        
        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_state = model.state_dict().copy()
        
        # Early stopping
        early_stopping(val_loss)
        
        if (epoch + 1) % 10 == 0 or epoch == 0:
            print(f"  Epoch {epoch+1}/{n_epochs}: "
                  f"Train Loss={train_loss:.4f}, Train Acc={train_acc:.4f}, "
                  f"Val Loss={val_loss:.4f}, Val Acc={val_acc:.4f}")
        
        if early_stopping.early_stop:
            print(f"  Early stopping at epoch {epoch+1}")
            break
    
    # Load best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
    
    return model, history


def create_model(model_type: str, n_features: int, n_classes: int, **kwargs):
    """Create model based on type."""
    if model_type == 'cnn_1d':
        return create_cnn_model(n_features, n_classes, model_type='standard', **kwargs)
    elif model_type == 'cnn_resnet':
        return create_cnn_model(n_features, n_classes, model_type='resnet', **kwargs)
    elif model_type == 'transformer':
        return create_transformer_model(n_features, n_classes, model_type='standard', **kwargs)
    elif model_type == 'tab_transformer':
        return create_transformer_model(n_features, n_classes, model_type='tab', **kwargs)
    else:
        raise ValueError(f"Unknown model type: {model_type}")


def main():
    parser = argparse.ArgumentParser(description='Train DL models for methylation classification')
    parser.add_argument('-f', '--features', required=True, help='Features file (TSV)')
    parser.add_argument('-l', '--labels', required=True, help='Labels file (TSV)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('-m', '--model', default='cnn_1d', 
                       choices=['cnn_1d', 'cnn_resnet', 'transformer', 'tab_transformer'],
                       help='Model type')
    parser.add_argument('--epochs', type=int, default=100, help='Number of epochs')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate')
    parser.add_argument('--patience', type=int, default=10, help='Early stopping patience')
    parser.add_argument('--cv_folds', type=int, default=0, help='Number of CV folds (0 for train/val/test split)')
    parser.add_argument('--no_gpu', action='store_true', help='Disable GPU')
    
    args = parser.parse_args()
    
    # Device
    device = torch.device('cuda' if torch.cuda.is_available() and not args.no_gpu else 'cpu')
    print(f"Using device: {device}")
    
    print("=== Deep Learning Training ===")
    print(f"Model: {args.model}")
    
    if args.cv_folds > 0:
        # Cross-validation mode
        print(f"\nCross-validation with {args.cv_folds} folds")
        
        data_module = CrossValidationDataModule(
            args.features, args.labels,
            n_folds=args.cv_folds,
            batch_size=args.batch_size
        )
        
        all_metrics = []
        all_predictions = np.zeros(data_module.n_samples)
        all_probabilities = np.zeros((data_module.n_samples, data_module.n_classes))
        
        for fold in range(args.cv_folds):
            print(f"\n--- Fold {fold + 1}/{args.cv_folds} ---")
            
            train_loader, val_loader = data_module.get_fold_loaders(fold)
            
            model = create_model(
                args.model,
                data_module.n_features,
                data_module.n_classes
            ).to(device)
            
            model, history = train_model(
                model, train_loader, val_loader,
                n_epochs=args.epochs,
                learning_rate=args.lr,
                device=device,
                patience=args.patience
            )
            
            # Evaluate
            val_metrics, preds, probs, labels = evaluate(
                model, val_loader, nn.CrossEntropyLoss(), device
            )
            
            all_metrics.append(val_metrics)
            
            # Store predictions for this fold
            _, val_idx = data_module.folds[fold]
            all_predictions[val_idx] = preds
            all_probabilities[val_idx] = probs
            
            print(f"  Fold {fold + 1} Results: "
                  f"Acc={val_metrics['accuracy']:.4f}, "
                  f"F1={val_metrics['f1']:.4f}, "
                  f"AUC={val_metrics['roc_auc']:.4f}")
        
        # Aggregate CV results
        final_metrics = {
            'accuracy': np.mean([m['accuracy'] for m in all_metrics]),
            'accuracy_std': np.std([m['accuracy'] for m in all_metrics]),
            'f1': np.mean([m['f1'] for m in all_metrics]),
            'f1_std': np.std([m['f1'] for m in all_metrics]),
            'roc_auc': np.mean([m['roc_auc'] for m in all_metrics]),
            'roc_auc_std': np.std([m['roc_auc'] for m in all_metrics]),
        }
        
        print(f"\n=== CV Results ===")
        print(f"Accuracy: {final_metrics['accuracy']:.4f} ± {final_metrics['accuracy_std']:.4f}")
        print(f"F1: {final_metrics['f1']:.4f} ± {final_metrics['f1_std']:.4f}")
        print(f"ROC-AUC: {final_metrics['roc_auc']:.4f} ± {final_metrics['roc_auc_std']:.4f}")
        
    else:
        # Train/val/test split mode
        print("\nUsing train/val/test split")
        
        data_module = MethylationDataModule(
            args.features, args.labels,
            batch_size=args.batch_size
        )
        
        model = create_model(
            args.model,
            data_module.n_features,
            data_module.n_classes
        ).to(device)
        
        print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")
        
        # Train
        print("\nTraining...")
        model, history = train_model(
            model,
            data_module.get_train_loader(),
            data_module.get_val_loader(),
            n_epochs=args.epochs,
            learning_rate=args.lr,
            device=device,
            patience=args.patience
        )
        
        # Test
        print("\nEvaluating on test set...")
        test_metrics, test_preds, test_probs, test_labels = evaluate(
            model, data_module.get_test_loader(), nn.CrossEntropyLoss(), device
        )
        
        final_metrics = test_metrics
        
        print(f"\n=== Test Results ===")
        print(f"Accuracy: {test_metrics['accuracy']:.4f}")
        print(f"F1: {test_metrics['f1']:.4f}")
        print(f"ROC-AUC: {test_metrics['roc_auc']:.4f}")
        print(f"\nClassification Report:")
        print(classification_report(test_labels, test_preds, 
                                   target_names=data_module.class_names))
        
        # Save model
        model_file = f"{args.output}.{args.model}.pt"
        torch.save({
            'model_state_dict': model.state_dict(),
            'model_config': {
                'model_type': args.model,
                'n_features': data_module.n_features,
                'n_classes': data_module.n_classes,
            },
            'label_encoder': data_module.label_encoder,
            'scaler': data_module.scaler,
        }, model_file)
        print(f"\nModel saved to: {model_file}")
        
        # Save history
        history_file = f"{args.output}.{args.model}.history.pkl"
        with open(history_file, 'wb') as f:
            pickle.dump(history, f)
        print(f"Training history saved to: {history_file}")
    
    # Save metrics
    metrics_file = f"{args.output}.{args.model}.metrics.json"
    with open(metrics_file, 'w') as f:
        json.dump(final_metrics, f, indent=2, default=float)
    print(f"Metrics saved to: {metrics_file}")
    
    print("\nDL training complete!")


if __name__ == '__main__':
    main()
