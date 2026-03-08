#!/usr/bin/env python3
"""
ML Classification Training Script
Trains multiple ML classifiers with nested cross-validation.
"""

import argparse
import numpy as np
import pandas as pd
import pickle
import json
from pathlib import Path

from sklearn.model_selection import StratifiedKFold, cross_val_predict, GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, confusion_matrix, classification_report,
    roc_curve, precision_recall_curve
)
import warnings
warnings.filterwarnings('ignore')


# Algorithm configurations
ALGORITHMS = {
    'random_forest': {
        'class': RandomForestClassifier,
        'params': {
            'n_estimators': [100, 200, 500],
            'max_depth': [5, 10, 20, None],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4],
            'class_weight': ['balanced', None]
        },
        'default_params': {'n_jobs': -1, 'random_state': 42}
    },
    'gradient_boosting': {
        'class': GradientBoostingClassifier,
        'params': {
            'n_estimators': [100, 200],
            'max_depth': [3, 5, 7],
            'learning_rate': [0.01, 0.1, 0.2],
            'subsample': [0.8, 1.0],
            'min_samples_split': [2, 5]
        },
        'default_params': {'random_state': 42}
    },
    'svm': {
        'class': SVC,
        'params': {
            'C': [0.1, 1, 10, 100],
            'gamma': ['scale', 'auto', 0.01, 0.1],
            'kernel': ['rbf', 'linear'],
            'class_weight': ['balanced', None]
        },
        'default_params': {'probability': True, 'random_state': 42}
    },
    'logistic_regression': {
        'class': LogisticRegression,
        'params': {
            'C': [0.01, 0.1, 1, 10, 100],
            'penalty': ['l1', 'l2'],
            'solver': ['liblinear', 'saga'],
            'class_weight': ['balanced', None]
        },
        'default_params': {'max_iter': 1000, 'random_state': 42}
    }
}


def load_data(features_file: str, labels_file: str) -> tuple:
    """Load feature matrix and labels."""
    X = pd.read_csv(features_file, sep='\t', index_col=0)
    y = pd.read_csv(labels_file, sep='\t', index_col=0).iloc[:, 0]
    
    # Align samples
    common = X.index.intersection(y.index)
    X = X.loc[common]
    y = y.loc[common]
    
    return X, y


def train_with_nested_cv(X, y, algorithm: str, outer_folds: int = 5, inner_folds: int = 3) -> dict:
    """
    Train a classifier using nested cross-validation.
    
    Outer loop: Estimate generalization performance
    Inner loop: Hyperparameter tuning
    """
    print(f"\n{'='*50}")
    print(f"Training: {algorithm}")
    print(f"Outer folds: {outer_folds}, Inner folds: {inner_folds}")
    
    algo_config = ALGORITHMS[algorithm]
    
    # Outer CV
    outer_cv = StratifiedKFold(n_splits=outer_folds, shuffle=True, random_state=42)
    inner_cv = StratifiedKFold(n_splits=inner_folds, shuffle=True, random_state=42)
    
    # Storage for results
    outer_scores = []
    all_predictions = np.zeros(len(y))
    all_probabilities = np.zeros((len(y), len(np.unique(y))))
    best_params_list = []
    feature_importances = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        print(f"\n  Outer fold {fold_idx + 1}/{outer_folds}")
        
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        
        # Inner CV for hyperparameter tuning
        base_model = algo_config['class'](**algo_config['default_params'])
        
        grid_search = GridSearchCV(
            base_model,
            algo_config['params'],
            cv=inner_cv,
            scoring='roc_auc',
            n_jobs=-1,
            verbose=0
        )
        
        grid_search.fit(X_train, y_train)
        
        best_params = grid_search.best_params_
        best_params_list.append(best_params)
        print(f"    Best params: {best_params}")
        
        # Predict on test fold
        best_model = grid_search.best_estimator_
        predictions = best_model.predict(X_test)
        probabilities = best_model.predict_proba(X_test)
        
        all_predictions[test_idx] = predictions
        all_probabilities[test_idx] = probabilities
        
        # Calculate fold score
        fold_score = roc_auc_score(y_test, probabilities[:, 1]) if len(np.unique(y)) == 2 else accuracy_score(y_test, predictions)
        outer_scores.append(fold_score)
        print(f"    Fold score: {fold_score:.4f}")
        
        # Get feature importances if available
        if hasattr(best_model, 'feature_importances_'):
            feature_importances.append(best_model.feature_importances_)
        elif hasattr(best_model, 'coef_'):
            feature_importances.append(np.abs(best_model.coef_).flatten())
    
    # Calculate final metrics
    results = calculate_metrics(y, all_predictions, all_probabilities)
    results['outer_scores'] = outer_scores
    results['mean_cv_score'] = np.mean(outer_scores)
    results['std_cv_score'] = np.std(outer_scores)
    results['best_params_per_fold'] = best_params_list
    
    # Average feature importances
    if feature_importances:
        avg_importance = np.mean(feature_importances, axis=0)
        results['feature_importances'] = dict(zip(X.columns, avg_importance))
    
    print(f"\n  Mean CV Score: {results['mean_cv_score']:.4f} ± {results['std_cv_score']:.4f}")
    
    return results


def train_final_model(X, y, algorithm: str, params: dict = None) -> tuple:
    """Train final model on all data."""
    algo_config = ALGORITHMS[algorithm]
    
    # Use provided params or default
    if params is None:
        # Quick grid search on all data
        base_model = algo_config['class'](**algo_config['default_params'])
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        
        grid_search = GridSearchCV(
            base_model,
            algo_config['params'],
            cv=cv,
            scoring='roc_auc',
            n_jobs=-1
        )
        grid_search.fit(X, y)
        params = grid_search.best_params_
    
    # Train final model
    final_params = {**algo_config['default_params'], **params}
    model = algo_config['class'](**final_params)
    model.fit(X, y)
    
    return model, params


def calculate_metrics(y_true, y_pred, y_prob) -> dict:
    """Calculate classification metrics."""
    metrics = {
        'accuracy': accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, average='weighted'),
        'recall': recall_score(y_true, y_pred, average='weighted'),
        'f1': f1_score(y_true, y_pred, average='weighted'),
    }
    
    # Binary classification specific
    if len(np.unique(y_true)) == 2:
        metrics['roc_auc'] = roc_auc_score(y_true, y_prob[:, 1])
        metrics['fpr'], metrics['tpr'], metrics['thresholds'] = roc_curve(y_true, y_prob[:, 1])
    else:
        metrics['roc_auc'] = roc_auc_score(y_true, y_prob, multi_class='ovr', average='weighted')
    
    metrics['confusion_matrix'] = confusion_matrix(y_true, y_pred).tolist()
    metrics['classification_report'] = classification_report(y_true, y_pred, output_dict=True)
    
    return metrics


def main():
    parser = argparse.ArgumentParser(description='Train ML classifiers with nested CV')
    parser.add_argument('-f', '--features', required=True, help='Features file (TSV)')
    parser.add_argument('-l', '--labels', required=True, help='Labels file (TSV)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('-a', '--algorithms', nargs='+', default=['random_forest'],
                       choices=list(ALGORITHMS.keys()), help='Algorithms to train')
    parser.add_argument('--outer_folds', type=int, default=5, help='Outer CV folds')
    parser.add_argument('--inner_folds', type=int, default=3, help='Inner CV folds')
    parser.add_argument('--train_final', action='store_true', help='Train final models on all data')
    
    args = parser.parse_args()
    
    print("=== ML Classification Training ===")
    
    # Load data
    X, y = load_data(args.features, args.labels)
    print(f"Data: {X.shape[0]} samples, {X.shape[1]} features")
    print(f"Labels: {y.value_counts().to_dict()}")
    
    # Train each algorithm
    all_results = {}
    best_model = None
    best_score = 0
    best_algorithm = None
    
    for algorithm in args.algorithms:
        results = train_with_nested_cv(
            X, y, algorithm,
            outer_folds=args.outer_folds,
            inner_folds=args.inner_folds
        )
        all_results[algorithm] = results
        
        if results['mean_cv_score'] > best_score:
            best_score = results['mean_cv_score']
            best_algorithm = algorithm
    
    # Print summary
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    
    for algo, results in all_results.items():
        print(f"\n{algo}:")
        print(f"  CV Score: {results['mean_cv_score']:.4f} ± {results['std_cv_score']:.4f}")
        print(f"  Accuracy: {results['accuracy']:.4f}")
        print(f"  F1: {results['f1']:.4f}")
        if 'roc_auc' in results:
            print(f"  ROC-AUC: {results['roc_auc']:.4f}")
    
    print(f"\nBest algorithm: {best_algorithm} (CV Score: {best_score:.4f})")
    
    # Train final models if requested
    final_models = {}
    if args.train_final:
        print("\n" + "="*50)
        print("Training final models on all data...")
        
        for algorithm in args.algorithms:
            # Use most common best params from CV
            from collections import Counter
            params_list = all_results[algorithm]['best_params_per_fold']
            params_str = [str(p) for p in params_list]
            most_common = Counter(params_str).most_common(1)[0][0]
            best_params = eval(most_common)
            
            model, params = train_final_model(X, y, algorithm, best_params)
            final_models[algorithm] = {'model': model, 'params': params}
            print(f"  {algorithm}: trained with {params}")
    
    # Save results
    # Metrics (without numpy arrays for JSON)
    metrics_file = f"{args.output}.metrics.json"
    metrics_output = {}
    for algo, results in all_results.items():
        metrics_output[algo] = {
            k: v for k, v in results.items() 
            if k not in ['fpr', 'tpr', 'thresholds', 'feature_importances']
        }
    
    with open(metrics_file, 'w') as f:
        json.dump(metrics_output, f, indent=2, default=str)
    print(f"\nMetrics saved to: {metrics_file}")
    
    # Full results (with numpy arrays)
    results_file = f"{args.output}.results.pkl"
    with open(results_file, 'wb') as f:
        pickle.dump(all_results, f)
    print(f"Full results saved to: {results_file}")
    
    # Final models
    if args.train_final:
        for algorithm, data in final_models.items():
            model_file = f"{args.output}.{algorithm}.model.pkl"
            with open(model_file, 'wb') as f:
                pickle.dump(data, f)
            print(f"Model saved to: {model_file}")
    
    # Feature importances
    for algorithm, results in all_results.items():
        if 'feature_importances' in results:
            imp_df = pd.DataFrame([
                {'feature': k, 'importance': v}
                for k, v in sorted(results['feature_importances'].items(), 
                                  key=lambda x: x[1], reverse=True)
            ])
            imp_file = f"{args.output}.{algorithm}.feature_importance.tsv"
            imp_df.to_csv(imp_file, sep='\t', index=False)
            print(f"Feature importances saved to: {imp_file}")
    
    print("\nML training complete!")


if __name__ == '__main__':
    main()
