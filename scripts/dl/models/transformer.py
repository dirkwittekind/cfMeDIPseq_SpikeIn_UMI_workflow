#!/usr/bin/env python3
"""
Transformer Model for Methylation Classification
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F


class PositionalEncoding(nn.Module):
    """Positional encoding for transformer."""
    
    def __init__(self, d_model: int, max_len: int = 10000, dropout: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        
        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        
        pe = torch.zeros(1, max_len, d_model)
        pe[0, :, 0::2] = torch.sin(position * div_term)
        pe[0, :, 1::2] = torch.cos(position * div_term)
        
        self.register_buffer('pe', pe)
    
    def forward(self, x):
        x = x + self.pe[:, :x.size(1)]
        return self.dropout(x)


class MethylationTransformer(nn.Module):
    """
    Transformer model for methylation-based classification.
    
    Treats methylation values as a sequence and uses self-attention
    to capture global patterns.
    """
    
    def __init__(
        self,
        n_features: int,
        n_classes: int,
        d_model: int = 128,
        n_heads: int = 4,
        n_layers: int = 2,
        d_ff: int = 256,
        dropout: float = 0.3,
        chunk_size: int = 100
    ):
        super().__init__()
        
        self.n_features = n_features
        self.n_classes = n_classes
        self.d_model = d_model
        self.chunk_size = chunk_size
        
        # Calculate number of chunks
        self.n_chunks = (n_features + chunk_size - 1) // chunk_size
        
        # Input projection
        self.input_projection = nn.Linear(chunk_size, d_model)
        
        # Positional encoding
        self.pos_encoder = PositionalEncoding(d_model, max_len=self.n_chunks, dropout=dropout)
        
        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_ff,
            dropout=dropout,
            batch_first=True
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)
        
        # Class token
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model))
        
        # Classifier
        self.classifier = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, n_classes)
        )
    
    def forward(self, x):
        batch_size = x.size(0)
        
        # Pad features to be divisible by chunk_size
        if self.n_features % self.chunk_size != 0:
            padding = self.chunk_size - (self.n_features % self.chunk_size)
            x = F.pad(x, (0, padding))
        
        # Reshape into chunks: (batch, n_chunks, chunk_size)
        x = x.view(batch_size, -1, self.chunk_size)
        
        # Project to d_model
        x = self.input_projection(x)
        
        # Add CLS token
        cls_tokens = self.cls_token.expand(batch_size, -1, -1)
        x = torch.cat([cls_tokens, x], dim=1)
        
        # Add positional encoding
        x = self.pos_encoder(x)
        
        # Transformer encoding
        x = self.transformer_encoder(x)
        
        # Use CLS token for classification
        cls_output = x[:, 0]
        
        # Classification
        output = self.classifier(cls_output)
        
        return output
    
    def get_attention_weights(self, x):
        """Get attention weights for interpretability."""
        batch_size = x.size(0)
        
        if self.n_features % self.chunk_size != 0:
            padding = self.chunk_size - (self.n_features % self.chunk_size)
            x = F.pad(x, (0, padding))
        
        x = x.view(batch_size, -1, self.chunk_size)
        x = self.input_projection(x)
        
        cls_tokens = self.cls_token.expand(batch_size, -1, -1)
        x = torch.cat([cls_tokens, x], dim=1)
        x = self.pos_encoder(x)
        
        # Get attention weights from first layer
        # Note: This is a simplified version; full attention extraction requires custom hooks
        attention_weights = []
        for layer in self.transformer_encoder.layers:
            attn_output, attn_weights = layer.self_attn(x, x, x, need_weights=True)
            attention_weights.append(attn_weights)
            x = layer(x)
        
        return attention_weights


class TabTransformer(nn.Module):
    """
    TabTransformer-style model for tabular methylation data.
    Uses attention across features rather than sequence positions.
    """
    
    def __init__(
        self,
        n_features: int,
        n_classes: int,
        d_model: int = 32,
        n_heads: int = 4,
        n_layers: int = 2,
        dropout: float = 0.2
    ):
        super().__init__()
        
        self.n_features = n_features
        self.n_classes = n_classes
        
        # Feature embeddings - each feature gets its own embedding
        self.feature_embeddings = nn.ModuleList([
            nn.Linear(1, d_model) for _ in range(n_features)
        ])
        
        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)
        
        # MLP for classification
        self.mlp = nn.Sequential(
            nn.Linear(n_features * d_model, 256),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, n_classes)
        )
    
    def forward(self, x):
        batch_size = x.size(0)
        
        # Embed each feature separately
        embeddings = []
        for i in range(self.n_features):
            feat = x[:, i:i+1]  # (batch, 1)
            emb = self.feature_embeddings[i](feat)  # (batch, d_model)
            embeddings.append(emb)
        
        # Stack embeddings: (batch, n_features, d_model)
        x = torch.stack(embeddings, dim=1)
        
        # Apply transformer
        x = self.transformer(x)
        
        # Flatten and classify
        x = x.view(batch_size, -1)
        output = self.mlp(x)
        
        return output


def create_transformer_model(n_features: int, n_classes: int, model_type: str = 'standard', **kwargs):
    """Factory function to create transformer models."""
    if model_type == 'standard':
        return MethylationTransformer(n_features, n_classes, **kwargs)
    elif model_type == 'tab':
        return TabTransformer(n_features, n_classes, **kwargs)
    else:
        raise ValueError(f"Unknown model type: {model_type}")
