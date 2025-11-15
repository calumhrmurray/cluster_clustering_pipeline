"""
Cluster Clustering Pipeline

A modular pipeline for computing weighted clustering around clusters
using galaxy catalogues.
"""

from .catalogue import CatalogueManager
from .filters import CatalogueFilter
from .clustering import ClusteringAnalysis, compute_weighted_clustering
from .shapes import ShapeCorrelationAnalysis, compute_cluster_lensing
from .weights import WeightCalculator, apply_systematic_weights

__version__ = '0.1.0'

__all__ = [
    'CatalogueManager',
    'CatalogueFilter',
    'ClusteringAnalysis',
    'ShapeCorrelationAnalysis',
    'WeightCalculator',
    'compute_weighted_clustering',
    'compute_cluster_lensing',
    'apply_systematic_weights',
]
