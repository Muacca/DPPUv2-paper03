"""
Parameter Space Scanning
========================

Systematic scanning of (V, η, θ) parameter space.
"""

import csv
import os
from datetime import datetime
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Optional

import numpy as np

from .potentials import POTENTIAL_FUNCTIONS
from ..action.stability import analyze_stability, StabilityType


def calculate_single_point(args) -> Dict:
    """Calculate stability for single parameter point."""
    topology, ny_variant, V_param, eta, theta = args
    V_param = round(V_param, 2)
    eta = round(eta, 2)
    theta = round(theta, 2)

    potential_func = POTENTIAL_FUNCTIONS[(topology, ny_variant)]
    r0, delta_V, stability_type = analyze_stability(
        potential_func, V_param, eta, theta
    )

    return {
        'topology': topology,
        'ny_variant': ny_variant,
        'V': V_param,
        'eta': eta,
        'theta': theta,
        'r0': r0,
        'delta_V': delta_V,
        'stability_type': stability_type.value if stability_type else 'type-III',
    }


def generate_parameter_grid(V_range, eta_range, theta_range, topologies, ny_variants):
    """Generate all parameter combinations."""
    grid = []
    for topology in topologies:
        for ny_variant in ny_variants:
            for V_param in V_range:
                for eta in eta_range:
                    for theta in theta_range:
                        grid.append((topology, ny_variant, V_param, eta, theta))
    return grid


def run_scan(
    V_points: int = 6,
    eta_points: int = 11,
    theta_points: int = 6,
    V_min: float = 0.0,
    V_max: float = 5.0,
    eta_min: float = -5.0,
    eta_max: float = 5.0,
    theta_min: float = 0.0,
    theta_max: float = 5.0,
    topologies: Optional[List[str]] = None,
    ny_variants: Optional[List[str]] = None,
    output_dir: str = 'output',
    output_filename: Optional[str] = None,
    n_workers: Optional[int] = None
) -> List[Dict]:
    """Run full parameter scan."""
    if topologies is None:
        topologies = ['S3', 'T3', 'Nil3']
    if ny_variants is None:
        ny_variants = ['FULL', 'TT', 'REE']
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    V_range = np.linspace(V_min, V_max, V_points)
    eta_range = np.linspace(eta_min, eta_max, eta_points)
    theta_range = np.linspace(theta_min, theta_max, theta_points)

    print(f"Generating grid: {V_points}×{eta_points}×{theta_points} × "
          f"{len(topologies)} × {len(ny_variants)}")

    grid = generate_parameter_grid(
        V_range, eta_range, theta_range, topologies, ny_variants
    )
    total_points = len(grid)
    print(f"Total points: {total_points:,}")

    os.makedirs(output_dir, exist_ok=True)

    results = []
    with Pool(n_workers) as pool:
        chunk_size = 10000
        for i in range(0, total_points, chunk_size):
            chunk = grid[i:i + chunk_size]
            chunk_results = pool.map(calculate_single_point, chunk)
            results.extend(chunk_results)
            print(f"Progress: {(i + len(chunk)) / total_points * 100:.1f}%")

    # Save results
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    fname = output_filename if output_filename is not None else f'dppu_scan_all_{timestamp}.csv'
    combined_file = os.path.join(output_dir, fname)

    with open(combined_file, 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=['topology', 'ny_variant', 'V', 'eta', 'theta',
                        'r0', 'delta_V', 'stability_type']
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"Results saved to: {combined_file}")
    return results
