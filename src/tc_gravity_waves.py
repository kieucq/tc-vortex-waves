#!/usr/bin/env python3
"""Command-line entry point for the TC vortex gravity-wave model."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from mod_spectral import load_config, run_model


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the pseudo-spectral linear shallow-water vortex model."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_name("config.yaml"),
        help="YAML configuration file (default: src/config.yaml)",
    )
    parser.add_argument(
        "--backend",
        choices=("numpy", "torch"),
        help="Execution backend; overrides compute.backend",
    )
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "cuda"),
        help="Torch device; overrides compute.device",
    )
    parser.add_argument(
        "--dtype",
        choices=("float32", "float64"),
        help="Torch precision; overrides compute.dtype",
    )
    parser.add_argument(
        "--ddp",
        action="store_true",
        help="Share FFT work across torchrun ranks (implies --backend torch)",
    )
    arguments = parser.parse_args()
    configuration_path = arguments.config.expanduser().resolve()
    configuration = load_config(configuration_path)
    compute = configuration.get("compute", {})
    use_distributed = arguments.ddp or bool(compute.get("distributed", False))
    backend = arguments.backend or str(compute.get("backend", "numpy"))
    if use_distributed:
        backend = "torch"

    if backend == "numpy":
        if int(os.environ.get("WORLD_SIZE", "1")) > 1:
            parser.error(
                "The NumPy backend cannot run under multiple torchrun ranks. "
                "Use --backend torch --ddp."
            )
        run_model(configuration, configuration_path.parent)
    else:
        from mod_spectral_ddp import run_torch_model

        run_torch_model(
            configuration,
            configuration_path.parent,
            device=arguments.device or str(compute.get("device", "auto")),
            dtype_name=arguments.dtype or str(compute.get("dtype", "float64")),
            require_distributed=use_distributed,
        )
