#!/bin/bash -l
#SBATCH -N 1
#SBATCH -J vortex-wave
#SBATCH -p gpu-debug --gpus 2
#SBATCH -t 00:59:00
##SBATCH -p gpu --gpus 4
##SBATCH -t 0:59:00
#SBATCH -A r00043
#SBATCH --mem=256G
conda deactivate
module load python/gpu/3.12.5
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
set -euo pipefail
set -x
#nvidia-smi

CONFIG_FILE="$SCRIPT_DIR/config.yaml"
MODEL_SCRIPT="$SCRIPT_DIR/tc_gravity_waves.py"

ddp="no"
maximum_wind_values=(30 40 50 60 70)
radius_values=(30 40 50 60 70 80 90 100)
#radius_values=(40 50 60)

for maximum_wind in "${maximum_wind_values[@]}"; do
    for radius_km in "${radius_values[@]}"; do
        sed -i -E "s/^([[:space:]]*maximum_wind_m_s:)[[:space:]].*$/\1 ${maximum_wind}.0/" "$CONFIG_FILE"
        sed -i -E "s/^([[:space:]]*radius_km:)[[:space:]].*$/\1 ${radius_km}.0/" "$CONFIG_FILE"

        echo "Running maximum_wind_m_s=${maximum_wind}.0, radius_km=${radius_km}.0"
        if [ "$ddp" == "yes" ]; then
            torchrun --standalone --nproc-per-node=2 src/tc_gravity_waves.py --backend torch --device cuda --dtype float32 --ddp
        else
            python "$MODEL_SCRIPT" --config "$CONFIG_FILE"
        fi 
    done
done
