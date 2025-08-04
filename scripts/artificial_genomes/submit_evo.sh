#!/bin/bash
#SBATCH --job-name=evo_gpu_generate
#SBATCH --output=logs/evo_generate_outputs_%x_%j.out
#SBATCH --error=logs/evo_generate_outputs_%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=80G
#SBATCH --account=izg5139_cr_default
#SBATCH --partition=standard
#SBATCH --gres=gpu:a100

source /storage/home/cpk5664/miniconda3/etc/profile.d/conda.sh
conda activate /storage/home/cpk5664/miniconda3/envs/biop

# Environment variables
export HF_HOME=/scratch/cpk5664
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export PYTHONPATH=$(pwd):$PYTHONPATH 

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <prompt_txt_file> <output_dir> <target_genome_size> <chunk_size>"
    exit 1
fi

prompt_file="$1"
output_dir="$2"
target_genome_size="$3"
chunk_size="$4"
mkdir -p "$output_dir"

echo "CUDA available: $(python -c 'import torch; print(torch.cuda.is_available())')"

# Read the initial prompt
initial_prompt=$(head -n 1 "$prompt_file")
current_genome_size=${#initial_prompt}

# Initialize final genome with the initial prompt
final_genome="$initial_prompt"

final_genome_file="${output_dir}/final_synthetic_genome.txt"

iteration=1
echo "Starting sequential generation with target size: $target_genome_size"
echo "Initial prompt size: $current_genome_size"
echo "Chunk size: $chunk_size"

while [ "$current_genome_size" -lt "$target_genome_size" ]; do
    echo "Iteration $iteration: Current genome size: $current_genome_size"
    
    # Calculate how many tokens we need for this iteration
    remaining_needed=$((target_genome_size - current_genome_size))
    
    # Use chunk_size unless we need fewer tokens to reach target
    if [ "$remaining_needed" -lt "$chunk_size" ]; then
        n_tokens="$remaining_needed"
    else
        n_tokens="$chunk_size"
    fi
    
    echo "Generating $n_tokens tokens..."
    
    # Read current prompt from file
    current_prompt=$(head -n 1 "$prompt_file")
    
    iteration_output_file="${output_dir}/iteration_${iteration}_${SLURM_JOB_ID}.txt"
    
    python scripts/generate.py \
        --model-name evo-1-8k-base \
        --prompt "$current_prompt" \
        --n-samples 1 \
        --n-tokens "$n_tokens" \
        --temperature 0.7 \
        --top-k 4 \
        --top-p 1.0 \
        --device cuda:0 \
        --cached-generation True \
        --batched True \
        --prepend-bos True \
        --verbose 0 > "$iteration_output_file"
    
    # Check if the generation was successful
    if [ $? -eq 0 ] && [ -f "$iteration_output_file" ] && [ -s "$iteration_output_file" ]; then
        echo "✓ Output saved to: $iteration_output_file"
        
        # Extract sequence after "Generated sequences:" from the saved file
        generated_sequence=$(sed -n '/Generated sequences:/,$ p' "$iteration_output_file" | tail -n +2 | tr -d '\n' | tr -d ' ')
        
        if [[ -n "$generated_sequence" ]]; then
            actual_length=${#generated_sequence}
            if [ "$actual_length" -gt "$n_tokens" ]; then
                generated_sequence="${generated_sequence:0:$n_tokens}"
                echo "Trimmed generated sequence from $actual_length to $n_tokens bases"
            fi
            
            echo "✓ Iteration $iteration completed, using ${#generated_sequence} bases"
            
            # Append to final genome
            final_genome="${final_genome}${generated_sequence}"
            current_genome_size=${#final_genome}
            
            # Update prompts.txt with the newly generated sequence as the next prompt
            echo "$generated_sequence" > "$prompt_file"
            
            # Save current state of final genome
            echo "$final_genome" > "$final_genome_file"
            
            echo "Updated prompts.txt with new sequence for next iteration"
            echo "Current final genome size: $current_genome_size"
            
            iteration=$((iteration + 1))
        else
            echo "✗ Iteration $iteration failed: Could not extract generated sequence from output file"
            echo "Content of output file:"
            cat "$iteration_output_file"
            break
        fi
    else
        echo "✗ Iteration $iteration failed: LLM generation failed or no output file created"
        if [ -f "$iteration_output_file" ]; then
            echo "Output file exists but may be empty. Content:"
            cat "$iteration_output_file"
        fi
        break
    fi
done

echo ""
echo "========== SEQUENTIAL GENERATION COMPLETED =========="
echo "Final genome size: ${#final_genome}"
echo "Target genome size: $target_genome_size"
echo "Total iterations: $((iteration - 1))"
echo "Chunk size used: $chunk_size"
echo "Final synthetic genome saved to: $final_genome_file"
echo "All iteration outputs saved in: $output_dir"