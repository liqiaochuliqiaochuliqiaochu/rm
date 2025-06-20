#!/bin/bash
#make REBASE database
#python3 makedb.py
#makeblastdb -in dna_seqs.fasta -dbtype nucl -out rebase_dna

ls gbk > list.txt

#!/bin/bash
#Run Final_input.py
mkdir In
echo "Running Final_input.py..."
python3 Final_input.py
# Run Final_blastn.py
echo "Running Final_blastn.py..."
python3 Final_blastn.py
mkdir BLASTN
# Check if BLASTN.sh exists before submitting
if [[ -f "BLASTN.sh" ]]; then
    echo "Submitting BLASTN.sh..."
    sbatch -n 12 BLASTN.sh
    wait
else
    echo "Warning: BLASTN.sh not found. Skipping."
fi

mkdir aaIn
# Run Final_GetaaInput.py
echo "Running Final_GetaaInput.py..."
python3 Final_GetaaInput.py
wait

# Run Final_filter.py
echo "Running Final_filter.py..."
python3 Final_filter.py
wait

# Run CD.sh and makeblastdb in parallel
echo "Submitting CD.sh..."
sbatch -n 12 CD.sh &

echo "Running makeblastdb..."
makeblastdb \
    -dbtype prot \
    -in aaIn/aadb.fasta \
    -input_type fasta \
    -parse_seqids \
    -out aaIn/aa.blastdb &

wait  # Wait for both tasks to finish
mkdir BLASTP
# Run Final_blastp.py
echo "Running Final_blastp.py..."
python3 Final_blastp.py

mkdir FinalRM
# Check if BLASTP.sh exists before submitting
if [[ -f "BLASTP.sh" ]]; then
    echo "Submitting BLASTP.sh..."
    job_id=$(sbatch -n 12 BLASTP.sh | awk '{print $4}')
    echo "BLASTP.sh submitted with Job ID: $job_id"

    # Wait for BLASTP.sh to finish before proceeding
    echo "Waiting for BLASTP.sh to complete..."
    srun --dependency=afterok:$job_id python3 Final_RMcount.py
else
    echo "Warning: BLASTP.sh not found. Skipping."
    echo "Running Final_RMcount.py directly..."
    python3 Final_RMcount.py
fi

mkdir sequence 
python3 Final_ExtractSequence.py

echo "Pipeline finished!"


