# generate tutorial data example input
python generate_example_inputs.py

# SynchDP alignment with reference sequence guided
python main.py  --mode reference  --querry_set ./querry_set.pkl  --reference_seq ./reference_seq_set.pkl  --date_info ./date_info.pkl  --omics_info ./omics_info.csv  --penalty 0.1 --pcut 0.65 --window 6 --ecut 3

# SynchDP alignment with De novo reference sequence
python main.py  --mode denovo  --querry_set ./querry_set.pkl  --date_info ./date_info.pkl  --omics_info ./omics_info.csv  --cover_ratio 95  --penalty 0.01 --pcut 0.7 --window 6 --ecut 3.5

