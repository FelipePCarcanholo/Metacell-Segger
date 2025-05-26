
#!/bin/bash
# should be put in the command line as /mnt/scripts/metacell_segger_anno_AllAtOnce.sh
#them just put a space in front and put all of the arguments in the correct order
#have to choose to download from cellxgene or to use your own single cell dataset
#if use yours, please put the firs argument as "" 
#second argument is the path to store all the outputs "/your/store/path"
#third argument is the path to xenium data
#4 argument is the name of the column to pass the information, usualy "cell_type"
#5 argument is the name of the path for single cell dataset if not downloading, 
#if you put the donwload one, the 5 argument will not be used 

# Define variables
path_to_download_cellxgene="$1"  # First argument
path_to_store="$2"  # Second argument
spatial_dir="$3"  # Third argument
cell_type_column_name="$4"  # 4 argument
path_to_single_cell="$5"  # 5 argument

# Check if path_to_download_cellxgene is empty
if [ -z "$path_to_download_cellxgene" ]; then
    # If empty, keep path_to_single_cell as it is
    echo "No download needed. Using existing path: $path_to_single_cell"
else
    # If not empty, download the file and update path_to_single_cell    
    echo "Downloading from $path_to_download_cellxgene to $path_to_store/single_cell.h5ad..."
    curl -o "$path_to_store/single_cell.h5ad" "$path_to_download_cellxgene"

    # Update path_to_single_cell to the new file path
    path_to_single_cell="$path_to_store/single_cell.h5ad"

    echo "Download complete. New path_to_single_cell: $path_to_single_cell"
fi

#creating folders
mkdir -p "$path_to_store/metacell" 
mkdir -p "$path_to_store/segger" 
segger_dir="$path_to_store/segger"
mkdir -p "$segger_dir/models" 
mkdir -p "$segger_dir/segment_res" 

#export PYTHONPATH=/home/fcarcanh/.local/lib/python3.10/site-packages:$PYTHONPATH
# export PYTHONPATH=/home/fcarcanh/.local/lib/python3.10/site-packages:/mnt/scratch1/Fcarcanholo/cellxgene_python/segger_dev:$PYTHONPATH

#metacell
{ time python /mnt/scripts/Metacell_CML_cellxgene.py \
    --single_cell_path "$path_to_single_cell" \
    --metacell_output_path "$path_to_store/metacell" \
    --name_column_cell_type "$cell_type_column_name" ; } 2> "$path_to_store/metacell/metacell_time.txt"


#segger



{ time python3 /workspace/segger_dev/src/segger/cli/create_dataset_fast.py \
   --base_dir "$spatial_dir" \
   --data_dir "$segger_dir" \
   --sample_type xenium \
   --scrnaseq_file "$path_to_store/metacell/metacell_output.h5ad" \
   --celltype_column "$cell_type_column_name" \
   --k_bd 3 \
   --dist_bd 15.0 \
   --k_tx 3 \
   --dist_tx 5.0 \
   --tile_width 200 \
   --tile_height 200 \
   --neg_sampling_ratio 5.0 \
   --frac 1.0 \
   --val_prob 0.1 \
   --test_prob 0.2 \
   --n_workers 16 ; } 2> "$segger_dir/creating_time.txt" 
# subst in accelerator cuda (GPU) to cpu for this computer
{ time python3 /workspace/segger_dev/src/segger/cli/train_model.py \
    --dataset_dir "$segger_dir" \
    --models_dir "$segger_dir/models" \
    --sample_tag first_training \
    --init_emb 8 \
    --hidden_channels 32 \
    --num_tx_tokens 500 \
    --out_channels 8 \
    --heads 2 \
    --num_mid_layers 2 \
    --batch_size 4 \
    --num_workers 2 \
    --accelerator cuda \
    --max_epochs 200 \
    --devices 1 \
    --strategy auto \
    --precision 16-mixed ; } 2> "$segger_dir/training_time.txt" 
#changing knn_method from cuda to kd_tree
{ time python3 /workspace/segger_dev/src/segger/cli/predict_fast.py \
    --segger_data_dir "$segger_dir" \
    --models_dir "$segger_dir/models" \
    --benchmarks_dir "$segger_dir/segment_res" \
    --transcripts_file "$spatial_dir/transcripts.parquet" \
    --batch_size 1 \
    --num_workers 1 \
    --model_version 0 \
    --save_tag segger_embedding_1001 \
    --min_transcripts 3 \
    --cell_id_col segger_cell_id \
    --use_cc false \
    --knn_method kd_tree \
    --file_format anndata \
    --k_bd 10 \
    --dist_bd 50 \
    --k_tx 10 \
    --dist_tx 50 ; } 2> "$segger_dir/prediction_time.txt"

# passing the label information and getting the output


{ time python /mnt/scripts/transfer_label_CML.py \
    --metacell_path "$path_to_store/metacell" \
    --name_column_cell_type "$cell_type_column_name" \
    --base_dir_segger "$segger_dir" ; } 2> "$segger_dir/transfer_label.txt"
