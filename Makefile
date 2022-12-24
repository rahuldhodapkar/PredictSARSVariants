#
# Simple orchestration for variant prediction experiment
#

.PHONY: clean finetune run

clean:
	python ./src/clean_fasta.py


finetune:
	python ./src/run_clm.py --model_name_or_path nferruz/ProtGPT2 \
					  --train_file ./calc/train.txt \
					  --validation_file ./calc/test.txt \
					  --tokenizer_name nferruz/ProtGPT2 \
					  --block_size 256 \
					  --do_train \
					  --do_eval \
					  --output_dir model_out \
					  --learning_rate 1e-06


run:
	python ./src/generate_rbd_predictions.py

