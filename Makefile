#
# Simple orchestration for variant prediction experiment
#

.PHONY: clean finetune run

clean:
	


finetune:
	python run_clm.py --model_name_or_path nferruz/ProtGPT2 \
					  --train_file training.txt \
					  --validation_file validation.txt \
					  --tokenizer_name nferruz/ProtGPT2 \
					  --do_train \
					  --do_eval \
					  --output_dir calc \
					  --learning_rate 1e-06


run:
	python ./src/generate_rbd_predictions.py

