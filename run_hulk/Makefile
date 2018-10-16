dolphin: clean
	cargo run -- -q /data/kyclark/dolphin/fasta -o /data/kyclark/dolphin/hulk-out -f --create_weighted_matrix --bin_dir ../scripts

mssa_mrsa: clean
	cargo run -- -q /data/kyclark/bugs/mssa_mrsa/fasta -o /data/kyclark/bugs/hulk/mssa_mrsa -f --create_weighted_matrix --bin_dir ../scripts

sap_pyo: clean
	cargo run -- -q /data/kyclark/bugs/sap_pyo/fasta -o /data/kyclark/bugs/hulk/sap_pyo -f --create_weighted_matrix --bin_dir ../scripts

ecoli_flex: clean
	cargo run -- -q /data/kyclark/bugs/ecoli_flex/fasta -o /data/kyclark/bugs/hulk/ecoli_flex -f --create_weighted_matrix --bin_dir ../scripts

ecoli_sap: clean
	cargo run -- -q /data/kyclark/bugs/ecoli_sap/fasta -o /data/kyclark/bugs/hulk/ecoli_sap -f --create_weighted_matrix --bin_dir ../scripts

bugs: mssa_mrsa sap_pyo ecoli_flex ecoli_sap

clean:
	rm -f *.log
