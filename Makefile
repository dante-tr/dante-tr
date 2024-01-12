

release:
	./target/release/remastr -f data/chromosomeX.fna -n data/mini_HGVS.txt -b data/mini2.bam

test:
	RUSTFLAGS=-Awarnings cargo test
