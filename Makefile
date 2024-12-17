

release:
	./target/release/remastr -f data/chromosomeX.fna -n data/mini_HGVS.txt -b data/mini2.bam

test:
	RUSTFLAGS=-Awarnings cargo test

release_windows:
	cargo build --target x86_64-pc-windows-gnu --release
	# wine target/x86_64-pc-windows-gnu/release/dante_cli.exe
