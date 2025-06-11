On linux/i3 run xsettingsd before dante
`xsettingsd`

sudo pacman -S mingw-w64-gcc
rustup target add x86_64-pc-windows-gnu
cargo build --target x86_64-pc-windows-gnu

# How to build for Windows
- compile python part:
    `cp ../../dante/dante_remastr_standalone.py ~/virtualization/shared/`
    from Windows:
    `mkdir dante_v0.11.0`
    `cp Z:/dante_remastr_standalone.py dante_v0.10.0/`
    `micromamba activate dante3`
    `cd dante_v0.10.0`
    `pyinstaller --onefile ./dante_remastr_standalone.py`
    `cp dist/dante_remastr_standalone.exe Z:`
    from Linux:
    `cp ~/virtualization/shared/dante_remastr_standalone.exe ./assets/dante_remastr_standalone.exe`

- compile to binary:
    `cargo build --target x86_64-pc-windows-gnu --release`
- test on Windows:
    `cp ../target/x86_64-pc-windows-gnu/release/dante_gui.exe ~/virtualization/shared/dante_gui.exe`
    from Windows:
    `cp Z:/dante_gui.exe ..\Desktop\`
    run

# How to build for Linux
- compile python part
    from dante directory
    `pyinstaller --onefile ./dante_remastr_standalone.py`
    `cp ./dist/dante_remastr_standalone ../remaSTR/dante_gui/assets/dante_remastr_standalone`
- compile to binary
    `cargo build --release`
- test
    `cp ../target/release/dante_gui ~/scratch/dante_testing/`
    run


