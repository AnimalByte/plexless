# Installation

Plexless uses flate2 with the pure-Rust zlib-rs backend for FASTQ gzip output.
Source builds also compile HTSlib for CRAM support and require a C toolchain,
`make`, `pkg-config`, Clang, libclang, and system zlib in addition to a recent
stable Rust toolchain and Cargo:

```bash
# Ubuntu/Debian
sudo apt install build-essential clang libclang-dev make pkg-config zlib1g-dev

# Fedora/RHEL
sudo dnf install clang clang-devel gcc make pkgconf-pkg-config zlib-devel

# Arch Linux
sudo pacman -S base-devel clang pkgconf zlib
```

## Install from crates.io

Install the version available on crates.io and verify the CLI:

```bash
cargo install plexless
plexless --help
```

`cargo install plexless` installs only the `plexless` executable into Cargo's
binary directory, normally `$HOME/.cargo/bin`.
GitHub releases and crates.io are published separately. For the v0.7.0 features,
use the GitHub binary or build the v0.7.0 source tag below.

## Download a prebuilt binary

Linux x86-64 archives are available from the public
[GitHub releases page](https://github.com/AnimalByte/plexless/releases). After
downloading the archive for the desired version:

```bash
tar -xzf plexless-vX.Y.Z-x86_64-unknown-linux-gnu.tar.gz
mkdir -p "$HOME/.local/bin"
install -m 0755 \
  plexless-vX.Y.Z-x86_64-unknown-linux-gnu/plexless \
  "$HOME/.local/bin/plexless"
plexless --help
```

The v0.7.0 Linux x86-64 binary is built on Ubuntu 22.04 and requires glibc 2.35
or newer and the zlib runtime library (`zlib1g` on Ubuntu/Debian and `zlib` on
Fedora/Arch). It is dynamically linked. Other platforms can build from source.
The release also provides the raw executable, `SHA256SUMS`, and an archive-specific
checksum file. Verify downloaded files with `sha256sum -c` before installing.
The archive includes licenses, documentation and `BUILD_INFO.txt` identifying
the source commit and build toolchain.

## Build from source

Clone the public repository and build an optimized binary:

```bash
git clone https://github.com/AnimalByte/plexless.git
cd plexless
git checkout v0.7.0
cargo build --release
```

The resulting executable is `target/release/plexless`. Alternatively, install
the current checkout through Cargo:

```bash
cargo install --path . --locked --force
```
