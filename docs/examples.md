# Examples

The programs below are the files under `examples/` in the repository. They are
complete: the pixel source, the write options and the data they carry are
defined in the code, and a top-level build compiles them, so an example that
stops compiling fails the build.

They are built with the library and left in the build tree:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/examples/write_hic out.hic 9
./build/examples/read_hic out.hic 10000
```

A program outside this repository builds against an installed hiccpp as
[the installation page](install.md) describes; the sources below need no change.

## Quickstart

`quickstart` opens a file and reads the records of one chromosome.

```cpp title="examples/quickstart.cpp"
--8<-- "examples/quickstart.cpp"
```

```
$ quickstart GM12878_combined_30.chr21_chr22.v7.hic 21 250000
10319 records on 21 at 250000 bp, 2.75158e+07 contacts
```

## Reading a file

`read_hic` prints the metadata of a file of version 6 to 9, reads the records of
one chromosome pair, walks the same matrix block by block on four threads, and
reads a normalization vector.

```cpp title="examples/read_hic.cpp"
--8<-- "examples/read_hic.cpp"
```

```
$ read_hic GM12878_combined_30.chr21_chr22.v6.hic 250000
version 6, genome hg19
resolutions: 2500000 1000000 500000 250000
2 chromosomes, first 21 of 48129895 bp
10319 records at 250000 bp, first 9250000 9250000 278
blocks hold 10319 records, 2.75158e+07 contacts
GW_KR vector of 193 values
```

## Writing a file

`write_hic` implements a `PixelSource` over a matrix held in memory and writes
it at two resolutions with the normalization vectors, then reads the file back.
A source that reads its pixels from elsewhere, a database or another file, hands
them over in as many batches as it likes; the writer never holds more than one
chromosome pair.

```cpp title="examples/write_hic.cpp"
--8<-- "examples/write_hic.cpp"
```

```
$ write_hic out.hic 9
wrote out.hic: version 9, resolutions 50000 10000
chr1 holds 5 records at 10000 bp
vector VC
vector VC_SQRT
vector SCALE
```
