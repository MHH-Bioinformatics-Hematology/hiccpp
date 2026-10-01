# hiccpp

A C++20 library that reads the Juicer `.hic` format, versions 6 to 9, and writes
versions 8 and 9, without Python or Java.

Reading reproduces hicstraw 1.3.1 record for record. Writing reproduces Juicer
tools `pre` and `addNorm`: version 8 files follow release 1.22.01, version 9
files release 2.20.00, including the expected values, the VC, VC_SQRT, KR and
SCALE normalization vectors, and the restriction fragment resolutions of
`pre -f`. Blocks are compressed on several threads, and the output does not
depend on the number of threads.

```cpp title="examples/quickstart.cpp"
--8<-- "examples/quickstart.cpp"
```

```
$ quickstart GM12878_combined_30.chr21_chr22.v7.hic 21 250000
10319 records on 21 at 250000 bp, 2.75158e+07 contacts
```

[The examples page](examples.md) carries two complete programs that build with
the library: `read_hic` reads a file and `write_hic` writes one.

What the library does not do is listed in [the deviations page](DEVIATIONS.md): writing versions 6
and 7 is refused, since no obtainable Juicer tools release writes them to
validate against.
