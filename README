`pngwolf` is a tool to minimize the size of PNG image files. There are
a number of factors that affect the size of PNG image files, such as
the number of colors in the image and whether the image data is stored
as RGBA data or in the form of references to a color palette. The main
factor is the quality of the Deflate compression used to compress the
image data, which is in turn affected by the quality of the compressor
and how well the data to be compressed is arranged.

The PNG format supports a number of scanline filters that transform the
image data by relating nearby pixels mathematically. Choosing the right
filters for each scanline can make the image data more compressible. It
is, however, infeasible for non-trivial images to find the best filters
so typical encoders rely on a couple of heuristics to find good filters.

`pngwolf` employs a genetic algorithm to find better filter combinations
than traditional heuristics. It derives a couple of filter combinations
heuristically, adds a couple of random combinations, and then looks how
well each combination compresses. Two very different combinations may
compress similarily well, for instance, one combination may be very good
for the first couple of scanlines, while the other may be very good for
the last couple of scanlines. So taking the beginning of one combination
and the tail of the other to make a new one may result in a combination
that compresses better then the original two.

That is, in essence, what `pngwolf` does, over and over again. Further,
the most widely used PNG encoders use the zlib library for compression.
The zlib library favours speed over compression ratio in some cases, so
whatever filters are selected to aid compression, the result with zlib
may not be the smallest possible. The 7-Zip library by Igor Pavlov has
a Deflate encoder that favours size over speed at certain settings. So,
`pngwolf` attempts to make use of both: a fast zlib setting is used to
estimate how well some filter combination aids compression, and when it
gets bored, it uses 7-Zip to generate the final result.

Doing this `pngwolf` is able to compress some images better than other
optimizers (like `OptiPNG`, `AdvanceCOMP`, `pngcrush`, and `pngout`),
either because it finds better filter combinations then they do, or be-
cause it uses 7-Zip's Deflate implementation (`AdvanceCOMP` uses that
aswell, although an older version which sometimes performs better and
sometimes worse, for reasons yet to be studied). It does not attempt to
make other optimizations, like converting indexed images to RGB format.

None of the tools mentioned, including `pngwolf` follow any kind of ho-
listic approach to PNG optimization, so to get the best results they
need to be used in combination (and sometimes applying them repeatedly
or in different orders provides the best results). As far as I can tell
most other tools do not try to preserve the filter combination in the
original image, so `pngwolf` should usually be used last or second-to-
last in the optimization process.

For images that are already optimized using all the other tools, there
is about `1%` further reduction to be expected from `pngwolf` for suit-
able images. Still, it should be rare to find images on the Web that
`pngwolf` cannot compress a little bit further.

The tool suffers from the lack of a Deflate encoder that makes it easy
store the results of data analysis (where are duplicate substrings in
the data) and combine them (if you recall the earlier example where it
takes the head of one combination and the tail of another, an encoder
would not have to analyze all of the two parts again, only where they
overlap). So it can often take a long time (as in minutes) to find the
best results.

Regardless of the performance deficiency `pngwolf` is well-suited as a
research tool to come up with better heuristics for filter selection,
or to extend the genetic algorithm approach to other aspects of PNG op-
timization (the main thing being considered is re-arranging the entries
in color palettes so the image data compresses better). The tool logs
extensive information in a YAML-based machine-readable format while it
attempts to optimize images which should aid in that.

It also addresses two (other) user-interface issues I had in using the
other tools, namely it allows you to make it stop trying to find better
optimizations at well-specified points (such as the total time used),
and if you start an optimization run but grow impatient and abort the
program, results should not get lost, but should be stored anyway.

To compile `pngwolf` you need three additional libraries:

  * GAlib http://lancet.mit.edu/ga/dist/
  * 7-Zip http://www.7-zip.org/download.html ("7-Zip Source code")
  * zlib  http://zlib.net/

Put these into `galib`, `7zip`, and `zlib` sub-directories into the
directory where pngwolf.cxx is located, and then either use the CMake
utility (http://www.cmake.org/) on the `CMakeLists.txt`, or simply
specify all the files specified in `CMakeLists.txt` as input to your
compiler. The latter would look like:

  % gcc -I7zip/CPP -Igalib pngwolf.cxx galib/ga/GA1DArrayGenome.C 
      galib/ga/GAAllele.C ... -lstdc++ -o pngwolf

If you are using 7-Zip 9.20 there are two bugs in 7-Zip that prevent
gcc building https://sourceforge.net/support/tracker.php?aid=3200655
it. To address that, apply the patch `sevenzip920.patch` like so:

  % patch -p 0 < sevenzip920.patch

I've done this successfully with Visual C++ 2010 and Cygwin gcc 4.3.4.
