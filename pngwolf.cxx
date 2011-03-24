////////////////////////////////////////////////////////////////////
//
// pngwolf - Optimize PNG file size by genetically finding filters
//
// Copyright (C) 2008-2011 Bjoern Hoehrmann <bjoern@hoehrmann.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// $Id$
//
////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#include <WinSock2.h>
#pragma comment(lib, "ws2_32.lib")
#else
#include <arpa/inet.h>
#endif

#include "Common/MyWindows.h"
#include "Common/MyInitGuid.h"
#include "7zip/IStream.h"
#include "7zip/Compress/ZlibEncoder.h"
#include "7zip/Common/FileStreams.h"
#include "7zip/Common/InBuffer.h"
#include "7zip/Common/StreamObjects.h"

#include <signal.h>
#include <stdlib.h>
#include <stdint.h>
#include <ga/ga.h>
#include <iostream>
#include <list>
#include <vector>
#include <stdio.h>
#include <ctime>
#include <iomanip>
#include <functional>
#include <numeric>
#include <zlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <map>
#include <bitset>

#ifdef _MSC_VER
#pragma warning(push, 4)
#pragma warning(disable: 4996)
#pragma warning(disable: 4100)
#endif

////////////////////////////////////////////////////////////////////
// Miscellaneous structures and types
////////////////////////////////////////////////////////////////////
typedef enum {
  None  = 0,
  Sub   = 1,
  Up    = 2,
  Avg   = 3,
  Paeth = 4
} PngFilter;

typedef struct PngChunk PngChunk;
struct PngChunk {
  uint32_t size;
  uint32_t type;
  uint32_t crc32;
  std::vector<char> data;
};

typedef struct IhdrChunk IhdrChunk;
struct IhdrChunk {
  uint32_t width;
  uint32_t height;
  uint32_t depth;
  uint32_t color;
  uint32_t comp;
  uint32_t filter;
  uint32_t interlace;
};

typedef GA1DArrayAlleleGenome<PngFilter> PngFilterGenome;

class PngWolf {
public:
  // IHDR data
  IhdrChunk ihdr;

  // Derived IHDR data
  size_t scanline_width;
  size_t scanline_delta;

  // The input image as list of chunks
  std::list<PngChunk> chunks;

  // Urfilter
  std::vector<PngFilter> original_filters;

  // Filters; TODO: who owns them?
  PngFilterGenome* ge_original;
  PngFilterGenome* ge_all_none;
  PngFilterGenome* ge_all_sub;
  PngFilterGenome* ge_all_up;
  PngFilterGenome* ge_all_avg;
  PngFilterGenome* ge_all_paeth;
  PngFilterGenome* ge_heuristic;
  PngFilterGenome* ge_experiment1;
  PngFilterGenome* ge_experiment2;
  PngFilterGenome* ge_experiment3;
  std::vector<PngFilterGenome*> best_genomes;

  // Command line options
  unsigned max_stagnate_time;
  unsigned max_time;
  unsigned max_evaluations;
  unsigned max_deflate;
  size_t population_size;
  const char* in_path;
  const char* out_path;
  bool verbose_analysis;
  bool verbose_summary;
  bool verbose_critters;
  bool exclude_singles;
  bool exclude_original;
  bool exclude_heuristic;
  bool include_experiment1;
  bool include_experiment2;
  bool include_experiment3;
  bool normalize_alpha;
  int zlib_level;
  int zlib_windowBits;
  int zlib_memLevel;
  int zlib_strategy;
  int szip_pass;
  int szip_fast;
  int szip_cycl;

  // User input
  bool should_abort;

  // Keeping track of time
  time_t program_begun_at;
  time_t search_begun_at;
  time_t last_improvement_at;
  time_t last_step_at;
  time_t done_deflating_at;

  // IDAT
  std::vector<char> original_inflated;
  std::vector<char> original_deflated;
  std::vector<char> original_unfiltered;

  //
  std::map<uint32_t, size_t> invis_colors;

  //
  unsigned nth_generation;
  unsigned genomes_evaluated;

  //
  std::vector<char> best_inflated;
  std::vector<char> best_deflated;

  // The genetic algorithm
  GAGeneticAlgorithm* ga;

  // Logging
  void log_analysis();
  void log_critter(PngFilterGenome* curr_best);
  void log_summary();
  void log_genome(PngFilterGenome* ge);

  // Various
  bool read_file();
  bool save_file();
  bool save_best_idat(const char* path);
  bool save_original_idat(const char* path);
  bool save_idat(const char* path, std::vector<char>& deflated, std::vector<char>& inflated);
  void init_filters();
  void run();
  void recompress();
  std::vector<char> refilter(PngFilterGenome& ge);

  // Constructor
  PngWolf() :
    should_abort(false),
    nth_generation(0),
    genomes_evaluated(0),
    ge_all_none(NULL),
    ge_all_avg(NULL),
    ge_all_sub(NULL),
    ge_all_up(NULL),
    ge_all_paeth(NULL),
    ge_heuristic(NULL),
    ge_experiment1(NULL),
    ge_experiment2(NULL),
    ge_experiment3(NULL),
    ge_original(NULL),
    done_deflating_at(0)
  {}

  ~PngWolf() {
    // TODO: This should probably delete the genomes.
  }

};

static const char PNG_MAGIC[] = "\x89\x50\x4E\x47\x0D\x0A\x1A\x0A";
static const uint32_t IDAT_TYPE = 0x49444154;
static const uint32_t IHDR_TYPE = 0x49484452;
static const uint32_t IEND_TYPE = 0x49454e44;

////////////////////////////////////////////////////////////////////
// Global PngWolf instance
////////////////////////////////////////////////////////////////////
static PngWolf wolf;

////////////////////////////////////////////////////////////////////
// PNG Scanline Filters
////////////////////////////////////////////////////////////////////

unsigned char paeth_predictor(unsigned char a, unsigned char b, unsigned char c) {
  unsigned int p = a + b - c;
  unsigned int pa = abs((int)(p - a));
  unsigned int pb = abs((int)(p - b));
  unsigned int pc = abs((int)(p - c));

  if (pa <= pb && pa <= pc)
    return a;

  if (pb <= pc)
    return b;

  return c;
}

void filter_row_none(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  while (xix < (row+1)*bytes)
    dst[xix++] = src[xix];
}

void filter_row_sub(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t aix = xix;
  xix += pwidth;
  while (xix < (row+1)*bytes)
    dst[xix++] -= src[aix++];
}

void filter_row_up(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  if (row == 0)
    return;
  while (xix < (row+1)*bytes)
    dst[xix++] -= src[bix++];
}

void filter_row_avg(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t aix;

  if (row == 0) {
    size_t aix = xix;
    xix += pwidth;
    while (xix < (row+1)*bytes)
      dst[xix++] -= src[aix++] >> 1;
    return;
  }

  aix = xix;
  for (; pwidth > 0; --pwidth)
    dst[xix++] -= src[bix++] >> 1;
  
  while (xix < (row+1)*bytes)
    dst[xix++] -= (src[aix++] + src[bix++]) >> 1;
}

void filter_row_paeth(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t aix, cix;

  if (row == 0) {
    size_t aix = xix;
    xix += pwidth;
    while (xix < (row+1)*bytes)
      dst[xix++] -= paeth_predictor(src[aix++], 0 , 0);
    return;
  }

  aix = xix;
  cix = aix - bytes;

  for (; pwidth > 0; --pwidth)
    dst[xix++] -= paeth_predictor(0, src[bix++] , 0);
  
  while (xix < (row+1)*bytes)
    dst[xix++] -= paeth_predictor(src[aix++], src[bix++], src[cix++]);
}

void unfilter_row_none(unsigned char* idat, size_t row, size_t bytes) {
  return;
}

void unfilter_row_sub(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t aix = xix;
  xix += pwidth;
  while (xix < (row+1)*bytes)
    idat[xix++] += idat[aix++];
}

void unfilter_row_up(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  if (row == 0)
    return;
  while (xix < (row+1)*bytes)
    idat[xix++] += idat[bix++];
}

void unfilter_row_avg(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t aix;

  if (row == 0) {
    size_t aix = xix;
    xix += pwidth;
    while (xix < (row+1)*bytes)
      idat[xix++] += idat[aix++] >> 1;
    return;
  }

  aix = xix;
  for (; pwidth > 0; --pwidth)
    idat[xix++] += idat[bix++] >> 1;
  
  while (xix < (row+1)*bytes)
    idat[xix++] += (idat[aix++] + idat[bix++]) >> 1;
}

void unfilter_row_paeth(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t aix, cix;

  if (row == 0) {
    size_t aix = xix;
    xix += pwidth;
    while (xix < (row+1)*bytes)
      idat[xix++] += paeth_predictor(idat[aix++], 0 , 0);
    return;
  }

  aix = xix;
  cix = aix - bytes;

  for (; pwidth > 0; --pwidth)
    idat[xix++] += paeth_predictor(0, idat[bix++] , 0);
  
  while (xix < (row+1)*bytes)
    idat[xix++] += paeth_predictor(idat[aix++], idat[bix++] , idat[cix++]);
}

void unfilter_idat(unsigned char* idat, size_t rows, size_t pwidth, size_t bytes) {
  size_t row;
  for (row = 0; row < rows; ++row) {
    switch(idat[row*bytes]) {
    case 0:
      break;
    case 1:
      unfilter_row_sub(idat, row, pwidth, bytes);
      break;
    case 2:
      unfilter_row_up(idat, row, pwidth, bytes);
      break;
    case 3:
      unfilter_row_avg(idat, row, pwidth, bytes);
      break;
    case 4:
      unfilter_row_paeth(idat, row, pwidth, bytes);
      break;
    default:
      assert(!"bad filter type");
    }
    idat[row*bytes] = 0;
  }
}

void filter_idat(unsigned char* src, unsigned char* dst, PngFilterGenome& filter, size_t pwidth, size_t bytes) {
  for (int row = 0; row < filter.size(); ++row) {
    switch(filter.gene(row)) {
    case 0:
      filter_row_none(src, dst, row, pwidth, bytes);
      break;
    case 1:
      filter_row_sub(src, dst, row, pwidth, bytes);
      break;
    case 2:
      filter_row_up(src, dst, row, pwidth, bytes);
      break;
    case 3:
      filter_row_avg(src, dst, row, pwidth, bytes);
      break;
    case 4:
      filter_row_paeth(src, dst, row, pwidth, bytes);
      break;
    default:
      assert(!"bad filter type");
    }
    dst[row*bytes] = (unsigned char)filter.gene(row);
  }
}

////////////////////////////////////////////////////////////////////
// Signal handlers
////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
BOOL WINAPI console_event_handler(DWORD Event) {
  switch (Event) {
  case CTRL_C_EVENT:
  case CTRL_BREAK_EVENT:
  case CTRL_CLOSE_EVENT:
  case CTRL_LOGOFF_EVENT:
  case CTRL_SHUTDOWN_EVENT:
    wolf.should_abort = true;
    return TRUE;
  }
  return FALSE;
}
#else
void sigint_handler(int signum) {
  wolf.should_abort = true;
}
#endif

////////////////////////////////////////////////////////////////////
// Genome Helpers
////////////////////////////////////////////////////////////////////
template <> PngFilter
GAAlleleSet<PngFilter>::allele() const {
  return (PngFilter)GARandomInt(lower(), upper());
}

std::vector<char> deflate_7zip(std::vector<char>& inflated, unsigned pass, unsigned fast, unsigned cycl) {

  NCompress::NZlib::CEncoder c;
  PROPID algoProp = NCoderPropID::kAlgorithm;
  PROPID passProp = NCoderPropID::kNumPasses;
  PROPID fastProp = NCoderPropID::kNumFastBytes;
  PROPID cyclProp = NCoderPropID::kMatchFinderCycles;

  PROPVARIANT v;
  v.vt = VT_UI4;

  // TODO: figure out what to do with errors here

  c.Create();

  NCompress::NDeflate::NEncoder::CCOMCoder* d = 
    c.DeflateEncoderSpec;

  v.ulVal = 1;
  if (d->SetCoderProperties(&algoProp, &v, 1) != S_OK) {
  }

  v.ulVal = pass;
  if (d->SetCoderProperties(&passProp, &v, 1) != S_OK) {
  }

  v.ulVal = fast;
  if (d->SetCoderProperties(&fastProp, &v, 1) != S_OK) {
  }

  v.ulVal = cycl;
  if (d->SetCoderProperties(&cyclProp, &v, 1) != S_OK) {
  }

  CBufInStream* in_buf = new CBufInStream;

  // TODO: find a way to use a a fixed buffer since we know
  // the maximum size for it and don't use more than one. It
  // might also be a good idea to keep the other objects for
  // all the passes through this to avoid re-allocations and
  // the possible failures that might go along with them.
  CDynBufSeqOutStream* out_buf = new CDynBufSeqOutStream;
  in_buf->Init((const Byte*)&inflated[0], inflated.size());
  CMyComPtr<ISequentialInStream> in(in_buf);
  CMyComPtr<ISequentialOutStream> out(out_buf);

  if (c.Code(in, out, NULL, NULL, NULL) != S_OK) {
  }

  std::vector<char> deflated(out_buf->GetSize());
  memcpy(&deflated[0], out_buf->GetBuffer(), deflated.size());

  return deflated;
}

std::vector<char> deflate_zlib(std::vector<char>& inflated) {
  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;

  // TODO: It might make sense to keep the zlib stream 
  // around and reinitialize it using deflateReset.
  if (deflateInit2(&strm, wolf.zlib_level, Z_DEFLATED,
    wolf.zlib_windowBits, wolf.zlib_memLevel,
    wolf.zlib_strategy) != Z_OK) {
    abort();
  }

  strm.next_in = (Bytef*)&inflated[0];
  strm.avail_in = inflated.size();

  size_t max = deflateBound(&strm, inflated.size());
  std::vector<char> new_deflated(max);

  strm.next_out = (Bytef*)&new_deflated[0];
  strm.avail_out = max;

  // TODO: aborting here probably leaks memory
  if (deflate(&strm, Z_FINISH) != Z_STREAM_END)
    abort();

  deflateEnd(&strm);

  new_deflated.resize(max - strm.avail_out);

  return new_deflated;
}

std::vector<char> PngWolf::refilter(PngFilterGenome& ge) {
  std::vector<char> copy(original_unfiltered.size());

  // TODO: this should be made redundant, but the filters
  // need to be modified so they copy all the data first.
  memcpy(&copy[0], &original_unfiltered[0], copy.size());

  filter_idat((unsigned char*)&original_unfiltered[0],
    (unsigned char*)&copy[0], ge,
    scanline_delta, scanline_width);

  return copy;
}

float Evaluator(GAGenome& genome) {
  PngFilterGenome& ge =
    (PngFilterGenome&)genome;

  // TODO: wolf should be user data, not a global
  if (wolf.should_abort)
    return FLT_MAX;

  wolf.genomes_evaluated++;

  std::vector<char> filtered = wolf.refilter(ge);
  std::vector<char> deflated = deflate_zlib(filtered);

  return float(deflated.size());
}

////////////////////////////////////////////////////////////////////
// Helper
////////////////////////////////////////////////////////////////////
unsigned sum_abs(unsigned c1, unsigned char c2) {
  return c1 + (c2 < 128 ? c2 : 256 - c2);
}

void PngWolf::log_genome(PngFilterGenome* ge) {
  for (int gix = 0; gix < ge->size(); ++gix) {
    if (gix % 72 == 0)
      printf("\n    ");
    printf("%1d", ge->gene(gix));
  }
}

void PngWolf::log_summary() {

  int diff = original_deflated.size() - best_deflated.size();

  if (verbose_summary) {
    printf("best filter sequence found:");
    log_genome(best_genomes.back());
    printf("\n");
    printf(""
      "best zlib deflated idat size: %0.0f\n"
      "total time spent optimizing:  %0.0f\n"
      "number of genomes evaluated:  %u\n"
      "size of 7zip deflated data:   %u\n"
      "size difference to original:  %d\n",
      best_genomes.back()->score(),
      difftime(time(NULL), program_begun_at),
      genomes_evaluated,
      best_deflated.size(),
      -diff);
  }

  if (diff >= 0)
    printf("# %u bytes smaller\n", diff);
  else
    printf("# %u bytes bigger\n", abs(diff));

  fflush(stdout);
}

void PngWolf::log_analysis() {

  printf("---\n"
    "# %u x %u pixels at depth %u (mode %u) with IDAT %u bytes (%u deflated)\n",
    ihdr.width, ihdr.height, ihdr.depth, ihdr.color,
    original_inflated.size(), original_deflated.size());

  if (!verbose_analysis)
    return;

  printf(""
    "image file path:    %s\n"
    "width in pixels:    %u\n"
    "height in pixels:   %u\n"
    "color mode:         %u\n"
    "color bit depth:    %u\n"
    "interlaced:         %u\n"
    "scanline width:     %u\n"
    "scanline delta:     %u\n"
    "inflated idat size: %u\n"
    "deflated idat size: %u\n"
    "chunks present:     ",
    this->in_path,
    this->ihdr.width,
    this->ihdr.height,
    this->ihdr.color,
    this->ihdr.depth,
    this->ihdr.interlace,
    this->scanline_width,
    this->scanline_delta,
    this->original_inflated.size(),
    this->original_deflated.size());

  std::list<PngChunk>::iterator c_it;

  for (c_it = chunks.begin(); c_it != chunks.end(); ++c_it) {
    // TODO: check that this is broken on bad endianess systems
    // Also, since no validation is performed for the types, it
    // is possible to break the YAML output with bad files, but
    // that does not seem all that important at the moment.
    printf("%c", (c_it->type >> 24));
    printf("%c", (c_it->type >> 16));
    printf("%c", (c_it->type >> 8));
    printf("%c", (c_it->type >> 0));
    printf(" ");
  }

  if (ihdr.color == 6 && ihdr.depth == 8) {
    printf("\ninvisible colors:\n");
    std::map<uint32_t, size_t>::iterator it;
    uint32_t total = 0;

    // TODO: htonl is probably not right here
    for (it = invis_colors.begin(); it != invis_colors.end(); ++it) {
      printf("  - %08X # %u times\n", htonl(it->first), it->second);
      total += it->second;
    }

    bool skip = invis_colors.size() == 1
      && invis_colors.begin()->first == 0x00000000;

    printf("  # %u pixels (%0.2f%%) are fully transparent\n",
      total, (double)total / ((double)ihdr.width * (double)ihdr.height));

    if (invis_colors.size() > 0 && !skip)
      printf("  # --normalize-alpha changes them into transparent black\n");
  } else {
    printf("\n");
  }

  printf("original filters:");
  log_genome(this->ge_original);

  printf("\n"
    "zlib deflated idat size for original filter: %0.0f\n"
    "zlib deflated idat size for None:            %0.0f\n"
    "zlib deflated idat size for Sub:             %0.0f\n"
    "zlib deflated idat size for Up:              %0.0f\n"
    "zlib deflated idat size for Avg:             %0.0f\n"
    "zlib deflated idat size for Paeth:           %0.0f\n"
    "zlib deflated idat size for experiment 1:    %0.0f\n"
    "zlib deflated idat size for experiment 2:    %0.0f\n"
    "zlib deflated idat size for experiment 3:    %0.0f\n"
    "zlib deflated idat size for basic heuristic: %0.0f\n"
    "basic heuristic filters:",
    this->ge_original->score(),
    this->ge_all_none->score(),
    this->ge_all_sub->score(),
    this->ge_all_up->score(),
    this->ge_all_avg->score(),
    this->ge_all_paeth->score(),
    this->ge_experiment1->score(),
    this->ge_experiment2->score(),
    this->ge_experiment3->score(),
    this->ge_heuristic->score());

  log_genome(this->ge_heuristic);
  printf("\nexperiment 1 filters:");
  log_genome(this->ge_experiment1);
  printf("\nexperiment 2 filters:");
  log_genome(this->ge_experiment2);
  printf("\nexperiment 3 filters:");
  log_genome(this->ge_experiment3);
  printf("\n");

  fflush(stdout);
}


void PngWolf::log_critter(PngFilterGenome* curr_best) {
  PngFilterGenome* prev_best = best_genomes.back();
  
  if (!this->verbose_critters) {
    printf(""
      "- zlib deflated idat size: %7u # %+5d bytes %+4.0f seconds since previous\n",
      unsigned(curr_best->score()),
      signed(curr_best->score() - prev_best->score()),
      difftime(time(NULL), last_improvement_at));
    return;
  }

  printf(""
    "  ##########################################################################\n"
    "- zlib deflated idat size: %7u # %+5d bytes %+4.0f seconds since previous\n"
    "  ##########################################################################\n"
    "  zlib bytes since previous improvement: %+d\n"
    "  zlib bytes since first generation:     %+d\n"
    "  seconds since program launch:          %+0.0f\n"
    "  seconds since previous improvement:    %+0.0f\n"
    "  current zlib level:                     %u\n"
    "  current zlib window bits:               %u\n"
    "  current zlib memory level:              %u\n"
    "  current zlib strategy:                  %u\n"
    "  current size of the population:         %u\n"
    "  current objective score sum:            %0.0f\n"
    "  current objective score average:        %0.0f\n"
    "  current objective score stddev:         %0.0f\n"
    "  current objective score variance:       %0.0f\n"
    "  current objective score minimum:        %0.0f\n"
    "  current objective score maximum:        %0.0f\n"
    "  current generation is the nth:          %u\n"
    "  number of genomes evaluated:            %u\n"
    "  best filters so far:",
    unsigned(curr_best->score()),
    signed(curr_best->score() - prev_best->score()),
    difftime(time(NULL), last_improvement_at),
    signed(curr_best->score() - prev_best->score()),
    signed(curr_best->score() - best_genomes.front()->score()),
    difftime(time(NULL), program_begun_at),
    difftime(time(NULL), last_improvement_at),
    zlib_level,
    zlib_windowBits,
    zlib_memLevel,
    zlib_strategy,
    ga->population().size(),
    ga->population().sum(),
    ga->population().ave(),
    ga->population().dev(),
    ga->population().var(),
    ga->population().min(),
    ga->population().max(),
    nth_generation,
    genomes_evaluated);

  log_genome(curr_best);

  printf("\n");
  fflush(stdout);
};

void PngWolf::init_filters() {

  GAAlleleSet<PngFilter> allele(None, Paeth);
  PngFilterGenome ge(ihdr.height, allele, Evaluator);

  // Copy the Urcritter to all the critters we want to hold
  // on to for anlysis and for the initial population.

  // TODO: Can clone fail? What do we do then?

  this->ge_all_avg = (PngFilterGenome*)ge.clone();
  this->ge_all_none = (PngFilterGenome*)ge.clone();
  this->ge_all_sub = (PngFilterGenome*)ge.clone();
  this->ge_all_up = (PngFilterGenome*)ge.clone();
  this->ge_all_paeth = (PngFilterGenome*)ge.clone();
  this->ge_original = (PngFilterGenome*)ge.clone();
  this->ge_heuristic = (PngFilterGenome*)ge.clone();
  this->ge_experiment1 = (PngFilterGenome*)ge.clone();
  this->ge_experiment2 = (PngFilterGenome*)ge.clone();
  this->ge_experiment3 = (PngFilterGenome*)ge.clone();

  for (int i = 0; i < ge.size(); ++i) {
    ge_original->gene(i, original_filters[i]);
    ge_all_avg->gene(i, Avg);
    ge_all_sub->gene(i, Sub);
    ge_all_none->gene(i, None);
    ge_all_paeth->gene(i, Paeth);
    ge_all_up->gene(i, Up);
  }

  std::vector<char> singles[] = {
    refilter(*ge_all_none),
    refilter(*ge_all_sub),
    refilter(*ge_all_up),
    refilter(*ge_all_avg),
    refilter(*ge_all_paeth)
  };

  for (int row = 0; row < ge.size(); ++row) {
    int best_sum = INT_MAX;
    int best_flt = 0;

    // "The following simple heuristic has performed well in
    // early tests: compute the output scanline using all five
    // filters, and select the filter that gives the smallest
    // sum of absolute values of outputs. (Consider the output
    // bytes as signed differences for this test.) This method
    // usually outperforms any single fixed filter choice." as 
    // per <http://www.w3.org/TR/PNG/#12Filter-selection>.

    // Note that I've found this to be incorrect, as far as
    // typical RGB and RGBA images found on the web go, using
    // None for all scanlines outperforms the heuristic in 57%
    // of the cases. Even if you carefully check whether they
    // should really be stored as indexed images, there is not
    // much evidence to support "usually". A better heuristic
    // would be applying the heuristic and None to all and use
    // the combination that performs better.
    for (int flt = 0; flt< 5; ++flt) {
      std::vector<char>::iterator scanline =
        singles[flt].begin() + row * scanline_width;

      int sum = std::accumulate(scanline + 1,
        scanline + scanline_width, 0, sum_abs);

      // If, for this scanline, the current filter is better
      // then the previous best filter, we memorize this filter,
      // otherwise this filter can be disregarded for the line.
      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = flt;
    }

    ge_heuristic->gene(row, (PngFilter)best_flt);
  }

  // As an experimental heuristic, this compresses each scanline
  // individually and picks the filter that compresses the line
  // best. This may be a useful clue for the others, but tests
  // suggests this might interfere in cases where zlib is a poor
  // estimator, tuning genomes too much for zlib instead of 7zip.
  // Generally this should be expected to perform poorly for very
  // small images. In the standard Alexa 1000 sample it performs
  // better than the specification's heuristic in 73% of cases;
  // files would be around 3% (median) and 4% (mean) smaller.
  for (int row = 0; row < ge.size(); ++row) {
    int best_sum = INT_MAX;
    int best_flt = 0;

    for (int flt = 0; flt < 5; ++flt) {
      std::vector<char>::iterator scanline =
        singles[flt].begin() + row * scanline_width;

      std::vector<char> line(scanline, scanline + scanline_width);
      int sum = deflate_zlib(line).size();

      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = flt;
    }

    ge_experiment1->gene(row, (PngFilter)best_flt);
  }

  // unigram heuristic
  for (int row = 0; row < ge.size(); ++row) {
    size_t best_sum = INT_MAX;
    int best_flt = 0;

    for (int flt = 0; flt< 5; ++flt) {
      std::bitset<65536> seen;
      std::vector<char>::iterator it;
      std::vector<char>::iterator scanline =
        singles[flt].begin() + row * scanline_width;

      for (it = scanline; it < scanline + scanline_width; ++it)
        seen.set(uint8_t(*it));

      size_t sum = seen.count();

      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = flt;
    }

    ge_experiment2->gene(row, (PngFilter)best_flt);
  }

  // bigram heuristic
  for (int row = 0; row < ge.size(); ++row) {
    size_t best_sum = INT_MAX;
    int best_flt = 0;

    for (int flt = 0; flt< 5; ++flt) {
      std::bitset<65536> seen;
      std::vector<char>::iterator it;
      std::vector<char>::iterator scanline =
        singles[flt].begin() + row * scanline_width;

      for (it = scanline + 1; it < scanline + scanline_width; ++it)
        seen.set((uint8_t(*(it - 1)) << 8) | uint8_t(*it));

      size_t sum = seen.count();

      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = flt;
    }

    ge_experiment3->gene(row, (PngFilter)best_flt);
  }

}

////////////////////////////////////////////////////////////////////
// Experiment
////////////////////////////////////////////////////////////////////
void PngWolf::run() {
  GAPopulation pop;

  // As initial population this uses, by default, the filters in the
  // original image, the filters derived by the heuristic proposed by
  // the PNG specification, the unary filter selections where every
  // scanline uses the same filter, and a couple of random ones re-
  // required to fill the population to the requested size. With some
  // Genetic Algorithms results vary a lot depending on the filters
  // in the initial population. For instance, improvements are hard to
  // find with some GAs if the heuristic filter selection is included.

  // TODO: for now this uses copies but there is no deallocator for
  // the originals, at some point I should figure out how to own the
  // critters. Maybe the Wolf should have a second population that
  // owns them, that way the default deallocator should handle them.

  if (!exclude_singles) {
    pop.add(*this->ge_all_none);
    pop.add(*this->ge_all_sub);
    pop.add(*this->ge_all_up);
    pop.add(*this->ge_all_avg);
    pop.add(*this->ge_all_paeth);
  }

  if (!exclude_original)
    pop.add(*this->ge_original);

  if (!exclude_heuristic)
    pop.add(*this->ge_heuristic);

  if (include_experiment1)
    pop.add(*this->ge_experiment1);

  if (include_experiment2)
    pop.add(*this->ge_experiment2);

  if (include_experiment3)
    pop.add(*this->ge_experiment3);

  // If all standard genomes have been excluded a randomized one has
  // to be added so the population knows how to make more genomes.
  if (pop.size() == 0) {
    PngFilterGenome clone(*ge_original);
    clone.initialize();
    pop.add(clone);
  }

  // This adds random critters to the initial population. Very low
  // values for the population size generally lead to insufficient
  // genetic diversity so improvements are rarely found, while with
  // very high values evolution takes too much time. I've found the
  // value here works okay-ish for reasonably sized images. There
  // is the option to make this configurable, but there is not much
  // evidence the value makes that much of a difference.
  pop.size(population_size);

  // This defines ordering by score (idat size in our case). Lower
  // idat size is better than higher idat size, setting accordingly.
  pop.order(GAPopulation::LOW_IS_BEST);

  // With what few samples I have used in testing, GAIncrementalGA
  // works very well with the other options I've used, it finds im-
  // provements fairly reliably while other Algorithms have trouble
  // to find improvements after a certain number of generations.
  GAIncrementalGA ga(pop);

  // There is no particular reason to use the tournament selector,
  // but the general concept seems sound for our purposes here, and
  // I've found this to work better than some others in my simple
  // tests, better than selecting by Rank for instance.
  ga.selector(GATournamentSelector());

  // I am not entirely sure I understand the scaling concept in the
  // GALib library, I initially did not realize fitness and score
  // weren't the same thing in GALib, so I turned scaling off to
  // make "fitness" the same as the "score". After gaining a better
  // understanding of "scaling" I tried a couple of other things,
  // but no scaling still seemed to work best for my samples.
  ga.scaling(GANoScaling());

  // This is currently not used and maybe should not be used as the
  // scoping slash ownership is unclear. It would work only during
  // run() which is a bit non-intuitive, so TODO: maybe remove this.
  this->ga = &ga;

  ga.crossover(PngFilterGenome::TwoPointCrossover);

  best_genomes.push_back((PngFilterGenome*)pop.best().clone());

  printf("---\n");

  if (ihdr.height == 1)
    goto after_while;

  while (!wolf.should_abort) {

    double since_start = difftime(time(NULL), program_begun_at);
    double since_last = difftime(time(NULL), last_improvement_at);
    size_t deflated = genomes_evaluated * original_inflated.size();

    if (max_time > 0 && max_time < since_start)
      break;

    if (max_evaluations > 0 && max_evaluations < genomes_evaluated)
      break;

    if (max_stagnate_time > 0 && max_stagnate_time < since_last)
      break;

    if (max_deflate > 0 && max_deflate < deflated / (1024*1024))
      break;

    nth_generation++;

    ga.step();
    last_step_at = time(NULL);

    PngFilterGenome& new_best =
      (PngFilterGenome&)ga.population().best();

    if (new_best.score() >= best_genomes.back()->score())
      continue;

    log_critter(&new_best);

    last_improvement_at = time(NULL);
    best_genomes.push_back((PngFilterGenome*)new_best.clone());
  }

after_while:

  // Since we intercept CTRL+C the user should get some feedback
  // on that as soon as possible, so the log header comes here.
  printf("---\n");
  fflush(stdout);
}

void PngWolf::recompress() {
  best_inflated = refilter(*best_genomes.back());
  best_deflated = deflate_7zip(best_inflated,
    szip_pass, szip_fast, szip_cycl);
  done_deflating_at = time(NULL);
}

bool PngWolf::read_file() {

  FILE* file;
  char fileMagic[8];
  std::vector<char> temp(65535);
  unsigned int iend_chunk_count = 0;
  size_t expected;

  file = fopen(this->in_path, "rb");

  if (file == NULL)
    return true;

  if (fread(fileMagic, 8, 1, file) != 1)
    goto error;

  if (memcmp(fileMagic, PNG_MAGIC, 8) != 0)
    goto error;

  while (!feof(file)) {
    PngChunk chunk;

    // TODO: This is not quite right if there are stray bytes
    // (less than the required four) at the end of the file.
    if (fread(&chunk.size, sizeof(chunk.size), 1, file) != 1)
      break;

    chunk.size = ntohl(chunk.size);

    if (fread(&chunk.type, sizeof(chunk.type), 1, file) != 1)
      goto error;

    chunk.type = ntohl(chunk.type);
    chunk.data.resize(chunk.size);

    if (chunk.size > 0)
      if (fread(&chunk.data[0], chunk.size, 1, file) != 1)
        goto error;

    if (fread(&chunk.crc32, sizeof(chunk.crc32), 1, file) != 1)
      goto error;

    chunk.crc32 = ntohl(chunk.crc32);

    // IHDR
    if (chunk.type == IHDR_TYPE && chunk.size == 13) {

      // TODO: This does not check that this is the first and only
      // IHDR chunk in the file even though only one IHDR is allowed.

      memcpy(&ihdr.width, &chunk.data.at(0), sizeof(ihdr.width));
      ihdr.width = ntohl(ihdr.width);

      memcpy(&ihdr.height, &chunk.data.at(4), sizeof(ihdr.height));
      ihdr.height = ntohl(ihdr.height);

      ihdr.depth = chunk.data.at(8);
      ihdr.color = chunk.data.at(9);
      ihdr.comp = chunk.data.at(10);
      ihdr.filter = chunk.data.at(11);
      ihdr.interlace = chunk.data.at(12);
    }

    // IDAT
    if (chunk.type == IDAT_TYPE)
      original_deflated.insert(original_deflated.end(),
        chunk.data.begin(), chunk.data.end());

    // IEND
    if (chunk.type == IEND_TYPE)
      iend_chunk_count++;

    chunks.push_back(chunk);
  }

  fclose(file);
  file = NULL;

  // We can't do anything if there is no image data in the input.
  if (original_deflated.size() == 0)
    goto error;

  // For simplicity, we rely on the image having only exactly one
  // IEND chunk (as mandated by the specification) so inserting a
  // new IDAT chunk is simple. Note that there are other possible
  // errors with the chunk arrangement that are not checked for,
  // but the worst that would happen is that a broken image is re-
  // written into a new similarily broken image, which is fine.
  if (iend_chunk_count != 1)
    goto error;

  // PNG does not allow images with zero height or width, and at
  // the time of writing, only filter mode zero was permitted.
  if (ihdr.width == 0 || ihdr.height == 0 || ihdr.filter != 0)
    goto error;

  // At the time of writing, compression level zero was the only
  // valid one, and color modes could not exceed six; interlaced
  // images are not supported, mainly because the author did not
  // bother to implement Adam7 when the goal is to minimize size.
  // Futher checks on the color mode are performed later. TODO:
  // since interlaced images are not supported, that may merit a
  // specific error message pointing that design decision out.
  if (ihdr.comp != 0 || ihdr.interlace != 0 || ihdr.color > 6)
    goto error;

  // PNG does not allow bit depths below one or above 16
  if (ihdr.depth == 0 || ihdr.depth > 16)
    goto error;

  // PNG bit depths must be a power of two
  if ((ihdr.depth - 1) & ihdr.depth)
    goto error;

  static const uint32_t channel_map[] = {
    1, 0, 3, 1, 2, 0, 4
  };

  // This validates generally permissable color modes. It does not
  // fully check whether the combination of bith depth and color
  // mode is permitted by the specification. TODO: Maybe it should.
  if (channel_map[ihdr.color] < 1)
    goto error;

  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.next_in = (Bytef*)&original_deflated[0];
  strm.avail_in = original_deflated.size();

  if (inflateInit(&strm) != Z_OK)
    goto error;

  do {
      strm.avail_out = temp.size();
      strm.next_out = (Bytef*)&temp[0];
      int ret = inflate(&strm, Z_NO_FLUSH);

      // TODO: going to `error` here probably leaks some memory
      // but it would be freed when exiting the process, so this
      // is mostly important when turning this into a library.
      if (ret != Z_STREAM_END && ret != Z_OK)
        goto error;

      size_t have = temp.size() - strm.avail_out;
      original_inflated.insert(original_inflated.end(),
        temp.begin(), temp.begin() + have);

  } while (strm.avail_out == 0);

  if (inflateEnd(&strm) != Z_OK)
    goto error;

  expected = ihdr.height * int(ceil(
    float(((ihdr.width * channel_map[ihdr.color] * ihdr.depth + 8) / 8.0f))
  ));

  if (expected != original_inflated.size())
    goto error;

  scanline_width = expected / ihdr.height;
  scanline_delta = channel_map[ihdr.color] *
    int(ceil(float(ihdr.depth / 8.0)));

  for (size_t ix = 0; ix < ihdr.height; ++ix) {
    unsigned char filter = original_inflated.at(ix * scanline_width);

    // Abort when the image uses an unsupported filter type
    if (filter > 4)
      goto error;

    original_filters.push_back((PngFilter)filter);
  }

  // TODO: copy properly here
  original_unfiltered.resize(original_inflated.size());
  memcpy(&original_unfiltered[0],
    &original_inflated[0], original_inflated.size());

  unfilter_idat((unsigned char*)&original_unfiltered[0],
    ihdr.height, scanline_delta, scanline_width);

  // Creates a histogram of the fully transparent pixels in the
  // image. It is apparently common for graphics programs to
  // keep the color values of fully transparent pixels around,
  // but this is rarely desired and makes compression harder, so
  // we tell users about that and offer to normalize the pixels.

  // TODO: Now that it is modified, original_unfiltered is the
  // wrong name for the attribute.

  if (ihdr.color == 6 && ihdr.depth == 8) {
    uint32_t pixel;
    uint32_t zero = 0;
    for (uint32_t row = 0; row < ihdr.height; ++row) {
      for (uint32_t col = 1; col < scanline_width; col += 4) {
        size_t pos = row * scanline_width + col;

        if (original_unfiltered[pos + 3] != 0)
          continue;

        memcpy(&pixel, &original_unfiltered[pos], 4);
        invis_colors[pixel]++;

        if (normalize_alpha)
          memcpy(&original_unfiltered[pos], &zero, 4);
      }
    }
  }

  return false;

error:
  if (file != NULL)
    fclose(file);

  return true;
}

////////////////////////////////////////////////////////////////////
// Save data
////////////////////////////////////////////////////////////////////
bool PngWolf::save_original_idat(const char* path) {
  return save_idat(path, original_deflated, original_inflated);
}

bool PngWolf::save_best_idat(const char* path) {
  return save_idat(path, best_deflated, best_inflated);
}

bool PngWolf::save_idat(const char* path, std::vector<char>& deflated, std::vector<char>& inflated) {

  // TODO: when there is a simple inflate() function, make this
  // assert when deflated and inflated to not mach each other.
  // Or alternatively make some IDAT struct that has both and
  // then require its use for this function.

  FILE* out = fopen(path, "wb");

  static const uint8_t GZIP_HEADER[] = {
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x03
  };

  if (out == NULL)
    return true;

  if (fwrite(GZIP_HEADER, sizeof(GZIP_HEADER), 1, out) != 1) {
  }

  if (fwrite(&deflated[2], deflated.size() - 6, 1, out) != 1) {
  }

  // TODO: endianess?

  uint32_t crc = crc32(0L, Z_NULL, 0);
  crc = crc32(crc, (Bytef*)&inflated[0], inflated.size());

  if (fwrite(&crc, sizeof(crc), 1, out) != 1) {
  }

  uint32_t size = inflated.size();
  if (fwrite(&size, sizeof(size), 1, out) != 1) {
  }

  if (fclose(out) != 0) {
  }

  return false;
}

bool PngWolf::save_file() {

  // Create new list of chunks with old IDATs removed, and when
  // IEND is seen, insert a new IDAT chunk and the IEND chunk.
  // Then proceed with serializing the whole new list to a file.

  std::list<PngChunk> new_chunks;
  std::list<PngChunk>::iterator i;

  for (i = chunks.begin(); i != chunks.end(); ++i) {
    if (i->type == IDAT_TYPE)
      continue;

    if (i->type != IEND_TYPE) {
      new_chunks.push_back(*i);
      continue;
    }

    PngChunk new_idat;
    new_idat.type = IDAT_TYPE;
    new_idat.data = best_deflated;
    new_idat.size = new_idat.data.size();

    uint32_t idat_type_network = htonl(new_idat.type);

    new_idat.crc32 = crc32(0L, Z_NULL, 0);
    new_idat.crc32 = crc32(new_idat.crc32,
      (const Bytef*)&idat_type_network,
      sizeof(idat_type_network));

    new_idat.crc32 = crc32(new_idat.crc32,
      (const Bytef*)&new_idat.data[0],
      new_idat.data.size());

    new_chunks.push_back(new_idat);
    new_chunks.push_back(*i);
  }

  // TODO: it might be nice to not overwrite existing files, but
  // as a rule, if there are separate parameters for in and out,
  // that might not provide the best usability for users.

  FILE* out = fopen(out_path, "wb");

  if (out == NULL)
    return true;

  if (fwrite(PNG_MAGIC, 8, 1, out) != 1) {
  }

  for (i = new_chunks.begin(); i != new_chunks.end(); ++i) {
    uint32_t size = htonl(i->size);
    uint32_t type = htonl(i->type);
    uint32_t crc32 = htonl(i->crc32);

    // TODO: does this merit handling write errors?

    if (fwrite(&size, sizeof(size), 1, out) != 1) {
    }

    if (fwrite(&type, sizeof(type), 1, out) != 1) {
    }

    if (i->data.size() > 0) {
      if (fwrite(&i->data[0], i->data.size(), 1, out) != 1) {
      }
    }

    if (fwrite(&crc32, sizeof(crc32), 1, out) != 1) {
    }
  }

  if (fclose(out) != 0) {
  }

  return false;
}

////////////////////////////////////////////////////////////////////
// Help!
////////////////////////////////////////////////////////////////////
void
help(void) {
  printf("%s",
    " -----------------------------------------------------------------------------\n"
    " Usage: pngwolf --in=file.png --out=file.png                                  \n"
    " -----------------------------------------------------------------------------\n"
    "  --in=<path.png>                The PNG input image                          \n"
    "  --out=<path.png>               The PNG output file (defaults to not saving!)\n"
    "  --original-idat-to=<path.gz>   Save original IDAT data in a gzip container  \n"
    "  --best-idat-to=<path.gz>       Save best IDAT data in a gzip container      \n"
    "  --exclude-singles              Exclude single-filter genomes from population\n"
    "  --exclude-original             Exclude the filters of the input image       \n"
    "  --exclude-heuristic            Exclude the heuristically generated filters  \n"
    "  --include-experiments          Include experimental heuristic               \n"
    "  --population-size=<int>        Size of the population. Defaults to 19.      \n"
    "  --max-time=<seconds>           Timeout after seconds. (default: 0, disabled)\n"
    "  --max-stagnate-time=<seconds>  Give up if no improvement is found (d: 5)    \n"
    "  --max-deflate=<megabytes>      Give up after deflating this many megabytes  \n"
    "  --max-evaluations=<int>        Give up after evaluating this many genomes   \n"
    "  --zlib-level=<int>             zlib estimator compression level (default: 5)\n"
    "  --zlib-strategy=<int>          zlib estimator strategy (default: 0)         \n"
    "  --zlib-window=<int>            zlib estimator window bits (default: 15)     \n"
    "  --zlib-memlevel=<int>          zlib estimator memory level (default: 8)     \n"
    "  --7zip-mfb=<int>               7zip fast bytes 3..258 (default: 258)        \n"
    "  --7zip-mpass=<int>             7zip passes 0..15 (d: 2; > ~ slower, smaller)\n"
    "  --7zip-mmc=<int>               7zip match finder cycles (d: 258)            \n"
    "  --verbose-analysis             More details in initial image analysis       \n"
    "  --verbose-summary              More details in optimization summary         \n"
    "  --verbose-critter              More details when improvements are found     \n"
    "  --verbose                      Shorthand for all verbosity options          \n"
    "  --normalize-alpha              For RGBA, make fully transparent pixels black\n"
    "  --info                         Just print out verbose analysis and exit     \n"
    "  --help                         Print this help page and exit                \n"
    " -----------------------------------------------------------------------------\n"
    " To reduce the file size of PNG images `pngwolf` uses a genetic algorithm for \n"
    " finding the best scanline filter for each scanline in the image. It does not \n"
    " implement any other optimization techniques (a future version may attempt to \n"
    " use a similar approach to find a good arrangement of color palette entries). \n"
    "                                                                              \n"
    " To approximate the quality of a filter combination it compresses IDAT chunks \n"
    " using `zlib` and ultimately uses the Deflate encoder in `7-Zip` to store the \n"
    " output image. It is slow because it recompresses the IDAT data fully for all \n"
    " filter combinations even if only minor changes are made or if two filter com-\n"
    " binations are merged, as `zlib` has no built-in support for caching analysis \n"
    " data. Send mail if you know of a freely available encoder that supports that.\n"
    " -----------------------------------------------------------------------------\n"
    " Output images should be saved even if you send SIGINT (~CTRL+C) to `pngwolf`.\n"
    " The machine-readable progress report format is based on YAML http://yaml.org/\n"
    " -----------------------------------------------------------------------------\n"
    " Uses http://zlib.net/ and http://lancet.mit.edu/ga/ and http://www.7-zip.org/\n"
    " -----------------------------------------------------------------------------\n"
    " http://bjoern.hoehrmann.de/pngwolf/ (c) 2008-2011 http://bjoern.hoehrmann.de/\n"
    "");
}

////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////
int
main(int argc, char *argv[]) {

  bool argHelp = false;
  bool argVerboseAnalysis = false;
  bool argVerboseSummary = false;
  bool argVerboseCritters = false;
  bool argExcludeSingles = false;
  bool argExcludeOriginal = false;
  bool argExcludeHeuristic = false;
  bool argIncludeExperiment1 = false;
  bool argIncludeExperiment2 = false;
  bool argIncludeExperiment3 = false;
  bool argInfo = false;
  bool argNormalizeAlpha = false;
  const char* argPng = NULL;
  const char* argOut = NULL;
  const char* argBestIdatTo = NULL;
  const char* argOriginalIdatTo = NULL;
  int argMaxTime = 0;
  int argMaxStagnateTime = 5;
  int argMaxEvaluations = 0;
  int argMaxDeflate = 0;
  int argPopulationSize = 19;
  int argZlibLevel = 5;
  int argZlibStrategy = 0;
  int argZlibMemory = 8;
  int argZlibWindow = 15;
  int arg7zipFastBytes = 258;
  int arg7zipPasses = 2;
  int arg7zipCycles = 258;

  bool argOkay = true;;

#ifndef _MSC_VER
  sig_t old_handler;
#endif

  // Parse command line parameters
  for (int ax = 1; ax < argc; ++ax) {
    size_t nlen;

    const char* s = argv[ax];
    const char* value;

    // boolean options
    if (strcmp("--help", s) == 0) {
      argHelp = 1;
      break;

    } else if (strcmp("--verbose-analysis", s) == 0) {
      argVerboseAnalysis = true;
      continue;

    } else if (strcmp("--verbose-summary", s) == 0) {
      argVerboseSummary = true;
      continue;

    } else if (strcmp("--verbose-critters", s) == 0) {
      argVerboseCritters = true;
      continue;

    } else if (strcmp("--exclude-original", s) == 0) {
      argExcludeOriginal = true;
      continue;

    } else if (strcmp("--exclude-singles", s) == 0) {
      argExcludeSingles = true;
      continue;

    } else if (strcmp("--exclude-heuristic", s) == 0) {
      argExcludeHeuristic = true;
      continue;

    } else if (strcmp("--verbose", s) == 0) {
      argVerboseAnalysis = true;
      argVerboseSummary = true;
      argVerboseCritters = true;
      continue;

    } else if (strcmp("--info", s) == 0) {
      argInfo = true;
      argVerboseAnalysis = true;
      continue;

    } else if (strcmp("--normalize-alpha", s) == 0) {
      argNormalizeAlpha = true;
      continue;

    } else if (strcmp("--include-experiments", s) == 0) {
      argIncludeExperiment1 = true;
      argIncludeExperiment2 = true;
      argIncludeExperiment3 = true;
      continue;

    }

    value = strchr(s, '=');

    if (value == NULL) {
      argOkay = false;
      break;
    }

    nlen = value++ - s;

    // --name=value options
    if (strncmp("--in", s, nlen) == 0) {
      argPng = value;

    } else if (strncmp("--out", s, nlen) == 0) {
      argOut = value;

    } else if (strncmp("--best-idat-to", s, nlen) == 0) {
      argBestIdatTo = value;

    } else if (strncmp("--original-idat-to", s, nlen) == 0) {
      argOriginalIdatTo = value;

    } else if (strncmp("--max-time", s, nlen) == 0) {
      argMaxTime = atoi(value);

    } else if (strncmp("--max-stagnate-time", s, nlen) == 0) {
      argMaxStagnateTime = atoi(value);

    } else if (strncmp("--max-deflate", s, nlen) == 0) {
      argMaxDeflate = atoi(value);

    } else if (strncmp("--max-evaluations", s, nlen) == 0) {
      argMaxEvaluations = atoi(value);

    } else if (strncmp("--population-size", s, nlen) == 0) {
      argPopulationSize = atoi(value);

    } else if (strncmp("--zlib-level", s, nlen) == 0) {
      argZlibLevel = atoi(value);
      argOkay &= argZlibLevel >= 0;
      argOkay &= argZlibLevel <= 9;

    } else if (strncmp("--zlib-memory", s, nlen) == 0) {
      argZlibMemory = atoi(value);
      argOkay &= argZlibMemory >= 1;
      argOkay &= argZlibMemory <= 9;

    } else if (strncmp("--zlib-window", s, nlen) == 0) {
      argZlibWindow = atoi(value);
      argOkay &= argZlibWindow >= 8;
      argOkay &= argZlibWindow <= 15;

    } else if (strncmp("--zlib-strategy", s, nlen) == 0) {
      argZlibStrategy = atoi(value);
      argOkay &= argZlibStrategy == Z_DEFAULT_STRATEGY
              || argZlibStrategy == Z_FILTERED
              || argZlibStrategy == Z_HUFFMAN_ONLY
              || argZlibStrategy == Z_RLE;

    } else if (strncmp("--7zip-mfb", s, nlen) == 0) {
      arg7zipFastBytes = atoi(value);
      argOkay &= arg7zipFastBytes >= 3;
      argOkay &= arg7zipFastBytes <= 258;

    } else if (strncmp("--7zip-mpass", s, nlen) == 0) {
      arg7zipPasses = atoi(value);
      argOkay &= arg7zipPasses >= 1;
      argOkay &= arg7zipPasses <= 15;

    } else if (strncmp("--7zip-mmc", s, nlen) == 0) {
      arg7zipCycles = atoi(value);

    } else {
      // TODO: error
      argHelp = 1;
    }
  }

  if (argHelp || argPng == NULL || !argOkay) {
    help();
    return EXIT_SUCCESS;
  }

  wolf.in_path = argPng;
  wolf.max_deflate = argMaxDeflate;
  wolf.max_evaluations = argMaxEvaluations;
  wolf.verbose_analysis = argVerboseAnalysis;
  wolf.verbose_critters = argVerboseCritters;
  wolf.verbose_summary = argVerboseSummary;
  wolf.exclude_heuristic = argExcludeHeuristic;
  wolf.exclude_original = argExcludeOriginal;
  wolf.exclude_singles = argExcludeSingles;
  wolf.include_experiment1 = argIncludeExperiment1;
  wolf.include_experiment2 = argIncludeExperiment2;
  wolf.include_experiment3 = argIncludeExperiment3;
  wolf.population_size = argPopulationSize;
  wolf.max_stagnate_time = argMaxStagnateTime;
  wolf.last_step_at = time(NULL);
  wolf.last_improvement_at = time(NULL);
  wolf.program_begun_at = time(NULL);
  wolf.max_time = argMaxTime;
  wolf.out_path = argOut;
  wolf.zlib_level = argZlibLevel;
  wolf.zlib_memLevel = argZlibMemory;
  wolf.zlib_strategy = argZlibStrategy;
  wolf.zlib_windowBits = argZlibWindow;
  wolf.szip_fast = arg7zipFastBytes;
  wolf.szip_pass = arg7zipPasses;
  wolf.szip_cycl = arg7zipCycles;
  wolf.normalize_alpha = argNormalizeAlpha;

  // TODO: ...
  try {
  if (wolf.read_file())
    goto error;
  } catch (...) {
    goto error;
  }

  if (argOriginalIdatTo != NULL)
    if (wolf.save_original_idat(argOriginalIdatTo))
      goto out_error;

  wolf.init_filters();
  wolf.log_analysis();

  if (argInfo)
    goto done;

  wolf.search_begun_at = time(NULL);

#ifdef _MSC_VER
  SetConsoleCtrlHandler(console_event_handler, TRUE);
#else
  old_handler = signal(SIGINT, sigint_handler);
#endif

  wolf.run();

  // Uninstall the SIGINT interceptor to allow users to abort
  // the possibly very slow recompression step. This may lead
  // to users unintentionally hit CTRL+C twice, but there is
  // not much to avoid that, other than setting a timer with
  // some grace perdiod which strikes me as too complicated.

#ifdef _MSC_VER
  SetConsoleCtrlHandler(console_event_handler, FALSE);
#else
  signal(SIGINT, old_handler);
#endif

  wolf.recompress();

  if (wolf.out_path)
    if (wolf.save_file())
      goto out_error;

  if (argBestIdatTo != NULL)
    if (wolf.save_best_idat(argBestIdatTo))
      goto out_error;

  wolf.log_summary();

done:

  return EXIT_SUCCESS;

error:
  if (wolf.ihdr.interlace) {
    printf("Interlaced images are not supported\n");
    return EXIT_FAILURE;
  }

  printf("Some error occured while reading the input file\n");
  return EXIT_FAILURE;

out_error:
  printf("Some error occured while writing an output file\n");
  return EXIT_FAILURE;

}

// TODO: There are probably integer overflow issues with really,
// really big image files. Really, really big image ones. Biiig.
