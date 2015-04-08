/*
 *  Copyright (c) 2012, Jyri J. Virkki
 *  All rights reserved.
 *
 *  This file is under BSD license. See LICENSE file.
 */

/*
 * Refer to bloom.h for documentation on the public interfaces.
 */

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "bloom.h"
#include "murmurhash2.h"


static int bloom_check_add(struct bloom * bloom,
                           const void * buffer, int len, int add)
{
  if (bloom->ready == 0) {
    (void)printf("bloom at %p not initialized!\n", (void *)bloom);
    return -1;
  }

  int hits = 0;
  register unsigned int a = murmurhash2(buffer, len, 0x9747b28c);
  register unsigned int b = murmurhash2(buffer, len, a);
  register unsigned int x;
  register unsigned int i;
  register unsigned int byte;
  register unsigned int mask;
  register unsigned char c;

  for (i = 0; i < bloom->hashes; i++) {
    x = (a + i*b) % bloom->bits;
    byte = x >> 3;
    c = bloom->bf[byte];        // expensive memory access
    mask = 1 << (x % 8);

    if (c & mask) {
      hits++;
    } else {
      if (add) {
        bloom->bf[byte] = c | mask;
      }
    }
  }

  if (hits == bloom->hashes) {
    return 1;                   // 1 == element already in (or collision)
  }

  return 0;
}


int bloom_init(struct bloom * bloom, int entries, double error)
{
  bloom->ready = 0;

  if (entries < 1 || error == 0) {
    return 1;
  }

  bloom->entries = entries;
  bloom->error = error;

  double num = log(bloom->error);
  double denom = 0.480453013918201; // ln(2)^2
  bloom->bpe = -(num / denom);

  double dentries = (double)entries;
  bloom->bits = (int)(dentries * bloom->bpe);

  if (bloom->bits % 8) {
    bloom->bytes = (bloom->bits / 8) + 1;
  } else {
    bloom->bytes = bloom->bits / 8;
  }

  bloom->hashes = (int)ceil(0.693147180559945 * bloom->bpe);  // ln(2)

  bloom->bf = (unsigned char *)calloc(bloom->bytes, sizeof(unsigned char));
  if (bloom->bf == NULL) {
    return 1;
  }

  bloom->ready = 1;
  return 0;
}


int bloom_init_fill(struct bloom * bloom, int entries, double error, unsigned char * data, int len)
{
  int init_success = bloom_init(bloom, entries, error);
  if (init_success != 0 || bloom->bytes != len) {
    return 1;
  }

  memcpy(bloom->bf, data, len);

  return 0;
}


int bloom_check(struct bloom * bloom, const void * buffer, int len)
{
  return bloom_check_add(bloom, buffer, len, 0);
}


int bloom_add(struct bloom * bloom, const void * buffer, int len)
{
  return bloom_check_add(bloom, buffer, len, 1);
}


void bloom_print(struct bloom * bloom)
{
  (void)printf("bloom at %p\n", (void *)bloom);
  (void)printf(" ->entries = %d\n", bloom->entries);
  (void)printf(" ->error = %f\n", bloom->error);
  (void)printf(" ->bits = %d\n", bloom->bits);
  (void)printf(" ->bits per elem = %f\n", bloom->bpe);
  (void)printf(" ->bytes = %d\n", bloom->bytes);
  (void)printf(" ->hash functions = %d\n", bloom->hashes);
}


void bloom_free(struct bloom * bloom)
{
  if (bloom->ready) {
    free(bloom->bf);
  }
  bloom->ready = 0;
}


double bloom_intersect_est(struct bloom * b1, struct bloom * b2) {
  if (b1->bits != b2->bits || b1->hashes != b2->hashes) {
    return -1.0;  // must be of same size and same # of hashes
  }

  double m = (float)b1->bits;
  double k = (float)b1->hashes;
  int bytes = b1->bytes;
  unsigned char *abunion = bitwise_or(b1->bf, b2->bf, bytes);

  double na = -1.0*m*log(1.0-(float)count_bits(b1->bf, bytes)/m)/k;
  double nb = -1.0*m*log(1.0-(float)count_bits(b2->bf, bytes)/m)/k;
  double nab = -1.0*m*log(1.0-(float)count_bits(abunion, bytes)/m)/k;
  return na + nb - nab;
}


int count_bits(unsigned char *cs, int bytes) {
  int i, bits;
  unsigned char c;

  (void)printf("bytes: %d\n", bytes);

  bits=0;
  for (i = 0; i < bytes; i++) {
    c = cs[i];
    // (void)printf("c: %d\n", c);
    for (; c; c >>= 1)
    {
      bits += c & 1;
    }
  }

  return bits;
}


unsigned char *bitwise_or(unsigned char * c1, unsigned char * c2, int len) {
  unsigned char *res = (unsigned char*)malloc(len);
  int i;
  for (i = 0; i < len; i++) {
    res[i] = c1[i] | c2[i];
  }
  return res;
}
