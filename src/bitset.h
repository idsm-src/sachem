/*
 * Based on the BitSet.java code
 */
#ifndef BITSET_H_
#define BITSET_H_

#include <postgres.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "pgchem.h"


#define ADDRESS_BITS_PER_WORD   6
#define BITS_PER_WORD           (1 << ADDRESS_BITS_PER_WORD)
#define WORD_MASK               ((uint64_t) 0xffffffffffffffffL)


typedef struct
{
    uint64_t *words;
    int length;
    int wordsInUse;
} BitSet;


void inline bitset_init_alloc(BitSet *restrict const bitset, int size)
{
    int length = (size >> ADDRESS_BITS_PER_WORD) + 1;
    uint64_t *words = palloc(length * sizeof(uint64_t));

    bitset->words = words;
    bitset->length = length;
}


void inline bitset_init_setted(BitSet *restrict const bitset, int length)
{
    int endWordIndex = length >> ADDRESS_BITS_PER_WORD;
    int wordsInUse = endWordIndex + 1;
    uint64_t *words = palloc(wordsInUse * sizeof(uint64_t));

    for(int i = 0; i < endWordIndex; i++)
        words[i] = WORD_MASK;

    words[endWordIndex] = WORD_MASK >> (-((int) length) & 0x3f);

    bitset->words = words;
    bitset->length = wordsInUse;
    bitset->wordsInUse = wordsInUse;
}


void inline bitset_init_from_array(BitSet *restrict const bitset, uint8_t *data, int length)
{
    int idx = length;

    while(idx > 0 && data[idx - 1] == 0)
        idx--;

    int wordsInUse = (idx + 7) / 8;
    uint64_t *words = palloc(wordsInUse * sizeof(uint64_t));

    int i = 0;

    while(idx >= 8)
    {
        uint64_t value = 0;

        for(int i = 0; i < 8; i++)
            value |= ((uint64_t) data[i]) << (8 * i);

        words[i++] = value;

        data += 8;
        idx -= 8;
    }

    for(int j = 0; j < idx; j++)
        words[i] |= (data[j] & 0xffL) << (8 * j);

    bitset->words = words;
    bitset->length = wordsInUse;
    bitset->wordsInUse = wordsInUse;
}


void inline bitset_copy(BitSet *restrict const bitset, const BitSet *restrict const source)
{
#if BITSET_ASSERT
    if(bitset->length != source->length)
        elog(ERROR, "bitset_copy(): incompatible bitset lengths");
#endif

    memcpy(bitset->words, source->words, source->length * sizeof(uint64_t));
    bitset->wordsInUse = source->wordsInUse;
}


inline void bitset_set(BitSet *restrict const bitset, int bitIndex)
{
    int wordIndex = bitIndex >> ADDRESS_BITS_PER_WORD;

#if BITSET_ASSERT
    if(wordIndex >= bitset->length)
        elog(ERROR, "bitset_set(): wrong bitset value");
#endif

    bitset->words[wordIndex] |= (1L << (bitIndex % BITS_PER_WORD));
}


inline void bitset_unset(BitSet *restrict const bitset, int bitIndex)
{
    int wordIndex = bitIndex >> ADDRESS_BITS_PER_WORD;

#if BITSET_ASSERT
    if(wordIndex >= bitset->length)
        elog(ERROR, "bitset_unset(): wrong bitset value");
#endif

    bitset->words[wordIndex] &= ~(1L << (bitIndex % BITS_PER_WORD));

    if(bitset->wordsInUse > 0 && bitset->words[bitset->wordsInUse - 1] == 0)
        bitset->wordsInUse--;
}


inline void bitset_merge(BitSet *restrict const bitset, const BitSet *restrict const other)
{
    if(bitset == other)
        return;

    while(bitset->wordsInUse > other->wordsInUse)
        bitset->words[--bitset->wordsInUse] = 0;

    for(int i = 0; i < bitset->wordsInUse; i++)
        bitset->words[i] &= other->words[i];

    while(bitset->wordsInUse > 0 && bitset->words[bitset->wordsInUse - 1] == 0)
        bitset->wordsInUse--;
}


inline int bitset_number_of_trailing_zeros(uint64_t i)
{
    if(i == 0)
        return 64;

    uint32_t x;
    uint32_t y;
    uint32_t n = 63;

    y = (uint32_t) i;

    if(y != 0)
    {
        n = n -32;
        x = y;
    }
    else
    {
        x = (uint32_t) (i>>32);
    }

    y = x <<16;

    if(y != 0)
    {
        n = n -16;
        x = y;
    }

    y = x << 8;

    if(y != 0)
    {
        n = n - 8;
        x = y;
    }

    y = x << 4;

    if(y != 0)
    {
        n = n - 4;
        x = y;
    }

    y = x << 2;

    if(y != 0)
    {
        n = n - 2;
        x = y;
    }

    return n - ((x << 1) >> 31);
}


inline int bitset_next_set_bit(const BitSet *restrict const bitset, int fromIndex)
{
    int wordIndex = fromIndex >> ADDRESS_BITS_PER_WORD;

    if(wordIndex >= bitset->wordsInUse)
        return -1;

    uint64_t word = bitset->words[wordIndex] & (WORD_MASK << fromIndex);

    while(true)
    {
        if(word != 0)
            return (wordIndex * BITS_PER_WORD) + bitset_number_of_trailing_zeros(word);

        if(++wordIndex == bitset->wordsInUse)
            return -1;

        word = bitset->words[wordIndex];
    }
}

#endif /* BITSET_H_ */
