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


inline void bitset_init(BitSet *restrict const bitset, uint64_t *words, int length)
{
    bitset->words = words;
    bitset->length = length;
    bitset->wordsInUse = length;
}


inline void bitset_init_alloc(BitSet *restrict const bitset, int size)
{
    int length = (size + BITS_PER_WORD - 1) >> ADDRESS_BITS_PER_WORD;
    uint64_t *words = (uint64_t *) palloc(length * sizeof(uint64_t));

    bitset->words = words;
    bitset->length = length;
    bitset->wordsInUse = 0;
}


inline void bitset_init_empty(BitSet *restrict const bitset, int size)
{
    int length = (size + BITS_PER_WORD - 1) >> ADDRESS_BITS_PER_WORD;
    uint64_t *words = (uint64_t *) palloc0(length * sizeof(uint64_t));

    bitset->words = words;
    bitset->length = length;
    bitset->wordsInUse = 0;
}


inline void bitset_init_setted(BitSet *restrict const bitset, int size)
{
    int length = (size + BITS_PER_WORD - 1) >> ADDRESS_BITS_PER_WORD;
    uint64_t *words = (uint64_t *) palloc(length * sizeof(uint64_t));

    for(int i = 0; i < length - 1; i++)
        words[i] = WORD_MASK;

    words[length - 1] = WORD_MASK >> (-((int) size) & 0x3f);

    bitset->words = words;
    bitset->length = length;
    bitset->wordsInUse = length;
}


inline void bitset_init_from_array(BitSet *restrict const bitset, uint8_t *data, int length)
{
    int idx = length;

    while(idx > 0 && data[idx - 1] == 0)
        idx--;

    int wordsInUse = (idx + 7) / 8;
    uint64_t *words = (uint64_t *) palloc(wordsInUse * sizeof(uint64_t));

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


inline void bitset_copy(BitSet *restrict const bitset, const BitSet *restrict const source)
{
#if BITSET_ASSERT
    if(bitset->length != source->length)
        elog(ERROR, "bitset_copy(): incompatible bitset lengths");
#endif

    memcpy(bitset->words, source->words, source->length * sizeof(uint64_t));
    bitset->wordsInUse = source->wordsInUse;
}


inline bool bitset_get(BitSet *restrict const bitset, int bitIndex)
{
    int wordIndex = bitIndex >> ADDRESS_BITS_PER_WORD;

#if BITSET_ASSERT
    if(wordIndex >= bitset->length)
        elog(ERROR, "bitset_set(): wrong bitset value");
#endif

    return (bitset->words[wordIndex] & (1L << (bitIndex % BITS_PER_WORD))) != 0;
}


inline void bitset_set(BitSet *restrict const bitset, int bitIndex)
{
    int wordIndex = bitIndex >> ADDRESS_BITS_PER_WORD;

#if BITSET_ASSERT
    if(wordIndex >= bitset->length)
        elog(ERROR, "bitset_set(): wrong bitset value");
#endif

    bitset->words[wordIndex] |= (1L << (bitIndex % BITS_PER_WORD));

    if(bitset->wordsInUse <= wordIndex)
        bitset->wordsInUse = wordIndex + 1;
}


inline void bitset_unset(BitSet *restrict const bitset, int bitIndex)
{
    int wordIndex = bitIndex >> ADDRESS_BITS_PER_WORD;

#if BITSET_ASSERT
    if(wordIndex >= bitset->length)
        elog(ERROR, "bitset_unset(): wrong bitset value");
#endif

    bitset->words[wordIndex] &= ~(1L << (bitIndex % BITS_PER_WORD));

    while(bitset->wordsInUse > 0 && bitset->words[bitset->wordsInUse - 1] == 0)
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


inline int bitset_word_cardinality(uint64_t i)
{
    /* based on the Long.java code */

    i = i - ((i >> 1) & 0x5555555555555555L);
    i = (i & 0x3333333333333333L) + ((i >> 2) & 0x3333333333333333L);
    i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0fL;
    i = i + (i >> 8);
    i = i + (i >> 16);
    i = i + (i >> 32);

    return (int) i & 0x7f;
}


inline int bitset_cardinality(BitSet *restrict const bitset)
{
    int cardinality = 0;

    for(int i = 0; i < bitset->wordsInUse; i++)
        cardinality += bitset_word_cardinality(bitset->words[i]);

    return cardinality;
}


inline int bitset_and_cardinality(BitSet *restrict const bitset, const BitSet *restrict const other)
{
    int wordsInUse = bitset->wordsInUse;

    if(wordsInUse > other->wordsInUse)
        wordsInUse = other->wordsInUse;

    int cardinality = 0;

    for(int i = 0; i < wordsInUse; i++)
        cardinality += bitset_word_cardinality(bitset->words[i] & other->words[i]);

    return cardinality;
}

#endif /* BITSET_H_ */
