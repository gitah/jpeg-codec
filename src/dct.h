#include <stdint.h>

#ifndef DCT_H
#define DCT_H

uint32_t butterfly (uint32_t i , uint32_t c);
uint32_t quadra(int16_t *in, int16_t *co);

void dct_1d(uint8_t *x, uint8_t *out);
void dct(uint8_t *in, uint8_t *out);
#endif
