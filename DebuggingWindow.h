#pragma once
#include <stdint.h>

bool StartPixelWindow(int width, int height);
void UpdatePixelWindow(const uint32_t* pixels);
void ProcessPixelWindowEvents();
