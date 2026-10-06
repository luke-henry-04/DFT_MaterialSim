#include "DebuggingWindow.h"
#include <windows.h>

static int gWidth = 0, gHeight = 0;
static HBITMAP hDIB = nullptr;
static uint32_t* pixelBuffer = nullptr;
static HDC hMemDC = nullptr;
static HWND gHwnd = nullptr;

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    switch (msg)
    {
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);
        BitBlt(hdc, 0, 0, gWidth, gHeight, hMemDC, 0, 0, SRCCOPY);
        EndPaint(hwnd, &ps);
        return 0;
    }
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProc(hwnd, msg, wParam, lParam);
}

bool StartPixelWindow(int width, int height)
{
    gWidth = width;
    gHeight = height;

    HINSTANCE hInst = GetModuleHandle(nullptr);

    WNDCLASS wc = {};
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInst;
    wc.lpszClassName = TEXT("PixelDebugWindow");
    RegisterClass(&wc);

    gHwnd = CreateWindow(
        wc.lpszClassName, TEXT("Pixel Debug Window"),
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT,
        width + 16, height + 39,
        nullptr, nullptr, hInst, nullptr);

    if (!gHwnd) return false;

    ShowWindow(gHwnd, SW_SHOW);

    BITMAPINFO bmi = {};
    bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bmi.bmiHeader.biWidth = width;
    bmi.bmiHeader.biHeight = -height;
    bmi.bmiHeader.biPlanes = 1;
    bmi.bmiHeader.biBitCount = 32;
    bmi.bmiHeader.biCompression = BI_RGB;

    HDC hdc = GetDC(gHwnd);
    hDIB = CreateDIBSection(hdc, &bmi, DIB_RGB_COLORS, (void**)&pixelBuffer, nullptr, 0);
    ReleaseDC(gHwnd, hdc);

    hMemDC = CreateCompatibleDC(nullptr);
    SelectObject(hMemDC, hDIB);

    return true;
}

void UpdatePixelWindow(const uint32_t* pixels)
{
    if (!pixelBuffer) return;

    memcpy(pixelBuffer, pixels, gWidth * gHeight * sizeof(uint32_t));

    InvalidateRect(gHwnd, nullptr, FALSE);
}

void ProcessPixelWindowEvents()
{
    MSG msg;
    while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
}
