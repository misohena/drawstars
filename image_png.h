#ifndef IMAGE_PNG_H
#define IMAGE_PNG_H

#include <png.h>

template<typename ImageT>
bool writePNGFile(const char *filename, const ImageT &image)
{
    png_bytep *rows = nullptr;

    File file(filename, "wb");
    if(!file){
        return false;
    }

    const unsigned int width = image.width();
    const unsigned int height = image.height();
    const unsigned int lineStride = 3 * width;
    std::unique_ptr<std::uint8_t[]> body(new std::uint8_t[lineStride * height]);
    for(int y = 0; y < height; ++y){
        auto dst = body.get() + y * lineStride;
        auto src = image.at(0, y);
        for(int x = 0; x < width; ++x, ++src, dst += 3){
            using VT = typename ImageT::PixelType::ValueType;
            dst[0] = static_cast<std::uint8_t>(clamp(255*src->r, VT(0), VT(255)));
            dst[1] = static_cast<std::uint8_t>(clamp(255*src->g, VT(0), VT(255)));
            dst[2] = static_cast<std::uint8_t>(clamp(255*src->b, VT(0), VT(255)));
        }
    }


    const int color_type = PNG_COLOR_TYPE_RGB; //PNG_COLOR_TYPE_RGB_ALPHA, PNG_COLOR_TYPE_GRAY

    // init png objects.
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png_ptr){
        return false;
    }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr){
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        return false;
    }

    // error handler.
    if(setjmp(png_jmpbuf(png_ptr))){
        if(rows){
            png_free(png_ptr, rows);
        }
        png_destroy_write_struct(&png_ptr, &info_ptr);
        return false;
    }

    // setup i/o function.
    //png_set_write_fn(png_ptr, stream, &stream_write, &stream_flush);
    png_init_io(png_ptr, file);

    png_set_IHDR(
        png_ptr,
        info_ptr,
        width,
        height,
        8,
        color_type,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT);

    /// @todo support palette
    //png_set_PLTE(png_ptr, info_ptr, 

    rows = (png_bytep *)png_malloc(png_ptr, sizeof(png_bytep) * height);

    for(unsigned int y = 0; y < height; ++y){
        rows[y] = const_cast<png_bytep>(body.get() + y * lineStride);
    }
    png_set_rows(png_ptr, info_ptr, rows);

    const int transform = PNG_TRANSFORM_IDENTITY; //PNG_TRANSFORM_BGR;

    png_write_png(png_ptr, info_ptr, transform, NULL);

    png_free(png_ptr, rows);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    return true;
}


#endif
