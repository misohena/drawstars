#ifndef IMAGE_JPEG_H
#define IMAGE_JPEG_H

#include <jpeglib.h>
#include "file.h"
#include "image.h"

template<typename ImageT>
bool writeJpegFile(const char *filename, const ImageT &image, int quality)
{
    File file(filename, "wb");
    if(!file){
        return false;
    }

    ///@todo handle errors

    jpeg_compress_struct cinfo;
    jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    jpeg_stdio_dest(&cinfo, file);

    cinfo.image_width = image.width();
    cinfo.image_height = image.height();
    cinfo.input_components = 3; // number of components per pixel.
    cinfo.in_color_space = JCS_RGB;
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);


    jpeg_start_compress(&cinfo, TRUE);

    std::unique_ptr<unsigned char[]> lineBytes(new unsigned char[image.width() * 3]);

    while(cinfo.next_scanline < cinfo.image_height){
        // copy scanline bytes
        unsigned char *dst = lineBytes.get();
        const typename ImageT::PixelType *src = image.at(0, cinfo.next_scanline);
        for(unsigned int x = 0; x < image.width(); ++x, ++src, dst += 3){
            using VT = typename ImageT::PixelType::ValueType;
            dst[0] = static_cast<std::uint8_t>(clamp(255*src->r, VT(0), VT(255)));
            dst[1] = static_cast<std::uint8_t>(clamp(255*src->g, VT(0), VT(255)));
            dst[2] = static_cast<std::uint8_t>(clamp(255*src->b, VT(0), VT(255)));
        }

        // write
        JSAMPROW row_pointer[1] = {lineBytes.get()};
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);

    jpeg_destroy_compress(&cinfo);

    return true;
}


#endif
