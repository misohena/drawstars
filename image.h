#ifndef IMAGE_H
#define IMAGE_H

#include <cstdint>
#include <memory>
#include <utility>

#include "vmath.h"
#include "file.h"

//
// Pixel Color
//

template<typename ValueT>
struct ColorRGB
{
    using ValueType = ValueT;
    ValueType r, g, b;
    ColorRGB() : r(), g(), b(){}
    ColorRGB(ValueType r, ValueType g, ValueType b) : r(r), g(g), b(b){}
    ColorRGB operator+=(ColorRGB rhs)
    {
        r += rhs.r; g += rhs.g; b += rhs.b; return *this;
    }
    ColorRGB &operator*=(ValueType rhs)
    {
        r *= rhs; g *= rhs; b *= rhs; return *this;
    }
};
template<typename T>
ColorRGB<T> operator+(const ColorRGB<T> &lhs, const ColorRGB<T> &rhs)
{
    return ColorRGB<T>(lhs.r+rhs.r, lhs.g+rhs.g, lhs.b+rhs.b);
}
template<typename T, typename U>
ColorRGB<U> operator*(T lhs, const ColorRGB<U> &rhs)
{
    return ColorRGB<U>(lhs*rhs.r ,lhs*rhs.g, lhs*rhs.b);
}
template<typename T, typename U>
ColorRGB<U> operator*(const ColorRGB<U> &lhs, T rhs)
{
    return ColorRGB<U>(lhs.r*rhs ,lhs.g*rhs, lhs.b*rhs);
}
template<typename ValueT>
ColorRGB<ValueT> pow(const ColorRGB<ValueT> &lhs, ValueT rhs)
{
    return ColorRGB<ValueT>(
        std::pow(lhs.r, rhs),
        std::pow(lhs.g, rhs),
        std::pow(lhs.b, rhs));
}


//
// Image
//

template<typename PixelT>
class Image
{
public:
    using PixelType = PixelT;
private:
    unsigned int width_;
    unsigned int height_;
    std::unique_ptr<PixelType[]> pixels_;
public:
    Image()
        : width_(), height_(), pixels_()
    {}
    Image(unsigned int w, unsigned int h)
        : width_(w), height_(h), pixels_(new PixelType[w*h]())
    {}
    Image(Image &&rhs)
    {
        width_ = rhs.width_;
        height_ = rhs.height_;
        pixels_ = std::move(rhs.pixels_);
    }
    Image &operator=(Image &&rhs)
    {
        width_ = rhs.width_;
        height_ = rhs.height_;
        pixels_ = std::move(rhs.pixels_);
        return *this;
    }
    void reset(unsigned int w, unsigned int h)
    {
        width_ = w;
        height_ = h;
        pixels_.reset(new PixelType[w * h]);
    }
    bool empty() const{return width_ == 0 && height_ == 0;}
    explicit operator bool() const{return !empty();}
    unsigned int width() const{return width_;}
    unsigned int height() const{return height_;}
    std::size_t size() const{return width_ * height_;}
    std::size_t sizeInBytes() const{return size() * sizeof(PixelType);}
    PixelType *begin(){return pixels_.get();}
    const PixelType *begin() const {return pixels_.get();}
    PixelType *end(){return pixels_.get() + size();}
    const PixelType *end() const {return pixels_.get() + size();}
    PixelType *at(int x, int y){return pixels_.get() + (y * width_ + x);}
    const PixelType *at(int x, int y) const {return pixels_.get() + (y * width_ + x);}
    bool includes(int x, int y) const{
        return x >= 0 && y >= 0 && x < width_ && y < height_;}

    template<typename F>
    void apply(F f)
    {
        for(PixelType &p : *this){
            p = f(p);
        }
    }
};

template<typename T, typename F>
Image<T> apply(const Image<T> &srcImage, F f)
{
    Image<T> rv(srcImage.width(), srcImage.height());
    auto srcEnd = srcImage.end();
    auto src = srcImage.begin();
    auto dst = rv.begin();
    for(; src != srcEnd; ++src, ++dst){
        *dst = f(*src);
    }
    return std::move(rv);
}

template<typename T, typename F>
Image<T> &&apply(Image<T> &&image, F f)
{
    for(auto &p : image){
        p = f(p);
    }
    return std::move(image);
}

template<typename T, typename V>
Image<T> operator*(const Image<T> &lhs, V rhs)
{
    return apply(lhs, [=](const T &p){return p * rhs;});
}
template<typename T, typename V>
Image<T> &&operator*(Image<T> &&lhs, V rhs)
{
    return apply(std::move(lhs), [=](const T &p){return p * rhs;});
}


//
// Image I/O
//

template<typename ImageT>
ImageT readRawFile(const char *filename)
{
    ImageT image;
    if(File file = File(filename, "rb")){
        unsigned int w = 0, h = 0;
        std::fread(&w, sizeof(w), 1, file);
        std::fread(&h, sizeof(h), 1, file);
        if(w > 0 && h > 0 && w <= 65536 && h <= 65536){
            image.reset(w, h);
            if(std::fread(image.begin(), image.sizeInBytes(), 1, file) != 1){
                image.reset(0, 0);
            }
        }
    }
    return std::move(image);
}

template<typename ImageT>
bool writeRawFile(const char *filename, ImageT &image)
{
    if(File file = File(filename, "wb")){
        const unsigned int w = image.width();
        const unsigned int h = image.height();
        std::fwrite(&w, sizeof(w), 1, file);
        std::fwrite(&h, sizeof(h), 1, file);
        std::fwrite(image.begin(), image.sizeInBytes(), 1, file);
        return true;
    }
    else{
        return false;
    }
}

template<typename ImageT>
bool writeBMPFile(const char *filename, const ImageT &image)
{
    // Header
    const unsigned int fileHeaderSize = 14;
    const unsigned int infoHeaderSize = 40;
    const unsigned int width = image.width();
    const unsigned int height = image.height();
    const unsigned int pixelSize = 3;
    const unsigned int lineSize = pixelSize * width;
    const unsigned int lineStride = (lineSize + 3) & ~3;
    const unsigned int imageSize = lineStride * height;
    const unsigned int offBits = fileHeaderSize + infoHeaderSize;
    const unsigned int fileSize = offBits + imageSize;
#define word(v) (unsigned char)((v)&255), (unsigned char)((v)>>8&255)
#define dword(v) word(v), word(v>>16)
    const unsigned char header[] = {
        // BITMAPFILEHEADER
        'B', 'M', dword(fileSize), 0, 0, 0, 0, dword(offBits),
        // BITMAPINFOHEADER
        dword(infoHeaderSize),
        dword(width), dword(height), word(1), word(8*pixelSize),
        dword(0),
        dword(imageSize),
        dword(0), dword(0), dword(0), dword(0)
    };
#undef dword
#undef word

    // Body
    std::unique_ptr<std::uint8_t[]> body(new std::uint8_t[imageSize]);
    for(int y = 0; y < height; ++y){
        auto dst = body.get() + y * lineStride;
        auto src = image.at(0, height - y - 1);//bottom-up
        for(int x = 0; x < width; ++x, ++src, dst += 3){
            using VT = typename ImageT::PixelType::ValueType;
            dst[0] = static_cast<std::uint8_t>(clamp(255*src->b, VT(0), VT(255)));
            dst[1] = static_cast<std::uint8_t>(clamp(255*src->g, VT(0), VT(255)));
            dst[2] = static_cast<std::uint8_t>(clamp(255*src->r, VT(0), VT(255)));
        }
    }

    if(File file = File(filename, "wb")){
        std::fwrite(header, sizeof(header), 1, file);
        std::fwrite(body.get(), imageSize, 1, file);
        return true;
    }
    else{
        return false;
    }
}

#endif
