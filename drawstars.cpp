#include <cstdio>
#include <ctime>
#include <memory>
#include <thread>
#include <algorithm>
#include <map>
#include <vector>
#include <filesystem>
#include <string>
#include <string_view>
#include "vmath.h"
#include "image.h"
#include "image_jpeg.h"
#include "image_png.h"

using ColorElementType = double;
using ColorType = ColorRGB<ColorElementType>;
using ImageType = Image<ColorType>;

unsigned int numErrors = 0;
template<typename... Args>
void error(const char *fmt, Args... args)
{
    std::fprintf(stderr, fmt, args...);
    std::fprintf(stderr, "\n");
    ++numErrors;
}

float mag2flux(float deltaMag)
{
    return std::pow(10.0f, -0.4f*deltaMag);
}


//
// Image Processing
//
ColorType saturateKeepColor(ColorType c)
{
    const ColorElementType maxValue = std::max(std::max(c.r, c.g), c.b);
    if(maxValue > ColorElementType(1.0)){
        const ColorElementType invMaxValue = ColorElementType(1.0)/maxValue;
        return ColorType(
            c.r * invMaxValue,
            c.g * invMaxValue,
            c.b * invMaxValue);
    }
    else{
        return c;
    }
}
/*
ImageType &&saturateKeepColor(ImageType &&image)
{
    return apply(std::move(image), [](ColorType c){return saturateKeepColor(c);});
}

ImageType &&correctGamma(ImageType &&image, float gamma)
{
    const ColorElementType invG = 1 / gamma;
    return apply(std::move(image), [=](ColorType &p) {return pow(p, invG);});
}
*/

bool writeImageFile(const char *filename, const ImageType &image, int jpegQuality)
{
    const std::string ext = std::filesystem::path(filename).extension();
    if(ext == ".jpg" || ext == ".jpeg"){
        return writeJpegFile(filename, image, jpegQuality);
    }
    else if(ext == ".png"){
        return writePNGFile(filename, image);
    }
    else if(ext == ".bmp"){
        return writeBMPFile(filename, image);
    }
    else{
        return writeBMPFile(filename, image);
    }
}


//
// Star Catalog
//

template<typename T>
struct Range
{
    const T *begin_;
    const T *end_;
    Range(const T *b, const T *e) : begin_(b), end_(e) {}
    Range(const char *pBytes, std::size_t numObj)
        : begin_(reinterpret_cast<const T *>(pBytes))
        , end_(reinterpret_cast<const T *>(pBytes) + numObj)
    {}
    const T *begin() const {return begin_;}
    const T *end() const {return end_;}
};

struct GaiaStar
{
    // Gaia DR2
    // ra, dec : 赤経赤偉
    //         min ~ mean ~ max
    //    ra : 3.8206E-07 ~ 227.7908 ~ 360.0000
    //   dec : -89.9929 ~ -18.0945 ~ 89.9901
    // g, bp, rp : 各色のベガ等級
    //     g : 1.7076 ~ 19:2559 ~ 23:4257
    //    bp : 2.9856 ~ 19.5468 ~ 25.3263
    //    rp : 1.8629 ~ 18.0366 ~ 24.6326
    // (bp_rp: -5.4892 ~ 1.5122 ~ 9.8001)
    float ra, dec, g, bp, rp, teff_val;

    float ra_deg() const {return ra;}
    float dec_deg() const {return dec;}
    float ra_rad() const {return ra*DEG2RAD;}
    float dec_rad() const {return dec*DEG2RAD;}

    float teff() const
    {
        if(teff_val >= 100 && teff_val < 50000){ // 3229.0 <= teff_val <= 9803.0
            return teff_val;
        }
        else{
            // see: https://github.com/langurmonkey/gaiasky/blob/cfab3a75051a0d521c69b54b5f937b8dafa38b02/core/src/gaiasky/data/group/DR2DataProvider.java#L290
            const float xp = bp - rp;
            return xp <= 1.5 ?
                // Gaia broad band photometry : 5.1.G_BP−G_RP as indicator ofTeﬀ / https://arxiv.org/abs/1008.0815
                std::pow(10, 3.999 + (-0.654+ (0.709 - 0.316*xp) * xp) * xp) :
                lint(xp, 1.5f, 15.0f, 3521.6f, 3000.0f);
        }
    }
    float mag() const
    {
        return g;
    }
};
using GaiaStarRange = Range<GaiaStar>;

struct HipStar
{
    float ra, dec; //rad
    float hp, bv;

    float ra_deg() const {return ra*RAD2DEG;}
    float dec_deg() const {return dec*RAD2DEG;}
    float ra_rad() const {return ra;}
    float dec_rad() const {return dec;}

    float teff() const
    {
        // https://en.wikipedia.org/wiki/Color_index
        return 4600 * ((1 / ((0.92 * bv) + 1.7)) + (1 / ((0.92 * bv) + 0.62)) );
    }
    float mag() const
    {
#if 0
        // Hp, B-V (-0.2~1.5) to G ( https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html )
        return clamp(-0.03704 + (-0.3915 + (0.01855 - 0.03239*bv)*bv)*bv + hp, -5.0, 20.0);
#else
        return hp;
#endif
    }
};
using HipStarRange = Range<HipStar>;


//
// Star Appearance
//

/*
inline ColorType unpackColor(unsigned int rgb)
{
    return ColorType((rgb>>16&255)/255, (rgb>>8&255)/255, (rgb&255)/255);
}

inline ColorType interpolateColors(double x, double xMin, double xStepInv, unsigned int *colors, unsigned int numColors)
{
    const double indexF = (x - xMin) * xStepInv;
    const double index = std::floor(indexF);
    const double alpha = indexF - index;
    if(index < 0){
        return unpackColor(colors[0]);
    }
    else if(index + 1 < numColors){
        const std::size_t indexI = static_cast<std::size_t>(index);
        const ColorType c0 = unpackColor(colors[indexI]);
        const ColorType c1 = unpackColor(colors[indexI + 1]);
        return ColorType(
            lerp(c0.r, c1.r, alpha),
            lerp(c0.g, c1.g, alpha),
            lerp(c0.b, c1.b, alpha));
    }
    else{
        return unpackColor(colors[numColors - 1]);
    }
}
*/

// Ref:
// - How to Convert Temperature (K) to RGB: Algorithm and Sample Code | tannerhelland.com
//   https://tannerhelland.com/2012/09/18/convert-temperature-rgb-algorithm-code.html
ColorType temperatureToRGB(double k)
{
    k = clamp(k, 1000.0, 40000.0) / 100.0;
    return ColorType(
        k <= 66.0 ? 1.0 : clamp(1.29293618606274509804 * std::pow(k - 60.0, -0.1332047592), 0.0, 1.0),
        k <= 66.0 ? clamp(0.39008157876901960784 * std::log(k) - 0.63184144378862745098, 0.0, 1.0) : clamp(1.12989086089529411765 * std::pow(k - 60.0, -0.0755148492), 0.0, 1.0),
        k >= 66.0 ? 1.0 : k <= 19.0 ? 0.0 : clamp(0.54320678911019607843 * std::log(k - 10.0) - 1.19625408914, 0.0, 1.0) );
}

struct StarAppearanceOptions
{
    float fluxOffset;
    float fluxMultiplier;
    float fluxGamma;
    float fluxMax;
    float fluxIncRadius;
    float radiusDefault;
    float radiusMax;
};

struct StarAppearance
{
    ColorType color;
    float radius;
};

template<typename Star>
StarAppearance calculateStarAppearance(const Star &star, const StarAppearanceOptions &options)
{
    const float flux = mag2flux(star.mag());
    const float radius =
        flux > options.fluxIncRadius ?
        clamp(std::sqrt(options.radiusDefault * options.radiusDefault * flux / options.fluxIncRadius), 0.0f, options.radiusMax) :
        options.radiusDefault;
    const float fluxAdjusted = options.fluxOffset + options.fluxMultiplier * flux;
    const float fluxSaturated = std::min(options.fluxGamma==1.0 ? fluxAdjusted : std::pow(fluxAdjusted, options.fluxGamma), options.fluxMax);

    return {fluxSaturated * temperatureToRGB(star.teff()), radius};
}

class PixelPlotter
{
    StarAppearance appearance;
public:
    explicit PixelPlotter(const StarAppearance &a): appearance(a){}
    void operator()(ColorType *p, float plotXRad, float plotYRad, float centerXRad, float centerYRad, float dist){
        const float RADIUS_DELTA = 0.1f;
        const float d = std::sqrt(std::max(0.0f, 1.0f - dist / (appearance.radius + RADIUS_DELTA)));
        *p += d * appearance.color;
    }
};


//
// Equirectangular Projection
//

/**
 * 球面上の円を正距円筒図法で描画します。
 */
template<typename F>
void fillCircleOnEquirectangular(
    ImageType &canvas,
    float ra, float dec, float radius,
    F funcPlot)
{
    const int canvasW = canvas.width();
    const int canvasH = canvas.height();
    const float canvasW_2 = 0.5f * canvas.width(); //=canvasH?

    const float px2rad = float(PI / canvasH);
    const float x2rad = float(PI / canvasW_2);
    const float y2rad = float(PI / canvasH);
    const float rad2px = float(canvasH / PI);
    const float rad2x = float(canvasW_2 / PI);
    const float rad2y = float(canvasH / PI);

    const float centerX = canvasW * (360 - ra) / 360;
    const float centerY = canvasH * (90 - dec) / 180;

    const float centerXRad = centerX * x2rad;
    const float centerYRad = centerY * y2rad - (float(PI/2));
    const float sinCenterYRad = std::sin(centerYRad);
    const float cosCenterYRad = std::cos(centerYRad);

    const float radiusRad = radius * px2rad;
    const float cosRadiusRad = std::cos(radiusRad);

    // y方向の描画範囲は、中心から上下にradiusの範囲内
    // (正距円筒図法では経線上のピクセル数と度数の比率は常に一定なので)
    // x方向についてはy(dec:赤偉)ごとに計算する。
    const float ymin = std::max(std::floor(centerY - radius), 0.0f);
    const float ymax = std::min(std::ceil(centerY + radius), float(canvasH));
    for(float y = ymin; y < ymax; ++y){
        // [y, y+1)のラインを描画する。

        // y+0.5のライン上(緯線上)でのxの範囲を求める。
        //
        // - 円の中心の方角 : 赤経cx, 赤偉cy [rad]
        // - 描画する点の方角 : 赤経px, 赤偉py [rad]
        // - 円の中心の方角と描画する点の方角とのなす角 : r [rad]
        // とすると次の関係が成り立つ。
        // r=acos(sin(cy)*sin(py) + cos(cy)*cos(py)*cos(cx - px))
        //
        // これをpxについて解くと次のようになる。
        // px = -2*PI*n - acos(-sec(cy)*sec(py)*(sin(cy)*sin(py) - cos(r))) + cx
        // (nは任意の整数、secは1/cos)
        const float plotYRad = (y + 0.5f) * y2rad - float(PI/2);
        const float sinPlotYRad = sin(plotYRad);
        const float cosPlotYRad = cos(plotYRad);
        const float MIN_COS = 1e-6f;
        const float cosRX =
            (std::abs(cosCenterYRad) < MIN_COS || std::abs(cosPlotYRad) < MIN_COS) ? -2 :
            -(1/cosCenterYRad)*(1/cosPlotYRad) * (sinCenterYRad*sinPlotYRad - cosRadiusRad);
        const float rx =
            cosRX > 1 ? 0 : //緯線上は全て円の外
            cosRX < -1 ? canvasW_2 : // 緯線上は全て円の中
            std::acos(cosRX) * rad2px;
        //printf("rx=%f cosRX=%f \n", rx, cosRX);

        const float xmin = std::floor(centerX - rx);
        const float xmax = std::ceil(centerX + rx);
        for(float x = xmin; x < xmax; ++x){
            // [ (x, y) - (x+1, y+1) )の点を描画する(代表点は(x+0.5, y+0.5))。
            const float plotXRad = (x + 0.5f) * x2rad;
            // 中心(centerX, centerY)と(x+0.5, y+0.5)との距離を求める。
            const float dist = std::acos(sinCenterYRad*sinPlotYRad + cosCenterYRad*cosPlotYRad*std::cos(centerXRad - plotXRad)) * rad2px;
            const int ix = static_cast<int>(x + canvasW) % canvasW;// xは左右端で繰り返し処理。
            const int iy = static_cast<int>(y);
            if(canvas.includes(ix, iy)){ //念のため
                funcPlot(canvas.at(ix, iy), plotXRad, plotYRad, centerXRad, centerYRad, dist);
            }
        }
    }
}

template<typename StarType>
void drawStarOnEquirectangular(ImageType &canvas, const StarType &star, const StarAppearanceOptions &starAppearanceOptions)
{
    const StarAppearance appearance = calculateStarAppearance(star, starAppearanceOptions);

    fillCircleOnEquirectangular(
        canvas,
        star.ra_deg(), star.dec_deg(),
        appearance.radius,
        PixelPlotter(appearance));
}

template<typename Range>
void drawStarsOnEquirectangular(ImageType &canvas, Range stars, float minMag, float maxMag, const StarAppearanceOptions &starAppearanceOptions)
{
    for(const auto &star : stars){
        const auto mag = star.mag();
        if(!(mag >= minMag && mag <= maxMag)){
            continue;
        }
        drawStarOnEquirectangular(canvas, star, starAppearanceOptions);
    }
}


//
// Perspective Projection
//

Vec4F dirRADec(float ra, float dec)
{
    // x+:6h, y+: North, z+:0h(Vernal Equinox)
    return Vec4F(
        std::cos(dec) * std::sin(ra),
        std::sin(dec),
        std::cos(dec) * std::cos(ra),
        1);
}
Vec4F dirAzEl(float az, float el)
{
    // x+:East, y+:Zenith, z+:South
    return Vec4F(
        std::cos(el)*std::sin(az),
        std::sin(el),
        std::cos(el)*std::cos(az),
        1);
}

std::pair<float, float> convertRADecToAzEl(float ra, float dec, const Mat4F &matEquToHor)
{
    const Vec4F dirEquatorial = dirRADec(ra, dec);
    const Vec4F dirHorizontal = matEquToHor * dirEquatorial;
    const float x = dirHorizontal.x;//for East
    const float y = dirHorizontal.y;//for Zenith
    const float z = dirHorizontal.z;//for South
    const float az = std::atan2(-x, z);
    const float el = std::asin(y);

    return {az*RAD2DEG, el*RAD2DEG};
}

double convertUnixSecondsToGreenwichMeanSiderealTime(double unixSeconds)
{
    ///@todo validation
    ///@todo add leap seconds
    const double jd00 = unixSeconds / SEC24H + (2440587.5 - 2451545.0); //2440587.5=Unix Epoch(in JD), 2451545.0=J2000.0(in JD)
    const double t = jd00 / 36525.0; //36525.0=Days per Julian century
    const double f = SEC24H * std::fmod(jd00, 1.0);
    const double A = 24110.54841  -  SEC24H / 2.0;
    const double B = 8640184.812866;
    const double C = 0.093104;
    const double D =  -6.2e-6;
    const double gmst = ((A + (B + (C + D * t) * t) * t) + f) * SEC2RAD; //[rad]
    const double gmstNormalized = std::fmod(gmst, 2*PI);
    return gmstNormalized < 0 ? (2*PI) + gmstNormalized : gmstNormalized;
}

Mat4F createEquatorialToHorizontalMatrix(double lat, double lng, double unixSeconds)
{
    const double st = convertUnixSecondsToGreenwichMeanSiderealTime(unixSeconds) + lng*DEG2RAD;
    return Mat4F::rotX((lat-90) * DEG2RAD) * Mat4F::rotY(-st);
}

Mat4F createViewMatrix(double az, double el, double roll)
{
    return Mat4F::rotZ(-roll*DEG2RAD) *Mat4F::rotX(-el * DEG2RAD) * Mat4F::rotY((180+az) * DEG2RAD);
}

Mat4F createProjectionMatrix(double screenW, double screenH, double fovY, double viewZ)
{
    Mat4F mat = Mat4F::perspective(fovY, screenW, screenH, 0.0125, 3.0);
    if(viewZ != 0){
        mat = mat * Mat4F::translate(0, 0, -viewZ);
    }
    return mat;
}

template<typename F>
void fillCircle(ImageType &canvas, float cx, float cy, float radius, F funcPlot)
{
    const float canvasW = canvas.width();
    const float canvasH = canvas.height();
    const int xmin = static_cast<int>(clamp(std::floor(cx - radius), 0.0f, canvasW));
    const int xmax = static_cast<int>(clamp(std::ceil(cx + radius), 0.0f, canvasW));
    const int ymin = static_cast<int>(clamp(std::floor(cy - radius), 0.0f, canvasH));
    const int ymax = static_cast<int>(clamp(std::ceil(cy + radius), 0.0f, canvasH));
    for(int py = ymin; py < ymax; ++py){
        for(int px = xmin; px < xmax; ++px){
            const float pxf = px + 0.5f;
            const float pyf = py + 0.5f;
            const float dx = pxf - cx;
            const float dy = pyf - cy;
            const float dist = std::sqrt(dx*dx + dy*dy);
            funcPlot(canvas.at(px, py), pxf, pyf, cx, cy, dist);
        }
    }
}

template<typename Range>
void drawStarsPerspective(ImageType &canvas, Range stars, Mat4F mat, float minMag, float maxMag, const StarAppearanceOptions &starAppearanceOptions)
{
    for(const auto &star : stars){
        const auto mag = star.mag();
        if(!(mag >= minMag && mag <= maxMag)){
            continue;
        }

        Vec4F pos = mat * dirRADec(star.ra_rad(), star.dec_rad());
        if(pos.w != 0){
            const float invW = 1 / pos.w;
            const float x = pos.x * invW;
            const float y = pos.y * invW;
            const float z = pos.z * invW;
            if(x >= -1 && x <= 1 && y >= -1 && y <= 1 && z >= -1 && z <= 1){
                const float screenX = canvas.width()/2*(1+x);
                const float screenY = canvas.height()/2*(1-y);
                if(canvas.includes(screenX, screenY)){
                    const StarAppearance appearance = calculateStarAppearance(star, starAppearanceOptions);
                    const ColorType color = appearance.color;
                    const float radius = appearance.radius;
                    //*canvas.at(screenX, screenY) += color;
                    fillCircle(
                        canvas, screenX, screenY, radius,
                        PixelPlotter(appearance));
                }
            }
        }
    }
}


//
// Threading
//

template<typename StarType, typename F>
bool drawStars(ImageType &canvas, const char *filename, F funcDrawStars)
{
    const unsigned int canvasW = canvas.width();
    const unsigned int canvasH = canvas.height();

    File inputFile = File(filename, "rb");
    if(!inputFile){
        return false;
    }

    // prepare threads
    const unsigned int numThreads = std::thread::hardware_concurrency();
    printf("numThreads=%u\n", numThreads);

    std::unique_ptr<std::thread[]> threads(new std::thread[numThreads]);
    std::unique_ptr<ImageType[]> canvasForThreads(new ImageType[numThreads]);
    for(unsigned int ti = 0; ti < numThreads; ++ti){
        canvasForThreads[ti] = ImageType(canvasW, canvasH);
    }

    // read objects and render it
    const std::size_t OBJECT_BYTES = sizeof(StarType);
    const std::size_t BUFFER_OBJECTS = 1000000;
    const std::size_t BUFFER_BYTES = BUFFER_OBJECTS * OBJECT_BYTES;
    const std::unique_ptr<char[]> buffer_bytes(new char[BUFFER_BYTES]);
    for(;;){
        // read astronomical objects
        const std::size_t numObjectsRead = std::fread(
            buffer_bytes.get(),
            OBJECT_BYTES,
            BUFFER_OBJECTS, inputFile);

        // divide objects for threads
        const std::size_t numObjectsPerThread = (numObjectsRead + numThreads - 1) / numThreads;

        std::size_t numObjectsUnassigned = numObjectsRead;
        for(unsigned int ti = 0; ti < numThreads && numObjectsUnassigned > 0; ++ti){
            ImageType &threadCanvas = canvasForThreads[ti];
            const std::size_t threadNumObj = std::min(numObjectsUnassigned, numObjectsPerThread);
            const char * const threadHeadObj = buffer_bytes.get() + (numObjectsRead - numObjectsUnassigned) * OBJECT_BYTES;
            numObjectsUnassigned -= threadNumObj;
            if(threadNumObj > 0){
                threads[ti] = std::thread([=, &threadCanvas](){funcDrawStars(threadCanvas, Range<StarType>(threadHeadObj, threadNumObj));});
            }
        }

        for(unsigned int ti = 0; ti < numThreads; ++ti){
            if(threads[ti].joinable()){
                threads[ti].join();
            }
        }

        if(numObjectsRead < BUFFER_OBJECTS){
            break;
        }
    }

    // accumurate results of all threads
    for(unsigned int ti = 0; ti < numThreads; ++ti){
        auto srcPtr = canvasForThreads[ti].begin();
        auto dstPtr = canvas.begin();
        for(std::size_t count = canvasW * canvasH; count; --count){
            *dstPtr++ += *srcPtr++;
        }
    }
    return true;
}


//
// Command Line Options
//

// Convert string to value

template<typename RV> RV fromStrView(std::string_view str);
template<> std::string_view fromStrView<std::string_view>(std::string_view str)
{
    return str;
}
template<> int fromStrView<int>(std::string_view str)
{
    return std::stoi(std::string(str));
}
template<> float fromStrView<float>(std::string_view str)
{
    if(str.size() > 1 && str[0] == 'M'){ //M6 means Magnitude 6
        return mag2flux(std::stof(std::string(str.begin() + 1, str.end())));
    }
    else{
        return std::stof(std::string(str));
    }
}
template<> bool fromStrView<bool>(std::string_view str)
{
    if(str == "true"){
        return true;
    }
    else if(str == "false"){
        return false;
    }
    else{
        return std::stoi(std::string(str)) != 0;
    }
}

// Option Value Container

struct Options; //define by project

// Option Value Type

struct OptionTypeBase
{
    virtual ~OptionTypeBase(){}
    virtual void setArgument(Options &options, std::string_view str) const = 0;
    virtual void setDefault(Options &options) const = 0;
};
template<typename T>
struct OptionType : OptionTypeBase
{
    T Options::*dst;
    T defaultValue;
    OptionType(T Options::*dst, T defaultValue)
        : dst(dst), defaultValue(defaultValue){}
    void setArgument(Options &options, std::string_view str) const override
    {
        options.*dst = fromStrView<T>(str);
    }
    virtual void setDefault(Options &options) const override
    {
        options.*dst = defaultValue;
    }
};
struct OptionTime : OptionTypeBase
{
    double Options::*dst;
    OptionTime(double Options::*dst)
        : dst(dst){}
    void setArgument(Options &options, std::string_view str) const override
    {
        int year, month, day, hour, min, sec;
        if(sscanf(str.data(), "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &min, &sec) == 6){
            std::tm lt = {sec, min, hour, day, month-1, year-1900, 0, 0, 0};
            std::time_t unixSeconds = std::mktime(&lt);
            if(unixSeconds != -1){
                options.*dst = unixSeconds;
                return;
            }
        }
        throw std::invalid_argument("Invalid time syntax");
    }
    virtual void setDefault(Options &options) const override
    {
        options.*dst = std::time(nullptr);
    }
};
template<typename ET>
struct OptionEnum : OptionTypeBase
{
    ET Options::*dst;
    ET defaultValue;
    std::vector<std::pair<std::string_view, ET>> candidates;
    OptionEnum(ET Options::*dst, ET defaultValue, std::vector<std::pair<std::string_view, ET>> &&candidates)
        : dst(dst), defaultValue(defaultValue), candidates(std::move(candidates))
    {}
    void setArgument(Options &options, std::string_view str) const override
    {
        for(const auto &nameVal : candidates){
            if(nameVal.first == str){
                options.*dst = nameVal.second;
                return;
            }
        }
        throw std::invalid_argument("Unknown value name");
    }
    void setDefault(Options &options) const override
    {
        options.*dst = defaultValue;
    }
};
template<typename T, typename ...Args>
std::shared_ptr<OptionTypeBase> makeType(T Options::*pm, Args&&... args)
{
    if constexpr (std::is_enum_v<T>) {
        return std::make_shared<OptionEnum<T>>(pm, std::forward<Args>(args)...);
    }
    else{
        return std::make_shared<OptionType<T>>(pm, std::forward<Args>(args)...);
    }
}

// Option Definition

struct Option
{
    std::vector<std::string_view> names;
    std::shared_ptr<OptionTypeBase> type;
    bool Options::*specified = nullptr;

    bool needsValue() const
    {
        return true;
    }
};
struct OptionCategory
{
    std::string_view description;
    std::vector<Option> options;
};


enum ProjectionMode {
    PROJ_EQUIRECTANGULAR,
    PROJ_PERSPECTIVE,
    PROJ_LOAD_FLOATS //load a previous result (apply post processing only)
};
struct Options
{
    std::string_view catalogDir;
    bool gaiaEnabled;
    float gaiaMinMag;
    bool hip2Enabled;
    float hip2MaxMag;
    float minMag;
    float maxMag;

    std::string_view outputFloatsFile;
    std::string_view outputImageFile;
    int jpegQuality;

    ProjectionMode projectionMode;
    int width;
    int height;
    float fluxOffset;
    float fluxMultiplier;
    float fluxGamma;
    float fluxMax;
    float fluxIncRadius;
    float radiusDefault;
    float radiusMax;

    float postOffset;
    float postMultiplier;
    float postGamma;
    bool postKeepColor;

    double unixSeconds;
    float lat;
    float lng;
    float viewRoll;
    float viewAz;
    float viewEl;
    float viewRA;
    float viewDec;
    float fovYDeg;
    bool hasViewAz;
    bool hasViewEl;
    bool hasViewRA;
    bool hasViewDec;
    float viewZ;
};

const std::vector<OptionCategory> optionCategories = {
    {"Input", {
        {{"--catalog-dir"}, makeType(&Options::catalogDir, ".")},
        {{"--gaia"}, makeType(&Options::gaiaEnabled, true)},
        {{"--gaia-min-mag"}, makeType(&Options::gaiaMinMag, 8)},
        {{"--hip2"}, makeType(&Options::hip2Enabled, true)},
        {{"--hip2-max-mag"}, makeType(&Options::hip2MaxMag, 8)},
        {{"--min-mag"}, makeType(&Options::minMag, -100)},
        {{"--max-mag"}, makeType(&Options::maxMag, 100)},
    }},
    {"Output", {
        {{"--floats"}, makeType(&Options::outputFloatsFile, "output.floats")},
        {{"--output", "-o"}, makeType(&Options::outputImageFile, "output.png")},
        {{"--jpeg-q"}, makeType(&Options::jpegQuality, 75)},
    }},
    {"Rendering", {
        {{"--projection", "--proj"},
         makeType(
             &Options::projectionMode,
             PROJ_EQUIRECTANGULAR,
             std::vector<std::pair<std::string_view, ProjectionMode>>{
                 {"EQUIRECTANGULAR", PROJ_EQUIRECTANGULAR},
                 {"ER", PROJ_EQUIRECTANGULAR},
                 {"PERSPECTIVE", PROJ_PERSPECTIVE},
                 {"PERS", PROJ_PERSPECTIVE},
                 {"LOAD", PROJ_LOAD_FLOATS},
             })},
        {{"--width", "--w"}, makeType(&Options::width, 4096)},
        {{"--height", "--h"}, makeType(&Options::height, 2048)},

        {{"--flux-offset"}, makeType(&Options::fluxOffset, 0.0f)},
        {{"--flux-multiplier"}, makeType(&Options::fluxMultiplier, 1.0f)},
        {{"--flux-gamma"}, makeType(&Options::fluxGamma, 1.0f)},
        {{"--flux-max"}, makeType(&Options::fluxMax, 1.0f)},
        {{"--flux-inc-radius"}, makeType(&Options::fluxIncRadius, 3.98107170554e-3f)}, //mag2flux(6)
        {{"--radius-default"}, makeType(&Options::radiusDefault, 0.6f)},
        {{"--radius-max"}, makeType(&Options::radiusMax, 4.0f)},
        }},
    {"Post Process", {
        {{"--post-offset"}, makeType(&Options::postOffset, 0.0f)},
        {{"--post-multiplier"}, makeType(&Options::postMultiplier, 1.0f)},
        {{"--post-gamma"}, makeType(&Options::postGamma, 1.0f)},
        {{"--post-keep-color"}, makeType(&Options::postKeepColor, false)},
    }},
    {"Perspective Projection", {
        {{"--time"}, std::make_shared<OptionTime>(&Options::unixSeconds)},
        {{"--lat"}, makeType(&Options::lat, 35.681236f)},
        {{"--lng"}, makeType(&Options::lng, 139.767125f)},
        {{"--roll"}, makeType(&Options::viewRoll, 0.0f)},
        {{"--az"}, makeType(&Options::viewAz, 0.0f), &Options::hasViewAz},
        {{"--el", "--al"}, makeType(&Options::viewEl, 0.0f), &Options::hasViewEl},
        {{"--ra"}, makeType(&Options::viewRA, 285.0f), &Options::hasViewRA},
        {{"--dec"}, makeType(&Options::viewDec, -25.0f), &Options::hasViewDec},
        {{"--fovy"}, makeType(&Options::fovYDeg, 100.0f)},
        {{"--view-z"}, makeType(&Options::viewZ, 0.0f)},
    }},
};

void setOptionsDefault(Options &options)
{
    for(const auto &cat : optionCategories){
        for(const auto &opt : cat.options){
            opt.type->setDefault(options);
        }
    }
}
const Option *findOption(std::string_view name)
{
    for(const auto &cat : optionCategories){
        for(const auto &opt : cat.options){
            const auto it = std::find(opt.names.begin(), opt.names.end(), name);
            if(it != opt.names.end()){
                return &opt;
            }
        }
    }
    return nullptr;
}
bool emitOption(Options &options, const Option *opt, std::string_view name, std::string_view value)
{
    try{
        opt->type->setArgument(options, value);
        if(opt->specified){
            options.*(opt->specified) = true;
        }
        return true;
    }
    catch(const std::logic_error &e){ //sto? throws invalid_argument or out_of_range
        error("Invalid option value %s (%s)",
            std::string(value).c_str(),
            std::string(name).c_str());
        return false;
    }
}
bool parseCommandLine(Options &options, int argc, char *argv[])
{
    const Option *unresolvedOpt = nullptr;
    std::string_view unresolvedName;
    for(int ai = 1; ai < argc; ++ai){
        const char * const arg = argv[ai];
        if(unresolvedOpt){
            emitOption(options, unresolvedOpt, unresolvedName, std::string_view(arg));
            unresolvedOpt = nullptr;
            unresolvedName = "";
        }
        else if(arg[0] == '-'){
            const char *p = arg + 1;
            for(;!(*p == '=' || *p == '\0'); ++p);
            const std::string_view name = std::string_view(arg, p - arg);
            const Option * const opt = findOption(name);
            if(!opt){
                error("Unknown option %s", std::string(name).c_str());
                return false;
            }
            if(opt->needsValue()){
                if(*p == '='){
                    emitOption(options, opt, name, std::string_view(p+1));
                }
                else{
                    unresolvedOpt = opt;
                    unresolvedName = name;
                }
            }
            else{
                emitOption(options, opt, name, std::string_view());
            }
        }
        else{
            //files_.push_back(arg);
        }
    }
    return true;
}



//
// Main
//

int main(int argc, char *argv[])
{
    Options options = {};
    setOptionsDefault(options);
    if(!parseCommandLine(options, argc, argv)){
        return -1;
    }

    const char * const GAIA_FILE = "gaia/gaia_ra_dec_g_bp_rp_teff.dat";
    const char * const HIP_FILE = "hip2/hip2_ra_dec_mag_bv.dat";
    const ProjectionMode projectionMode = options.projectionMode;
    // Input Options
    const std::string catalogDir(options.catalogDir);
    const std::string gaiaFile = catalogDir + "/" + GAIA_FILE;
    const std::string hipFile = catalogDir + "/" + HIP_FILE;

    // Rendering Options
    const unsigned int canvasWidth = options.width;
    const unsigned int canvasHeight = options.height;
    const StarAppearanceOptions starAppearanceOptions= {
        options.fluxOffset,
        options.fluxMultiplier,
        1.0f / options.fluxGamma,
        options.fluxMax,
        options.fluxIncRadius,
        options.radiusDefault,
        options.radiusMax
    };

    // Output Options
    std::string outputFloatsFile(options.outputFloatsFile);
    std::string outputImageFile(options.outputImageFile);

    // prepare rendering target
    ImageType canvas(canvasWidth, canvasHeight);

    if(0){
        // 描画テスト
        for(int y = -900; y <= 900; y += 100){
            drawStarOnEquirectangular(canvas, HipStar{float(180*DEG2RAD), float((y+0)*PI/canvasHeight), 0, 0}, StarAppearanceOptions{1.0, 0.5, 2.0});
            drawStarOnEquirectangular(canvas, HipStar{float(185*DEG2RAD), float((y+0.25)*PI/canvasHeight), 0, 0}, StarAppearanceOptions{1.0, 0.5, 2.0});
            drawStarOnEquirectangular(canvas, HipStar{float(190*DEG2RAD), float((y+0.5)*PI/canvasHeight), 0, 0}, StarAppearanceOptions{1.0, 0.5, 2.0});
            drawStarOnEquirectangular(canvas, HipStar{float(195*DEG2RAD), float((y+0.75)*PI/canvasHeight), 0, 0}, StarAppearanceOptions{1.0, 0.5, 2.0});
            drawStarOnEquirectangular(canvas, HipStar{float(200*DEG2RAD), float((y+1.0)*PI/canvasHeight), 0, 0}, StarAppearanceOptions{1.0, 0.5, 2.0});
        }
    }
    else if(projectionMode == PROJ_LOAD_FLOATS){
        canvas = readRawFile<ImageType>(outputFloatsFile.c_str());
        if(!canvas){
            printf("Failed to load %s\n", outputFloatsFile.c_str());
            return -1;
        }
        printf("%s loaded\n", outputFloatsFile.c_str());
    }
    else if(projectionMode == PROJ_EQUIRECTANGULAR){
        if(options.gaiaEnabled && !drawStars<GaiaStar>(
               canvas, gaiaFile.c_str(),
               [=](ImageType &canvas, GaiaStarRange stars){
                   drawStarsOnEquirectangular(canvas, stars, std::max(options.gaiaMinMag, options.minMag), options.maxMag, starAppearanceOptions);
               })){
            return -1;
        }
        if(options.hip2Enabled && !drawStars<HipStar>(
               canvas, hipFile.c_str(),
               [=](ImageType &canvas, HipStarRange stars){
                   drawStarsOnEquirectangular(canvas, stars, options.minMag, std::min(options.maxMag, options.hip2MaxMag), starAppearanceOptions);
               })){
            return -1;
        }
    }
    else{
        // calculate matrix
        printf("unixSeconds=%lf\n", options.unixSeconds);
        const Mat4F matEquToHor = createEquatorialToHorizontalMatrix(options.lat, options.lng, options.unixSeconds);
        float viewAz = options.viewAz;
        float viewEl = options.viewEl;
        if((options.hasViewRA && options.hasViewDec) || !(options.hasViewAz || options.hasViewEl)){
            const float viewRA = options.viewRA;
            const float viewDec = options.viewDec;
            const std::pair<float, float> viewAzEl = convertRADecToAzEl(options.viewRA*DEG2RAD, options.viewDec*DEG2RAD, matEquToHor);
            viewAz += viewAzEl.first; //deg
            viewEl += viewAzEl.second; //deg
        }
        printf("az=%lf el=%lf roll=%lf\n", viewAz, viewEl, options.viewRoll);
        const Mat4F mat = createProjectionMatrix(canvasWidth, canvasHeight, options.fovYDeg, options.viewZ)
            * createViewMatrix(viewAz, viewEl, options.viewRoll)
            * matEquToHor;
        // draw stars
        if(options.gaiaEnabled && !drawStars<GaiaStar>(
               canvas, gaiaFile.c_str(),
               [=](ImageType &canvas, GaiaStarRange stars){
                   drawStarsPerspective(canvas, stars, mat, std::max(options.gaiaMinMag, options.minMag), options.maxMag, starAppearanceOptions);
               })){
            return -1;
        }
        if(options.hip2Enabled && !drawStars<HipStar>(
               canvas, hipFile.c_str(),
               [=](ImageType &canvas, HipStarRange stars){
                   drawStarsPerspective(canvas, stars, mat, options.minMag, std::min(options.maxMag, options.hip2MaxMag), starAppearanceOptions);
               })){
            return -1;
        }
    }

    //
    // output results
    //
    writeRawFile(outputFloatsFile.c_str(), canvas);

    const ColorType postOffsetColor(options.postOffset, options.postOffset, options.postOffset);
    const double invG = 1.0/options.postGamma;
    ImageType image = apply(
        canvas,
        [&](ColorType p){
            return pow(postOffsetColor + options.postMultiplier * p, invG);
        });
    if(options.postKeepColor){
        image.apply([&](ColorType p){return saturateKeepColor(p);});
    }
    writeImageFile(outputImageFile.c_str(), image, options.jpegQuality);

    return 0;
}
