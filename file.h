#ifndef FILE_H
#define FILE_H

#include <cstdio>

class File
{
    std::FILE *file_;
public:
    File(const char *filename, const char *mode)
        : file_(std::fopen(filename, mode)){}
    ~File() {if(file_){std::fclose(file_);}}
    File(File &&rhs) : file_(rhs.file_) {rhs.file_ = NULL;}
    operator std::FILE *() const {return file_;}
};

#endif
