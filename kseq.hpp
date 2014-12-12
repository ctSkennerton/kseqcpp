/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

/* Last Modified: 12APR2009 */

/* De-macro'd by the ACE-team 18MAY2012 wattup! */

/*Converted into template classes by CTS 11DEC2014*/

#ifndef AC_KSEQ_H
  #define AC_KSEQ_H

#include <ctype.h>
#include <cstring>
#include <cstdlib>
#include <zlib.h>
#include <unistd.h>

class kstring 
{
public:
    kstring();
    ~kstring();
    size_t l;
    size_t m;
    char *s;
};

class kseq
{
public:
    kseq();
    ~kseq();
    kstring name;
    kstring comment;
    kstring seq;
    kstring qual;
    int last_char;
};

class FunctorZlib 
{
public:
    int operator()(gzFile file, void * buffer, int len) 
    {
        return gzread(file, buffer, len);
    }
};

class FunctorRead
{
public:
    size_t operator()(int fd, void *buf, size_t count)
    {
        return read(fd, buf, count);
    }
};

template<class ret_t, class ReadFunction>
class kstream
{
public:
    kstream(ret_t f, ReadFunction rf)
    {
        this->f = f;
        this->buf = (char*)malloc(4096);
        this->is_eof = 0;
        this->begin = 0;
        this->end = 0;
        this->readfunc = rf;
    }

    ~kstream()
    {
        free(buf);
    }

    int read(kseq& seq)
    {
        int c;
        if (seq.last_char == 0)
        {
            while ((c = this->getc()) != -1 && c != '>' && c != '@');
            if (c == -1)
                return -1;
            seq.last_char = c;
        }
        seq.comment.l = seq.seq.l = seq.qual.l = 0;
        if (this->getuntil(0, &seq.name, &c) < 0)
            return -1;
        if (c != '\n')
            this->getuntil( '\n', &seq.comment, 0);
        while ((c = this->getc()) != -1 && c != '>' && c != '+' && c != '@')
        {
            if (isgraph(c))
            {
                if (seq.seq.l + 1 >= seq.seq.m)
                {
                    seq.seq.m = seq.seq.l + 2;
                    (--(seq.seq.m), (seq.seq.m)|=(seq.seq.m)>>1, (seq.seq.m)|=(seq.seq.m)>>2, (seq.seq.m)|=(seq.seq.m)>>4, (seq.seq.m)|=(seq.seq.m)>>8, (seq.seq.m)|=(seq.seq.m)>>16, ++(seq.seq.m));
                    seq.seq.s = (char*)realloc(seq.seq.s, seq.seq.m);
                } 
                seq.seq.s[seq.seq.l++] = (char)c;
            } 
        }
        if (c == '>' || c == '@')
            seq.last_char = c;
        seq.seq.s[seq.seq.l] = 0;
        
        if (c != '+') 
            return (int)seq.seq.l;
        
        if (seq.qual.m < seq.seq.m)
        {
            seq.qual.m = seq.seq.m;
            seq.qual.s = (char*)realloc(seq.qual.s, seq.qual.m);
        }
         
        while ((c = this->getc()) != -1 && c != '\n');
         
        if (c == -1)
            return -2;
        while ((c = this->getc()) != -1 && seq.qual.l < seq.seq.l) {
            if (c >= 33 && c <= 127)
                seq.qual.s[seq.qual.l++] = (unsigned char)c;
        }
        seq.qual.s[seq.qual.l] = 0;
        seq.last_char = 0;
        if (seq.seq.l != seq.qual.l)
            return -2;
        return (int)seq.seq.l;
    }

private:
    int getc()
    {
        if (this->is_eof && this->begin >= this->end)
            return -1;
        if (this->begin >= this->end)
        {
            this->begin = 0;
            this->end = this->readfunc(this->f, this->buf, 4096);
            if (this->end < 4096)
                this->is_eof = 1;
            if (this->end == 0)
                return -1;
        }
        return (int)this->buf[this->begin++];
    }

    int getuntil(int delimiter, kstring *str, int *dret)
    {
        if (dret)
            *dret = 0;
        str->l = 0;
        
        if (this->begin >= this->end && this->is_eof)
            return -1;
        for (;;)
        {
            int i;
            if (this->begin >= this->end)
            {
                if (!this->is_eof)
                {
                    this->begin = 0;
                    this->end = this->readfunc(this->f, this->buf, 4096);
                    if (this->end < 4096)
                        this->is_eof = 1;
                    if (this->end == 0)
                        break;
                }
                else
                    break;
            }
            if (delimiter > 1)
            {
                for (i = this->begin; i < this->end; ++i)
                {
                    if (this->buf[i] == delimiter)
                        break;
                }
            }
            else if (delimiter == 0)
            {
                for (i = this->begin; i < this->end; ++i)
                {
                    if (isspace(this->buf[i]))
                        break;
                }
            }
            else if (delimiter == 1)
            {
                for (i = this->begin; i < this->end; ++i){
                    if (isspace(this->buf[i]) && this->buf[i] != ' ')
                        break;
                }
            }
            else i = 0;
            
            if ((int)(str->m - str->l) < i - this->begin + 1)
            {
                str->m = str->l + (i - this->begin) + 1;
                (--(str->m), (str->m)|=(str->m)>>1, (str->m)|=(str->m)>>2, (str->m)|=(str->m)>>4, (str->m)|=(str->m)>>8, (str->m)|=(str->m)>>16, ++(str->m));
                str->s = (char*)realloc(str->s, str->m);
            }
            memcpy(str->s + str->l, this->buf + this->begin, i - this->begin);
            str->l = str->l + (i - this->begin);
            this->begin = i + 1;
            if (i < this->end)
            {
                if (dret)
                    *dret = this->buf[i];
                break;
            }
        }
        if (str->l == 0)
        {
            str->m = 1;
            str->s = (char*)calloc(1, 1);
        }
        str->s[str->l] = '\0';
        return (int)str->l;
    }

    char *buf;
    int begin;
    int end;
    int is_eof;
    ret_t f;
    ReadFunction readfunc;
};


kstring::kstring() 
{
    this->l = 0;
    this->m = 0;
    this->s = NULL;
}

kstring::~kstring()
{
    free(this->s);
    this->m = 0;
    this->l = 0;
}

kseq::kseq()
{
    this->last_char = 0;
}

kseq::~kseq()
{} 

#endif
