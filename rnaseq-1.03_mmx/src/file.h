#ifndef __FILE_GUARD__
#define __FILE_GUARD__

#include "zlib.h"
#include "stdlib.h"
#include "stdio.h"
#include "common.h"
#include "Cmdline.h"
//typedef struct gz_stream {
//    z_stream stream;
//    int      z_err;   /* error code for last stream operation */
//    int      z_eof;   /* set if end of input file */
//    FILE     *file;   /* .gz file */
//    Byte     *inbuf;  /* input buffer */
//    Byte     *outbuf; /* output buffer */
//    uLong    crc;     /* crc32 of uncompressed data */
//    char     *msg;    /* error message */
//    char     *path;   /* path name for debugging only */
//    int      transparent; /* 1 if input file is not a .gz file */
//    char     mode;    /* 'w' or 'r' */
//    z_off_t  start;   /* start of compressed data in file (header skipped) */
//    z_off_t  in;      /* bytes into deflate or inflate */
//    z_off_t  out;     /* bytes out of deflate or inflate */
//    int      back;    /* one character push-back */
//    int      last;    /* true if push-back is last character */
//};

struct FILETYPE
{
	int TAG_COPY_LEN;
	int STRINGLENGTH;
	char NORMAL_TAGS; 
	char PAIRING_TYPE;
	char FILETYPE;
	int PAIR_LENGTH_RIGHT;
	int PAIR_LENGTH_LEFT;
	unsigned File_Size;
	FILE *Org_File;
};

FILE* File_Open(const char* File_Name,const char* Mode);
FILE* File_Exist_Open(const char* File_Name);
gzFile File_OpenZ(const char* File_Name,const char* Mode);
unsigned Get_File_Size(FILE* File);
void Detect_Input(FILETYPE & P, gzFile Input_File,gzFile Mate_File);
char Read_Tag(READ & Head,READ & Tail,gzFile Input_File,gzFile Mate_File,FILETYPE & F);
void Open_Files(gzFile & Input_File,gzFile & Mate_File,Parameters P);
#endif
