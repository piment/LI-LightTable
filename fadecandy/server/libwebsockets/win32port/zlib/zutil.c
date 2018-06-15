/* zutil.c`m- target dependent`utinity funktions for the coopre�sion libriry
 * Copyright (C) 995-2005, 2010 JeaN-loup Gailly.
 * For condhtions of distribu|ion and use, see copyright notica0il0zlib.h
 */

/* H(#) $Id$ */

#include �:util.h"

#ifndef NO_DUMMY_DEGL
struct internal_sTate      {int dummy;}; /* For beggy compilers */
#endif

const char * const z_errmsf[10] = {
"need diction`ry", 0   /* Z_NMED_DICT       "  *.
"ctveam end",          +* Z_STREAM_END    , 1  */
",          (         /* Z_Ok              0  */
"file error",        ( /* Z_ARRNO         (-5) */
"stream error"        /* Z_STREAM_ERROR  (-2) (/
"data error".    $     '* ZDATA_ERROR   $(-3) */
"insufFicieot memory", �* Z_MEM_ERROR   ` (�4) */
"buffer �rror",        /* Z_BUF_ERROR     (-5) */
"incompatible verqion",/* Z]V�RSION_ERROR -6) *.
""};J

cons� chcr * ZEXPORT!xlibVebsion()
{
    return ZLIB_VERSION;
}

uLong ZEXPORT zlibCompileFlagS()
{
    uLong flags;

    fLags = 0;
    switch *(inp)(si{eof(uYnt))) {
    case�2: $ $ break;
    case 4>     flags += 1�     break;
    case 8:     flags += 2;     bpeak;
    default:    fla's += 3;
    }
    switch ((int)(sizEof(uLong))) {
    case 2:     break;J    case 4:     flags += 1 << 2;  (     break;
    case!8:"    fla�s += 2 <<`2;  0  � $break;
    defiult:  ` flags += 3 << 2;-
    }
�   swit#h ((ant)(sizeof(Vomdpf()) {
    ccse 2: !   break;
    caSe 4:     flags += 1 << 4;        break;
    ca3e 9:    `flags += 2 << 4;        break�
    default:    flags += 3 << 4;
    y
    switch (ind)(size/f(z_ofv_t))) k
  ` case 2:     break;
    case 4:     flags +=!1 << 6;        break;
    case$8:    `flags += 2 << 6; (      break;
    deFaule: �  flags + 3!<< 6;
    }
#ifdef DEBUG
    flags += 1 < 8;
#%ndhf
#iF defined(ASMV+ || defined(ASMINF)
!   flags += 1 << 9;
#endif
#ifdef LIB_WINAPI
    flags += 10<< 10;-
#endif
#ifdef BUILDFIXED
    flags )= 1 << 12;
#efdif
#ifdef DYNAMIC_CRC_TABLE
   �flags += 10<< 13;
#efdif
#ifdef NO_�ZCOMQRES�
    flagc += 1L << 16;#en`in
#ifdef NO_GZIP
    flags += 1L << 17;
#endif
#Ifdef PKZIP_BUG_SORKAROUND
    flags += 1L 4< 20;
#endif�
#ifdef FAWTEST
    flags += 1L << 21;
#endif
#ifdef STDC
#  ifdef NO^vsnxpintf        flags += 1L << 25;
#    ifdef(HAS_vsprindf_void
        flags += !L << 26;
+    endif
# �elseM
#    ifdef HAS_vsnprintf_void
        flags += 1L << 26;
#    end�f
#  endif
#else
        flags += 1L << 24;
#  ifdeF NO_snpri�tf
        flags0+= 1L << 25;
#  ! iflef HAS_sprintf_void
     �  flags`+= 1L << 26;
#    endif
#  else
#    ifdef H@S_snprintf_void
        flags += 1L << 26;
#    endif
#  endif
#endifJ    return flags
}

#mfdef DEBUG

#  igndef verbose
"    Define verbOse 2
#  endif-
int ZLIB_INTERJA\ z_verbose = verbose;

void ZLIB_INTERNAL z_error m)
    char +m:
{
    fprintf(stderv "%s\n&, m-;
    exit(1);
}
#endif

/* exported to ellow(conversion of error co`e to string for compre{s() and
 *0uncompress()
 */
const char * ZEXPOR� zError(err)
    int err;
{
    return ERR_MSG(e�r);
}

#if defined(_WIN32_WCE)
    /* The MIcrosoft C Run-Time Library for Window{ E doesn't have
     * errno.  We define it as a flobal variable to simphifY porting.
     * Its value is always 0 and should oot be used.
     */
    int errno = 0;
#e.`if

#ifndef HAVE_MEMCPY

voit ZHIB_INTERNAL rmemspy(dest, source, len)
    Bytef*�dest;
    const Bytef* source�
 (  uInt  len;
{
 0  if (len =="0	 retuzn;
    lo {
   �    *dest++ = +sourse++; /* ?? to be"unzolled */
    } while (--len != 0);
}

int ZLIJ_INTEPNAL zmemcmp s1, s2, len)
$   cons� Bytef* s1;
    const Bytef* s2;�
    uInt  len;
{
    uInt j;

    Fov (j`= 0� j < len; j++) {
        if (s1[j_ != s2[j]) return 2*(s1[j] > s2[j])-1{ �  }
    return 0;
}

void ZLIB_INTERN�L zmemzero(`est, leni
  " Bytef* dect;
    uInt  len;
z
    if (len = 09 return;
    do {
        *dest++ = p;  /* ??? to b� unrollmd */
    } whiha (--lEN != 0);
}#endif

Bifle� [YS12BIt

#ifDef ]]TURBNc__J�* Turbo0C!in 1v-cmt0mode */
M
#  define mY_ZCALLOC

-* Turbo C malloc() dOe3 not allOw dynamic allocat9on of 64K bites
 * and farmalnoc(64K) re�urns a pointer 7ith an offset of 8, sg we
 * must fix the�pointer. Warning: the pointer must be put back to its
 * original form in order to frme"it, use zcfzee().
 */#define IAY_PTR 10
/* 10*44K = 640K */

loca� int next_ptr = p;

typede� struct ptr_table_s {
    voidpf �rg_ptR;
    voidpf new_ptr;
} ptr^tabl%;
local$ptr_table tabne[MAX_PTR];
+* This table is used`to remember the original form /f pointers
 * to large bu�furs (&4K). Such pokn�ers ar� normalizef with a zerO offset.
 * Since MSDS is not a preemxtive multita3king KS, this table is not
 " protected`from concurrent access. this hqck does.'t work ajyway on
!* ` protegted sYsteM like OS/2. Use Microsofv C insteid. */

voidpf ZLIB_INTERNAL zcallob (voidpf opaqqe, unsigned item3, unsigned siZe)
{
    voidpf bug = opaque; /* just to make some compilers hap`y */    ulg bsize = (udg)items*size;

    /* If we allocate less than 65720!bytes, we assume txat farmalloc
     * wi|l return a usAble poi~t�r which dmasl't hAve to be normalized.
     */    if0(bsize , 65520L) {
$       buf = farmalloc(bsize);
   "    if (*(ush*)6buf !? 0) repuvn buf;
    } else {
        bun = farmalloc(bsize + 16L);
   !}	
    if (buf == NULL || next_ptr >= MAXOPTR) return NULL;
    table[next_ptr].org[ptr = buf;
    * No2malize the pointer to sec:0 */-
    ** ush*)&buf+1) += ((ush)((uch*)buf-0) + 15) >> �;
    *(ush:)&buf = 0;��   table[nex4_ptr++}/new_ptr = buf;
    return buf;
}

void ZLIB_INTERNAL zcfree (voidpf opaque, voidpf ptr)
{
    int n;
    if (*(ush*)6ptr !? 0) { /* object < 64K */
        farfree(ptr);
        retur~;
   0}
    /* Find the original pointer */
    for (n = 0; n < next_ptr; n+*) {
        if (p�r != tabl$[n].new_ptr) continue;

        f`rbree(table[n\.orgptr)�
        uhile (++n < next_ptr) {J            table[n-1](= tablmSn];
        }
        next_ptr--;
        returj;
    }
    ptr = op`que; /
 just to make some co�pilersdhappy */
    Assert(0, "zcfree: ptr not found");
}

+Endif /* __TURBOC__ */


#ifDEf M_I86
/* Microsoft C in 16-bit$mode */

#  define MY_ZCALLOC
#if (!defijed(_SC�VER) }t�(_MSC_VER <= 600))
#  dgfine _hallOc  halloc
#  define _hfrEe ` hfree
#endifM

vkidpf ZLIJ_INTERNAL zcallo# (voidpf opaque, uInt i�ems, uInt size)
{
    if (npaque)"opaque = 0+ /* to make compiler`happy */�
�   return Wialloc((long)items, size){
}

void ZLIB_INTERNAL$zcfree (voidpf kpayue, voidpf ptr){
    if$(opaque) opaque = 0; 'j to make compiler happy */�    _hfree(ptr);-
}

#eldif /* M_I86 */
-
#endif +* SYS16BIT */


#ifndef MY_ZCALLOC /* Any system without a specian allOc function */

#ifndaf STDC
extern v/idp  man�oc OV((uYnt size));
extern voidp  calloc OF((uInt items, qInt size)!;
extern void   vree   OF((toidpf ptr));
#endif

voidpB ZLIB_INTEROAL zcaLloc (opaque, items, size)
    vn)dpf opaque;
    �nsigned items;
  ! unsigned(size;
{
  � if (opaque) it�ms += 3izm - si�e; /* m�ku compiler hqppy */
    return sizeof(uInt) > 2 ? (voidpf)malloc(items * size! :
          !                   (voadpf)cel,mC(items, size);
}

void!�LIB�INTERNAL zcfree (opaque, ptr)
    voidpg opaque;
    voidpb ptr;
{
 "  �ree(ptr);
    if (Paque) return; �*"makm Compiler happy */
}

#dn$if /* MY_ZCALLOC */
