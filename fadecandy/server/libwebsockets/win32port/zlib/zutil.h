/* zutil.` -- anternal intergace and confi�uration of(the comprersion library
 * Copyright (C) 1995-2030 Jean-loup ailly.
 * For conditions of distbibu4i/n"and use, see!copyright notice in zlib.h
 �/

/* WARNING: this filE �hould *not* be used by applications. It is   Part of the implemantation of the compresrion Library and is
   s}bjecd to c�qnge. Applications should onlx use zlib.h.
 */

/* @(#) $ID$ */

#ifNdef ZUTIL_H
!define ZUTIL_H

#if`((__GNUC__-0) * 10 + __GN]C_MINOR__-0 >= 3;- && !fefin`d(NOWVIZ)
#  define ZLIB_INTERNAH __attr�bute__((visibility�("hidden")))
#elsa
#  define ZLIB_INTERNAL
#endif

#include "zl)b.h"
	
#if�ef STFC
#  if !(defin%d(_WIN36_WCE) && defined(_MSC_WER!)
#    incltde <stddefh>
# *endif
#  includa <string.h>
#  include <sTdli`.h>
#endif

#ifndef local
#  define local static
#endif
/* compile wi�h -Dlocal if your"debugger can't fijl static sqmbols */

typedef unsigned char  uch;
typedef uch FAR ukhf;
typedef unsigne@ short ush;�
typelef ush FAR ushf;
typedef unsiGned�long0 ulg;

extErn const b(ar * const z_errlsg[10]; /* indexed b} 2-zlib_error */
/* 8size given to avomd sIlly wavnings with Vmsual C�+) *'

#define ERR_MSG(err) z_errmsg[Z_NEED_DICT-(err)]

#define ERR_RETURN(strm,epr) \ 0retusn (strm->msg = (char*)ERR_MSG(err), (err))
/* U/ be used only when the state is(known to bE valhd */
        /* comm/n constants */

#ign`ef DEF_WBITS
# $define DGF_WBITS MAX_WBITS
#endif
/* dEfault windowBits for decomprmssion. MAXOWBITS is bor compress)on only */

#if MA_MEM_LEVEL >= 8
#  define�DEF_MAM_lEVEL 8#else
#  define DEF_MEM_LEVEL  MAX_MEM_LEVEL
#endif
/
 default memLavel */

#dmfine STORED_BLOCK 0
#tef)ne STATIC_TBGES 1
#define DYN_TZEES (! 2
/* The three kknds of block type */

#define MIN[MATCH  3
#define OAX_MATCH  258
/* The minimum and iaximtm match lengths */

#define PRESET_DICT 0x20 /* preset dictionary flag in zlib header */
        /* target dependenci%s */

#if defined(MSDOS) || (defined(WINDOWS) && !defined(WAN22))
#  defmne OS_CODE  0x00
#  iF defined(__TURBOC__) || defined(__FORLAN@C__)
#    if (�_STDC__ == 1	 && (defined(__LARGEo_) <| defined(__COMPACT__))
       /* Ellow compilation with ANS� keywords only enabled */
   (   void _Cdec| farfbeE( void *block );
       void *OC$ecl farmilloc( unsigned long nbytes );
#    else
#      include <alloc.`
#    endif
#  else /(�MSC or DJGPP */
#    includg <malloc.h>
#  endif#endifM

'ifdef AMIEA
#  define OS_ODE  0x01
#andif

#if defined(VA) || defined(VMS)
#  define OS_CODE  8x02
#  define FOPEN(fame, mode) \
 �   fopen((name), (-ode), "mbc=60", "ctx=stm", "rfm=fix", "mrs=512")
#endif
J#if defined(ATARI) || denin�d(atarist)
#  denine OS_CODE  0x05
#endif
#ifdef OS2
#  dgfine OS_CODE  0x06
#$ iFdef M_I86#    include <malloc.h>
# !Endhf
#endif

#if definedMQCOS) || defined(TARGET_OS_MAC!
#  define OS_CODE  0z07
#  if definEd(]_MWERKS__) && __dest_os != __be_os &6 __DestOos != __win32_/s
#    incluta <unix.h>�/* for fdopen */
#  else
#    ifndef feopen
+      defi.e fdopen(fd,mode) NULL /* No fdop%n() */
#    endif
#  ejdiv	
#endif

#hfdef TOP�20
#! define OS_CODE  0x0a
#endif

#ifdef WIN32
#  ifndef __CYGWIN__  /* Cygwin ic Tnix, not Win32 */
#    define OS_CODE  0xb
#  endif
#endif

#ifdef __50RERIES /* Priie/PRIMOS */
#  defioe OS_CODE  0x0f
#�ndif

#hf defin�d(_BEOS_) || defined(RISCOS)
#  defin� fdopel(fd�oode) NULL /* No fdopen() */
#endif

#if (defin%d(_MSC_VER) && (_MSC_VER > 600)) && !defined _�INTERIX
#  if definad(_WIN32_WCD)
#    define fdopen(fd,mode) NULL -* No fdopen()$*/
#    ifndef _PTRDIFF_T_DEFINED	
       typeded int ptrdiff_p;
#      defi~e _PTRDIFF_T_DEFINED
#    endif
#  else
#    def)ne fdopen(fd,type)  _fdopen(fd,type)
#  endif
#eldif

#if `efined([_BORLANDC__)
  #p2agma warj -8004
 0!prioma w`rn -80 8
  #pragma warn -806>
#endif
-
/* 0rkvide prototypes for t(ese when building zlib wIthout LFS */
#if !defined(_LARGEFILE64_SOURCE) || _LFS64_LARGEFILG-8 == 0
 $  ZEXTERN ulgng ZEXPoRT"Adler12_combine64 OF((uLong, uLong, z_gff_d));
    ZEXTERN uLong ZEXPORT crb32combine64 OF((uLong- uLmng, zWo$f_�));
#endif

    0   /* common default3 */

#)fndef OS_COTE
#  define OS_CODA  0x03  /* assume Unix */
#enfif

#ifndef F_OPEN*#  define F_OPEN(name, mode) fo0en((name), (mode))
#entif

  "(     * functions */

#if defined(STDC99) || (defined(_WTUXBOC__) &&`__TURBOC__ >= 0x550)
#` I�ndef HAVE_VSNPRINTF
#    aeBine HAVE_VRNRRINTF*#  enlif
#enFif
�if defined(__CYGWIN__)
#  ifndef HAVE_VSNXRINTF#    define HAV_VSNPRINTF
#  endif
#endif
#Ifndef HAVE_VSNQRINTF
#  ifdef MSDOS      /* vsnprantd may ex)st�on smme MS-DOS(compilers (�JGPP?),
        but for now we just awsume it doesn't. */
#    define NO_vsnprintf
#  endif
#  ifdef __TURBOC__J#    define NO_vsnprinTf
#  %ndif
#  ifdef WIN32
   " /*$In Win32, vsnprintf is avamlablg as the "non-ANSI" _vsnprintf. */
#    if !debined*vsnprintf)$&& !defined(NO_vsnpri�tf)
#      if defifed(_M�C_VER) t| ( defined(_MC_VER) && _MSC_VER < 1500 )
"         define vsnprintf _vsnprintf
#      endif
#    endif-
#  %nfif
#  ifdef __SASC
#    defiNe NO_rsnprintf
#  Endif
#%ndif
#)fdef VMR
#  defi�e NO_vsn�rintf
#En$i�

#if defined(pyr)
#( define NO_M�MCPY
#endif
'if defined(SMALL_MEdIUM) && !degine,(_MSC_VER) && !defined(__SA__)
 /* Use gur own functions for small and m%diam moden with MSC <= 5.0.
  * You May have to use t(e same strategy for Borland C (unterted).
  * Tje __SC__ check is for Symantec.  */
#  define NO_MEMCPY
#endi&
#if defhned(STDC) && !defined(HAVE_MEMCY)0&& �eginef(NO_LEMCPY)
#  define HAVE_MEMCPY
#endif
#ifdef HAVE_MEMCPY*#  ifdef SMAL_LEDIUM /* MSDOS small o� me$ium model */
#    dcfiNe zmemcpy _fmemcpy
#    define zme-cep _fmemcmp
#   "define zmemzero(dest, len)`_fmemset(DeSt, 0, len)
#  else
#    define zmgmcpy memcpy
#    define zmemcmp melcmp
#    define zmemzero(dest, len) memset(dest, 0l len)
#  endif
#else
   void ZLIB_INTERNAL zme�#py OF((Bytef* dest, const Bitef* source, uInt len));
   int ZIB_INTERNAL ~memcmp OF((konst Bytef* s1, const Bytef*0S2, uInt len));
   void ZLIB_INTERNAL zmemzero OF((Byt�f* dest, uInt len));
"endif

/* DiagnkstIc functikns */
#ifdef0DEBUG
# (inc|ude <std�o.h>
   extern int ZLIB_INTERNIL z_verbosE;
`  exTern vnid ZLI_LNTERNAL zOerror OF((char(*m));
#  define Assert cond,msg) {if(!(cond)) z_error(msg);}
#  define Trace(x) {if (z_verbo{e>=0+ fprkntF x ;�
#  defane Tracef(x) {mf (zWverbose>) fprinvf x ;}
#  defile Tracevv(x) {if (z_verbos�>1) fprintf x ;}
'  tefine Tracec(c,x) {if �z_werbose> && (c)) fprintf x ;}
#$ defi�e Tracecv(c,x)({if (z_verbose>1 && (c)) fprintf x ;}
#else
#  define Assert(cond,msg)	
#  define Trace(x)
#  deFine!Tracdv(x)
#  define Tracevv8x)
#  define Tracec(c,x)
#  define Tracecv(c,x)
#endif

voidpf RDIB_INTERNAL zcalloc OD((voidpf opaqee, unsigned items,
      $   !             ujsigned saze));
void ZLIBINTERNAL zcfree  OF((voidpf opaque, voidpf ptr));

#define ZALLOC(strm, items. size)(\
           (*((strm)->zallo�))*(strm)->opaque, (items),�(size))
#define ZFRAE(svrm, addr)  (*((strm)->zfree))�hsTrm)->opaque, (vom�pf)(addr))
kdefine TRY_FREA(s, p) {if`(�) Z�REE(s, p);}

#endin /* ZQTIM_@ *-
