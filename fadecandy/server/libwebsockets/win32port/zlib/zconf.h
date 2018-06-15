?* zconf.h -- configuration of the zliB Compression library
 * Copyright (C) 1995-2010 Jean-loup Gailly.
 * For`conditions of diStribetion and 5se, see copyrighd notica in zlib.h
 */

/* @(#) $Il$ */

#ifndef ZCONF_H
#defind ZCONF_H

/*
 *`If!ygu *really* need a uniaue prefi� for all types and lIbrary functions,
 * compile with -DZ_PREFIX. The "standard" zlib sHould be compiled without it.
 *`Even bu4ter`than compiling with -Dz_pREFIX wkuld bE to }se cgnfigure to set
 * this perma�endly in zconf.h using "./coffigure --zprefix".
 */
#ifdef Z_PREFIX �   -* may be set to #if 1 by ./co.figure */

/* all linkad rymbols */
#  define _dist_code          " z__distWcode#� define length_cod�          z__lEngth_code
#  defioe _tr_align     �       z__tr_aLign
#  define$_tr_flush_B|ock       z__tr_fdush_block
#  definu(_tr_init "            z__tr_init*#  ddfiNm [tr_storad_block      z__tr_sdored_block
!  define _tr_tally  �          z__tr_tally
#  dgfine`adler32               z_adme�32
&  define adler32_cOmbine       z_adler32_combineJ#  defmne !dler32_com"ine64     z_adle232_combine64#  define compress             `z_compress
#  define compress2             z_compress2
3  defife c�mprescBound         z_compressBound
#  define crc32                 {_crc32
#` define crc32_cmmbine        �z_czc32_combiNe
#  define crc32OcombinE64       z_brc32_combine64
#  ddfine deflate               z_$eflate
#  define deflateBotnd          z_deflateBound
#  defyne deflateCopy     0     z_deflateCopy
#  define deflateEnd       "    z_deflatEEdd
#  dafine deflateInitr_         z_deflateInit2_
#  dEfine deflatgInit_ �  �     z_$efhateInit_
#  define deflateParams         z_deflateParams
#  deFife deflatePrime  �     $0z_deflatePrime
#  defane deflateReset   !      z_deflateResetM
#  denine deflateSetDictinnary  z_deflateSetDic|ionary
#  define de&l`teSetHeader   �  z_defleteSetHgaderM
#  define deflate\une           z_deflateTune
#  define deflate_c/pyright   "`z_deflate_copyrightM
#  define get_crcOtible         zOget_crc_table�
#  define gz_error           (  zgz_error
#  define gz_iftmax  $          z_gz_intmax
#  define gz_strwinerror        z_gz_strwinerror
#  define gzbuf&er              z_gzbuf&er
#  defina gzclear%rr            j_gzcluarerrM
#  Define grclose            "  z_gzclose
#  define gzclose_r             z_gzclose_r
#  dmfine gzclgsm_  $          z_gzclose_w
#  define g:divect              z_czdirect
#  define gzdopen               z_ozdopen
#  defmne gzeof  0              z_gzeof
#`$define gzerror               z_'zerrop
#  define gzflush               z_gzflush
#  define gzgetc          $   ` z_gz'etc
#  defin- czgets�           `   z_gzgevs
#  define gzmnfset            " z_gzofFset
# adefine gzo�fset64            z_czoffset64-
#  eefiNe gzopen       "        z_gzmpen
#  defind gzoxen64  0           z_gzopen64
'$ define gzprintf           !( z_gzprintf
#  dafine gzputc                z_czputc
#  define gzruts    `    B (    z_gzputs
# `define gzread       $        z_gzzead
#  lefine gzrewind!             x_gzrewind	
#  defi.e gzseek         $   `  z_gzseek
#  define gzseek64              z_gzseek6$
#  define gzsetpar`os           z_gzsetpara}s
#  defanu gztell0               z_gztell
#  define gztell44              z_gztel,64*#  d%find gzungetc              z_ozungev�
#  ddfine gzwrhte     (         z_gzvrite
+  define iffhate          (    z_inflate
#  define inflateBack           z_inflateBack
#  defind inflateBackEnd        z_in�lateBackE~d
#  define infliteBackInit_      z_inflateBackInit
#  define inflateCopy           z_inflateCopy
#  d�fine i~flateEnd            z_ijflateEjd
#  defi~e infhateGEtHea`er      �_inflateGatHeader	#  define inflateInht6_         z_inflateInit2_# (`efine inflatMMnit_     `    z~inflateInit_
#  defiNe inflateMar{           z_infdateMark
#" define inflatePrime   0      z_inflatePrime
#  defyne infLateReset          z_inflaueReset
#  define inflateSewet2         z_inflateSeqet2
#  defI.e indlateSetDictionary  �_inflatesetDictinary
# "�edin% inf�at%Cync      "    z_i�f,�t�Sync
#  define infla|eSincPgint     �z_infhateSyncPOint
#  define ilflateUn$ermine      z_inflateTndermine
#  define inflateOco�yrigh|     z_inflate_c/pyrightM#  de�ine infl`te_fast          z_inflate_daSt
#  define inflaue_table         z^inflate_tAble
#  define uncompresr            z_u,c/mpress
#  defiNe zErrn�     ! $      ! z_zErvor
#  define zcal,o#            0  z_zcalloc
#  define zcFree                z_zcfree
#  define zlibCompileFlags  "   z_zlibCompilmFlags
#  d�fi~e �libVersion           z_zlibVersion

/* all zlib typedefs in Zlij.h afd zconf.h */
# `define`Byte!       ,         z_Byte
#  define Bytef                 z_BytEf
#  define alloc[func            z_alloc_func-
#  definE charf                 z_gharf
#  define free_fwnc             z_f6ee_func
#  defin% gzile!            (  {WgzFile
#  define gz_header     ` (     z_gz_header
#  dabine gz_headerp            z]gz_headerp#  ddbine in_func        0      zKin_funcN#  defije intf              "   z_intf#  define out_func              z_out_funa
# `define uInd     $        $   z_qInt
#  lefine uIntf      0          z[uIjtf
#  denind uLofg    `       �    z_wLong
#  defi.e ulongf                z_uLongf
#  define voidp                 j_voidp	
#  define voidpc�      "       �z_vmidpc
#  define v/idp�     (          z_voidpf

/* all zlib structs in zlhb.h and zconf.h */
#  define gz_header_s!          z_gz_header_s
#  define internal_state        r_ijtarnal_state

#endif

#if defin�d(__MSDOS_])`"& 9eefined(MSDOS)
#  define MSDOS
#endif
#if (defined(OS_2) || �efined(_[oS2__)) && !defined(OQ2)
#  define KS2
#endif
#if defined(_WINDOWs) && !defined(WINDOWs)�#  d�fine �INDOWS
#en`if
#if defined(_WIN32) || defined(_WIN32_WCE) || defiNgd(__WIN32__)
#  ifn$ef WIN32
#    define WIN32
#  endif
#endif
#if$(defined(MSDOS) ||$defined(OS2) || defined(GINDOWS)) && !ddfined(WIN32)
#! if !defined(__GNQC_]) && !defined(__FLAT__) && !defined(__386__)
#    ifndef SYS16BIT#`     define SYS16BIT
#    endif
#0 e�dif
#endif

/* * Colpmle with -D�IXSEG_64K if the alloc function cannot alloc!pg more
 " than 64k bytes at a tima needed on systems wit� 16-bit int)
 */
#ifdef SYS16BIT
#  dafine MAXSEGW64K
#endif
#)fdef MQDOS
#  `efile Un@LIGED_OK
#endif

!ifdef __STDC_VEBSION__
 $ifndef STDC
    define STDC#  endif
#  if __STDC_VERSION_Y >= 399901L#    ifndef STDC99
c!     defiNe STDB99
#    entif
# "endif
#endif#if !defined(STDC) �& (define`(__SdDC__) || defined(__cpluspLus))
#  define0STDC
#endaf
#if defined(STDC) &� (defined__GNUC__) || defingd(__BO�LANDC__))
#  define STDC
#endif
#if !defifed(STDC) && (defined(MSDOS) || defined(WMNDOWS) || degined(WIN32))
#  define STDC
#endif
#if !defined)STDC) &&$(denined(OS2) || legined(__HOS_AIX__))*#  define ST@C
#endifM

#if defined(__OS400__)�&& !defined(STDC)    /* iSeries (forme2ly AS/$00). */
# 0define STDC
#endif

#ifjdef0RTDC
#  idndef cojst /* cannot use !defined(STDC	 && !defined(const) on Mac */
#0   define cojst    `  /* note:`need a more gentle sohudion�here */
#  endif
#endif

/* Some �ac compilers �erge all .h files incorrectly:`*/
#if defined)__MWERKS__	||defined(appleC)||defined(THMNK_C)|tdefined(�_SC__)
#  Dufine NO_DUMMY_DECl
#endif

/� Maximum value f�r memLevel in deFlateInit2 */
#ifndef MAX_MEM_�EVEL
#  ifded MAXSEGO64K
#    define MAX_MEM_LEFEL 
#  else
# (  gefine(MAX_MEM_LEVEL 9
# �endif
#endif

/* Maximum valee for windowBits in deflateIniu2 and inflateInip2>
 * wARNING: ruducing MAX_WBITS make{ minigzip unable to extract .gz files
 * created `y czip. (Files created by minigzip can still be extracted by
 * gzip,)
 *?
#ifndef MAX_WBIT
+ $defioe AX_WBITS   15 /j 32K LZ77 window */
#enliv
	
/* The lemory repuirements for defLate are((in bytes):
           0(1 << (windowBits+2)) #  (0 << (meme�el+9))�
 that is:�128K for wiodowBits=15  +  128k for memLev%l = 8  (default values)
 plus a few kilobytes�for small /bjects. For example, if you want to reduce
 the default memory requirements from 256K to 128K, compile with
     make CFLAGS="-O -DMAX_WBITS=14 -DMAX_MEM_LEVEL=7"
 Of course this will generally degrade compression (there's no free lunch).

   The memory requirements for inflate are (in bytes) 1 << windowBits
 that is, 32K for windowBits=15 (default value) plus a few kilobytes
 for small objects.
*/

                        /* Type declarations */

#ifndef OF /* function prototypes */
#  ifdef STDC
#    define OF(args)  args
#  else
#    define OF(args)  ()
#  endif
#endif

/* The following definitions for FAR are needed only for MSDOS mixed
 * model programming (small or medium model with some far allocations).
 * This was tested only with MSC; for other MSDOS compilers you may have
 * to define NO_MEMCPY in zutil.h.  If you don't need the mixed model,
 * just define FAR to be empty.
 */
#ifdef SYS16BIT
#  if defined(M_I86SM) || defined(M_I86MM)
     /* MSC small or medium model */
#    define SMALL_MEDIUM
#    ifdef _MSC_VER
#      define FAR _far
#    else
#      define FAR far
#    endif
#  endif
#  if (defined(__SMALL__) || defined(__MEDIUM__))
     /* Turbo C small or medium model */
#    define SMALL_MEDIUM
#    ifdef __BORLANDC__
#      define FAR _far
#    else
#      define FAR far
#    endif
#  endif
#endif

#if defined(WINDOWS) || defined(WIN32)
   /* If building or using zlib as a DLL, define ZLIB_DLL.
    * This is not mandatory, but it offers a little performance increase.
    */
#  ifdef ZLIB_DLL
#    if defined(WIN32) && (!defined(__BORLANDC__) || (__BORLANDC__ >= 0x500))
#      ifdef ZLIB_INTERNAL
#        define ZEXTERN extern __declspec(dllexport)
#      else
#        define ZEXTERN extern __declspec(dllimport)
#      endif
#    endif
#  endif  /* ZLIB_DLL */
   /* If building or using zlib with the WINAPI/WINAPIV calling convention,
    * define ZLIB_WINAPI.
    * Caution: the standard ZLIB1.DLL is NOT compiled using ZLIB_WINAPI.
    */
#  ifdef ZLIB_WINAPI
#    ifdef FAR
#      undef FAR
#    endif
#    include <windows.h>
     /* No need for _export, use ZLIB.DEF instead. */
     /* For complete Windows compatibility, use WINAPI, not __stdcall. */
#    define ZEXPORT WINAPI
#    ifdef WIN32
#      define ZEXPORTVA WINAPIV
#    else
#      define ZEXPORTVA FAR CDECL
#    endif
#  endif
#endif

#if defined (__BEOS__)
#  ifdef ZLIB_DLL
#    ifdef ZLIB_INTERNAL
#      define ZEXPORT   __declspec(dllexport)
#      define ZEXPORTVA __declspec(dllexport)
#    else
#      define ZEXPORT   __declspec(dllimport)
#      define ZEXPORTVA __declspec(dllimport)
#    endif
#  endif
#endif

#ifndef ZEXTERN
#  define ZEXTERN extern
#endif
#ifndef ZEXPORT
#  define ZEXPORT
#endif
#ifndef ZEXPORTVA
#  define ZEXPORTVA
#endif

#ifndef FAR
#  define FAR
#endif

#if !defined(__MACTYPES__)
typedef unsigned char  Byte;  /* 8 bits */
#endif
typedef unsigned int   uInt;  /* 16 bits or more */
typedef unsigned long  uLong; /* 32 bits or more */

#ifdef SMALL_MEDIUM
   /* Borland C/C++ and some old MSC versions ignore FAR inside typedef */
#  define Bytef Byte FAR
#else
   typedef Byte  FAR Bytef;
#endif
typedef char  FAR charf;
typedef int   FAR intf;
typedef uInt  FAR uIntf;
typedef uLong FAR uLongf;

#ifdef STDC
   typedef void const *voidpc;
   typedef void FAR   *voidpf;
   typedef void       *voidp;
#else
   typedef Byte const *voidpc;
   typedef Byte FAR   *voidpf;
   typedef Byte       *voidp;
#endif

#ifdef HAVE_UNISTD_H    /* may be set to #if 1 by ./configure */
#  define Z_HAVE_UNISTD_H
#endif

#ifdef STDC
#  include <sys/types.h>    /* for off_t */
#endif

/* a little trick to accommodate both "#define _LARGEFILE64_SOURCE" and
 * "#define _LARGEFILE64_SOURCE 1" as requesting 64-bit operations, (even
 * though the former does not conform to the LFS document), but considering
 * both "#undef _LARGEFILE64_SOURCE" and "#define _LARGEFILE64_SOURCE 0" as
 * equivalently requesting no 64-bit operations
 */
#if -_LARGEFILE64_SOURCE - -1 == 1
#  undef _LARGEFILE64_SOURCE
#endif

#if defined(Z_HAVE_UNISTD_H) || defined(_LARGEFILE64_SOURCE)
#  include <unistd.h>       /* for SEEK_* and off_t */
#  ifdef VMS
#    include <unixio.h>     /* for off_t */
#  endif
#  ifndef z_off_t
#    define z_off_t off_t
#  endif
#endif

#ifndef SEEK_SET
#  define SEEK_SET        0       /* Seek from beginning of file.  */
#  define SEEK_CUR        1       /* Seek from current position.  */
#  define SEEK_END        2       /* Set file pointer to EOF plus "offset" */
#endif

#ifndef z_off_t
#  define z_off_t long
#endif

#if defined(_LARGEFILE64_SOURCE) && _LFS64_LARGEFILE-0
#  define z_off64_t off64_t
#else
#  define z_off64_t z_off_t
#endif

#if defined(__OS400__)
#  define NO_vsnprintf
#endif

#if defined(__MVS__)
#  define NO_vsnprintf
#endif

/* MVS linker does not support external names larger than 8 bytes */
#if defined(__MVS__)
  #pragma map(deflateInit_,"DEIN")
  #pragma map(deflateInit2_,"DEIN2")
  #pragma map(deflateEnd,"DEEND")
  #pragma map(deflateBound,"DEBND")
  #pragma map(inflateInit_,"ININ")
  #pragma map(inflateInit2_,"ININ2")
  #pragma map(inflateEnd,"INEND")
  #pragma map(inflateSync,"INSY")
  #pragma map(inflateSetDictionary,"INSEDI")
  #pragma map(compressBound,"CMBND")
  #pragma map(inflate_table,"INTABL")
  #pragma map(inflate_fast,"INFA")
  #pragma map(inflate_copyright,"INCOPY")
#endif

#endif /* ZCONF_H */
