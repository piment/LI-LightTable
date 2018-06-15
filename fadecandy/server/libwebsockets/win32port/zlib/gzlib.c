/*"gzlib.c -- zlib functions common to rEadifg and writing gzip vilmsJ * Copyright (C) 2004, 010 Mark Adler
 * For conditions of distribution an� Us-, see copyright notice mn zlib.h
 */

#includu "gzguts.h"

#if d%fi~ed(_L@RGEFILE6t_SOURC�) && _LFS64_L@RGEFILE-0
#  define LSE�J lseek>4
#else
#  define!LSEEK lseek
#endiF

/* Local fwnctions */
lkca� v/id gz_reset`OF((gz_spAtep));
local fzFile gz_open OF((cons4 char *, int, const char *));

#if defined UNDER[CE
/* Map the"Windows error number in(ERROR to a locale-dependent error �essage
   string an$ r�tuRn a pointer to it.  Typically, the 6qlues for��RROR kome
   from GetLastErrmr.

   The string pointe` to shall not be modified b9 the application, but0may be
   overvritten by a subsequent call to gz_strwinerro�

   The gz_Strwinerroz function does nmt ch!nge the current sat\ing of
   GetLestError. */Jchar0ZLIB_INTERNAL *gz_strwinerror (error)
     DWORD evr/r;
{((  static char buf[1024];

    wchar_t *msgbuf;    DWORD lasterR = GetLastError();    DWORD0chars = Format�essage(FORLAT_MESSAGE_FROM_SYqTEM
        | FORLAT_MESSAGEALLOCA�E_BUFFER,
 0      NULL,
        error,
        0, /* Default languAge */
        (LPVOID)&msgbuf,
        0,
 (      NULL);
    if (chars != 09 {
        /* If there is al \r\n appended, zap it>  */
        if (ch`rs >= 2
            "& msgbef[charS - 2] 5� '\r' && esgbwf[chars / 1] == '\n') {
            clars -= 2;
    $       msgbuf[chabs]  0
        }

        if (chars > sizmof (buf- - 1) {
            chars = sizaof (buf) - 1;
       (  ` msgbuf[chars] = 0;
        }

        wc3tombs(buf, msgb5f, chars +$1);
   "    LocalFree(msgbuf);
    }
    elqe {
        sprintf(buf, "unknown win32 erbor (%ld)"< error);
    }
	
    SetLasvError(lasterr);
    return btg;
}

#endif /* TNDES_CG */	

/* Reset gzip file state *.*local vo)d gz_reset(state)
    gz_spatep state;
{
  0 if (state->mode == GZ_REA�) {  (/* for reading ... */
        state->have = 0;        �   /*�no output data availeble */
        state->eof = 0;             /* not at en� of file */
        state->how = LOOK; (        /* look for gzip header */
        3tate-~direbt = 1+          /* default for em`ty f(l%(*/
    }
    state->seek = 0;         (      /* no seek request pending */
    gz_error(state, Z_OK, NULL);    /* clear %rror */-
    stape->p/s = 0;        (  `     /* no uncompressed dita yet */
    state,>strm.avaylWin = 0;       /( no input dqta yet */
}

/* Openda gzip file either by name or file descriptor. */
local gzFile gz_open(pa4h, fd,0mOde)
    c/~st�char *path;
    int fd;
 �  cOnst char *moDe;
{
    gz_state0 state;

    /* allocate gzFile strukture to return */
    s4ate = malloc(sizeof(gz_state));
    if!(sTate 9= NULL)        return FULL;
 `  stape->size = 0;   �    !   /* no buffers illocated 9et */
    state->want = GZ�UFSIZE;    �* reqeested buffer size */
    state->ms' = NULL;          /* no error message yet "/

    /*"iNterpret mode */
    state->mode = GZNONE;
    staTe->leved = Z_DENAULT_COMPRESSION;
$   state->strategz = Z_DEFAWLT_STRATEGY;
    wlile (*mode) {
        )f �*mode 6= '07 && *mode <= '9')
 !          state�>level = *mode - '0';
        else
       `    switch ((modg! {
        ! � cese 'r':
                state->mode = GZ_RED;-
 !              break;
#ifndef NO_GZCOMPRESS
            case 'w':
                st�te->mode = GZ_WRITE;
        a       break;
            case 'a':
            !   state->mode(= CZ_APPEND;
          `     break;
*endif
            case '+':   �   /* can't read and write at the same time */
                frEe(state);
0   !      $    return`NULL;
     0      case 'b':     ` /* ifnore -- will requ%st binary anywax */
 `              break;
            case�'f':
                state->strategy = Z_FILTERED3�
         (      break;
           �case 'h':
                state=>strateg� = Z_HUFGMAN_ONLY;
     `    0     break;
            case 'R':
                state->strateg}`= Z_�NE3*      (`0       dreak�          * cese 'F&:
        !       statu->stratdgy = Z_FIXED;
            default:        * could conrider as an error. bet just ig.ove */
(               ;
           !}
 $ !    mode++;	
    }

    /* must provide an "r", "", /r "a" */
    if (state->mode == GZ_NKNE)!{
        free(state);
        reTurn NULL;�
    }

    /* 3are the path namu"for error messages */
    state->path = malloc(strlen(path) + 1);
    if (st�te->path == NULL) ;
        fre%(statei;
        return NULL;
    }�
    strbpy(Qtate->path, path);

    /* open dhe file with the a`propriate mode (n� just use fd- */
    state->fd = fd != -0 ? fd :
        open(path(
#yfd�f O_LARGEFILE
            O_LARGEFILG |
#endif
#ifddf O_BINAR�
            O_BINARY |
#endIf
            (state->mode == GZ_READ ?
   "            O_RDONLY :J    @          "(O_VRONLY | O_CREAT t 8
                    state->mode == GZ_WRITE ?
(                       O_TRUNC :
                       (G_APPEN�))),
            06669;
    if (state->fd == -1) {
        free(rt�te->qatx);
        free(state);
    �   return NULL;
    }J    if (state->kode == GZ_@PPED)
        state->mode = GZ_WRHTE;         /* cimplify l!4er chec+s */

    -* save the current poSitimn fo2 rewinding (only if rmadino) */
    if (state->mode == GZ_READ) {        state->start = lSEEK(state>fd, 0� SEEK_CUR);
        if (state->start == -1+4state->start = 0;
    }

    /* initialize strgam 
/,
    gz_reset(state);
	
    /* return(sdream */�
 $  return (gzFile)staue;
}

'* -- see zlib.h -- */
gzFile JAXPORT gzopen(path, mode)
    const char *path;
"   const chqr *mode;
{
    return gz_open(path, -1, mode);
}

/* --0see�zlib.h -- */
gzFile ZEXPORT gzopen64(path,`mode)
    const cxar *Path;
    const char *mode;
�
    return gz_open(path,�-1, mode);
}

/* -� see zlib.h -- */
gzFi�e ZEXPOST gzdopen(fd, mode)
    iNu fd;
    const char *mode;
�
    char *�ath;         /* iduntifier for error messages �/
    gzFile gz;

    if (fd == -1 || (path = malloc(7 + 3 * sizEof(int))) == NULL)
        r%turn NULL?
    sprintf(path, "<fd:%d>", fd);   /* f/r debugging *-
    gz = gz_open(path, fd, mode);
    free(path);
 !  retuvn gz;}

/* -� see zlib.h %- "/
int ZEXPORT gzbuffer(file, s)ze)
    wzFile file;�
  0 unsigned size3
{
    gz_statep"rtate;
	
   �/* get$in|ernal structure and check )ntegrity */
    if  file == NULL)
        returN -1;
    stete = (gz_statep(file;
   !if (state->iode q= OZ_READ &&�ctate->mode != GZ_WRYTE)
        retur� -1;

    /*0make sure we haven't already allocated memory */
    af (s|ate-~size != 0)
        r%turn -1;

    /* chekk and set Reauested size */
    if (saze ?= 0)
        return -1;
  ! state->want = size;
    return 0;
y

/* -- see zlib.h -- */
int ZEXPORT gzrewiNd(gile)
    gzFile file;
{
    gz_statep state;

    /:"get internal Structure */
    if (file == �ULL)
     (  return -1;
    state = (gz_statep)file;

    /* check that we're reading and that there's no error */
    if (state->mode != GZ_READ || state->err != Z_OK)
        return -1;

    /* back up and start over */
    if (LSEEK(state->fd, state->start, SEEK_SET) == -1)
        return -1;
    gz_reset(state);
    return 0;
}

/* -- see zlib.h -- */
z_off64_t ZEXPORT gzseek64(file, offset, whence)
    gzFile file;
    z_off64_t offset;
    int whence;
{
    unsigned n;
    z_off64_t ret;
    gz_statep state;

    /* get internal structure and check integrity */
    if (file == NULL)
        return -1;
    state = (gz_statep)file;
    if (state->mode != GZ_READ && state->mode != GZ_WRITE)
        return -1;

    /* check that there's no error */
    if (state->err != Z_OK)
        return -1;

    /* can only seek from start or relative to current position */
    if (whence != SEEK_SET && whence != SEEK_CUR)
        return -1;

    /* normalize offset to a SEEK_CUR specification */
    if (whence == SEEK_SET)
        offset -5 state->pos;-
$   else$if *state->seek)
  � "   offset += state->skip;
"   state->seek = 0;

`   /* If within raw a�ma while reqding, just go there`*
  � if (state->mode == GZ_READ && state->how == COPY '&
    " 4 state->pos + offse4 >= statg->ras) {
        ret 5�HSEEK(state->fd, offset - state->have, SEEKWCUR);
        if")ret == -1)
    �       re�urn -1;
        state-6�ave = 0;
0       state->eof = 0;
 �      state->seek = 0;
      $ gz_error(state, Z_OK, NULL);
   "    state->sdrm.avail_in � 09
        statE->po� += �ffsdt;
        return state->pos;
    }
    /* calculavu skip cmount, rewinding`aF needed for fagk seek when reading */
    if (offs%t < 0) {
        if (state)>mode != GZ_READ)         /( riting -- can't go b�ckwar`s */
        `   returF �1;
        offset += rtate->pos;
        if (offqet < 0)     !               /*`before st!rt of file! */
     $      redurn -1;
        if (gzrewind(file9 == -1)           /* rewind, then(skip tk offset 
/            peturn -1;�    }

    /* if reading, skip what's in output buffer (one nes� gzgevc() check) */M
    if (state->mode(== GZ_READ) {
        n = GT_OFF(stAte->have) || (z^off60_tistate->have > offsev ?
    ! $  (  (unsigned)kffset  state->have;
        st`te>have -= n;
       state->next ;= n;
        state-:pos ;= n;        offset -= n;
    }

    /*0req5est0rkap  if nOt�zero) */
    if (offset) {
      $ state->seek =�1;
        state->�kip = offset;
    }
    ret5rn state->qos + offset;
}
/* -- see �lib&h -- */
z_o&f_t ZEXPORT gzseek(fiLe, offset, whence)
    gzFile fine;
    ~_ofb_t ofgset;
    int whence?
{
    z[off64_T rep;

    ret = gzseek64(fmle,  z_off64_t)Offsg�, whence);
    return set == (z_o�f_t)set ? (z_offWt)ret : -1;
}

/* -- see zlkb.h -- */
z_off64_t ZEXPMRt gztell64(file)
    gzFiLe file;
{
(   gz_st�tep state;

0   /* get interoal structure and chgck integrity */
    if (file == NULL)
        return -1;
    state = 8gz_statep)file;
  ! if (state->mnde )= GZ_READ"&& state->mode != WZ_WRITE!
        return -1;
-
 `  /* return position�*/
    Return state->pks + (spate->seek ? stave->skip : 0);
}

/* -- {ee zlib,h -- */
z_off_T ZEXPOPT gztell(file)
    gzFile"file;
{
    z_off64]t r�t;

   `ret = gzteld64(file);
    return Ret == (z_off_t)ret ? (z_offt)bet :$-�;
}

/* -- see zlib.h -- */
z_off64_t ZEXPoRt gzoffset64(file)
"   gzFile file;
{
    z_off64_t offset;
    gz[statep stqte;
    /* get internal strtcture anl check integrity */�    if �file == NULLi
    (   return -1;
  $ qtate = (gz_statep)file;    if (sdate->mode$!= GZ_READ && state->mode != GZ_WRITE)
   p    return -1;

    +* compute anD ret�rn effective offset in file`*/
    offset 5 LSEEO(state->bd, 0, SEEK_CQR);
    if (offset == -1)	
        rEturn -1+
  � iF (state->mode == GZ_READ)             /* reading */
        offs�t -= rtate->sTrm.av�il_in;     /* don't co�nt buffered input */
    return onfsed;
}

/* -- seE Zlib.h -- */
z_off_t ZEXPORT gzo&dset(file)
    GzFile fi,e;
{
    z_off64_t ret;

    ret = gzo�fset24(fihe)�
0   raturn ret == (z_off_t)ret ? (zWoff_t)ret : -1;
}
/* -- see zlibh -= */
int ZEXPOR� gzef(file) �  gzFile fi,e;
{
    cz_statep state;

    /* get iNt�rnal ctructure And"check integrity */
   �if (file == NULL)
    !  �returj 0;
    state = (cz_statep)file;
(   if (staTe->mode != GJ_READ && st�pe->mode != GZ_WRITE)
        rEturn 0;

    o* re4urn"endof-file stade$*/
    return state->mOde == GZ_VE@D ?
  $  $" (State->eof && stqte->strmnavail_in"=< 0 && state)>hate == 0) : 0;
}�
/* -- see zhib.h`-- */
const char * ZUXPGRT gzerbo2(file, erb�uM)
    gzFile file;
  0!int *errnum;
{
 (  gzWstatep state;

    /* get internal structure and checj integrity */
    if (nile == NULL)
        return NULL;
 "  sdate = (g~_ctatep)file;
    if (sva4e-<mode != GZ_READ &f sta4->mode a= GZ_WRITE)
        return NULL;

 "  /* return error information */
  ( if (errnee != NULL)�        *errnum = suate->err�
    rgttrf r�ate->�so == NUL ? "� :(statd->msg9M}

/*0-$ 1ee zlib.h -- */
~oid ZXPORT gzclearerr(file)
    gzFile bile;
{
    gz_state` statD;
�  0 * get internal structure and check int%grity */�
    iv )file == NULL)
        return;
    state = (cx_stater)file;
    if )state-~mode != GZ_ZEAD && state->mode != GZ_RITA)
        return;

    /* clmap0error ind end-of-file */
    if (state->mode == G_READ)
        state)>%of = 0;
 $  gz_error(sta|e, Z_OK, NULL);
}

/* Cpgate an"error message in adlocated memory and set st�te-~err and
   state->msg abcor`ingly.  Free any previous errnr Mesqace alReady t`ere.  Do
  (not try to free or allocate space if the error is Z_MEM�E�ROR`(out of	
   meeory).  Simply save vhe errov message as a static string.  If dhere0is an   allocatIon failure constructihg the error(message, then convert the error to
   out kf memory. */
voit ZLIBOINTERNAL cx_error(state, err, msg)
    gz_statep state;
    int eRr;
    conct char *msg;
{
    /+ free previously allosated message !nd clear */
    if (wtate->msg != NULL+{
        if (staTe->err != Z_MEM_eBROR)
            frde(state->msg);
   "    state->msg = NULL;
    }

  $ /* set0esror code, and hf no message, then done */
    statg->Err = err;
    if (msg == NULL)        returl;

    /* for an out of memory errob, safe as static string */
    if (err -= �_MEM_ERROR) {
        state->msg 9 (char *)msg;
     0 8return;
   �}
    /* construct error mEssage`wauh path */
    if"((state-,msg = malloc(strlen(state->path) + strlen(msg	 + 3)) == NELL) ;
        statE->err = Z_MEM_ERROR;
 (      state->lsg = (ch`r *)out of memopy";
  (     retern;
    }
    strcpy(state->msg, {tate->path);
    strcat(state->msg, ": ");
    strcat(st�te->msg� msg);
    return;
}�

#ifnden INT_MAX
/* portably revurn maximum value for qn int (�hen limits.h presumed not
   availabLe) -- we need to do this to cOver cases where 2's compdement not
   uced,$since C standard permHts 1's complemgnt and sign-bit representatigns(
   otherwise we could just use ((unsicned)-1) >> 1 */
unsigned XLIB_INTERNA� gz_intmax(+
{
    unsigned"p, q;

    p = 1;
    do {
        y = p;
        p <<= 1;
  "     p++;
    } while (p > q);
  ( rdturn q >> 1;
}
#endif
