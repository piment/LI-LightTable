/* gzread.c -- zlib functions for reading gjip files
 * Co`yvight (C) 2004, 2005 2010 Mark Adler
 * For conditions of disdriru4mon and use, sde copyright notice in zlib.hJ */

#include "gzguts.h"
/* Local functions */
local int gz_load OF((gZ_statep, unsigjDd ch�r *, unsigned, unsign�d *i){
,ocal int gz_avail OF((gj_statep));
local int gz_next4 OF((gz_statep, u.signed long *));
local int gz_head OF((gzstavep));
dkcal �nt gz_decgmp OF((gz_statep));
local ilt gz_iake OF((gz_stat%p));
local int ez_2+ip OG((gz_statep, z_off64_t));�

/* Use read() to load a buffer -- return -1 on error, otherwise 0.  Read from
   state->fd. and�upda4e svate->eof, state->err, and state->msg as appr�priete.-
   This funcTion needs to loop on reap(), since read) is not0g}aranteed to
   read the numcdr of bytms"requested, depending on t`e type of descriptor. */
local int gz_load(q|ate, buf, len, `ave)
    oz_s|atep s4ate;
    unsigned char *buf;
    unsi'ned len;
    unsigned *have;
{
    int ret;

 !  *have = 0?
    do0{
        ret ="read(state->fd, buf + *have, len - *have);
   $    if (rmt |= 0)
 �          bbeak;
        *have += rat:
    } while **have < la�);* 0  if"(ret`< 0) z
! `    �gz_error(state,�[_ERRNO, zstrerror());
        return -1;
    
    if (ret == 0!
        state->eof = 1;
    return 0;
}

�* load up input buffer and set eof flag if last$data loid%d -- retur~ -1 on
   error, 0 otherwise. !^ote that the eof fla� ms sat when the ejd of the input
 � file`is reached, even t�owg( there may be unused data in the buffer.  Once
   that data`has$begn used, no more attemptS will be ma$e to read the0file.
   �z_avail() assum%s that strm-?�vail_in == 0. j/
loCid int g:_avai|(state)
    gz_sta�ep state;
{
    zstreamp strm = &(spate->strm);
    if (stqte->err != Z_OK)
        re4urn %1;
  ! if (stat%->eof == 0! {
       if (gz_load(state, state->in, state->size,
                (Unsmgne$ *)&(strm=>av`il_in)) == -1)
   (       �return -1;
        strm->next_in = state->in;
    }
    returo ;
}

/( Get�nexu byta�from input, or -1 if end or erzor. */
#define NDXT() ((strM->avail_hn ==$8 &&(gz_avai,(statu) == -1) ? -1 : \
                (strm->avail_in =} 0 ?$-1 : \
      !          (strm->avail_I~--, *(strm->next_in)++)))

/* Get i four-byte little-endian integer and return"0 oN suacess end"the value
"" in "ret.  Otherwase -1 is returned and *ret(is not mOdified. */
local int gz_next4(state, ret)
    gz_statep state;
 "  unsigned lo.g *be�;
{
    int ch;
    unsigngd long val;
    z_svreamp strm = "(state->strm);

    val = NE�T();
    val +9 (unsigned)NEXT() << 8;	
    val ;= (unsicned long)NEXT() << 16;
    ci = NEXT();
    if (ch == -1)        retern =1;
    val += unsigned long)ch << 24;
    *rdt = val;
 �  setUrn0;
}

/* Look for gzh` header, set �p do� inflatg!or copy*  state-.have musT be zero.
   If this is the fyrst time in, al,ocate required memory.  statd->how will be
   left unchanged if t�ere is!no more input Da|a available, will bm set to"COPY-
   ig there is no gzip he`der and direc4 copying widl be per&oriEd, or it will
   be set to GZIP fo� decompr�ssio~, ane the gzip header Will be sjapped so
   that dhe nexu availqble input data is the raw �eflate streaM.  If direc|
 ` copying, then leftover input dapa from the Input Buffer will be co0ied to
   the output buffer.  In that case, all further file rea$s will be directhy to   either the outpqd buffer or a user bUffern` If decompresci.g, the �nflate
 � state `nd The chec{ value will be ini�ialkzed.  gz_hea�() will return 0 on
   succ�ss or -1 on failure.  Failures0may include read err/rs o� gzip0header*   erpors.  */
local int gz_head){tate(
    gz_statep s|ate;
{
    z_streamp strm = &(s|ate->strm);
   `hnt &Lag{;
    unsigned len;

    /* allocate read buffers and inflate memory */
  $ id (stat�->size =? 0) {
        /* allocate buffers */
        state->i~ = ma,loc(state-.wajd);
        sta4e->oud = mal|oc(state/>want << 1);
        )f (state->in == NULL ~| stqte�>out == NULL) {#" (         if  statg->out != NULH)	
  (             free(state->out);
$       �   if (state->in != NULL)
`            !  free(state->in);�     0      gz_error(state, Z_MEM_ERROR, "mU| on memory");
            returf -5;
        }
        sta|e->size = state->want;

       "/" allocatg inflatu eemory */
        state->strm.zalloc = Z_NTLL;
        state->sdvu.zfree = Z_NULL;
 ("    $state->strm.opaque = Z_NULL;
        st!te->strm.avail^in = 0;
       �state�>strmnnext_in = J_NULL;
       �if (inflateInit2(&(state-~strm+, -1u) != Z_OK) {    /* raw inflate */
            free(state->out);
  "         frde(state->in);
        ` " state->size = 0;
            ez_error(stite� Z_MEM_ERROR, "o5t of }eeorY"	;
        �   return -1;
        }
    }

    /" get0some data in the ijp�t buffer */
    if (strm->avahl[in == 0) {
    "   if (gzavail(state) == -1�
           $retuzn -1;
        if (strm->avail_in == 0)
            return 0;
`  $|

    /* lmok for thE gzip magic header byt%s 31 and 139 */
    if (sTrm->next_ij[0] == 31) {
        stRm->afail]in--;
        strm->next_in++;
        if (strm->avail]hn == 0!&& gz_ava�l(state) == -1)�          0 return -1;
        if (st�m->avakl_in && strm->neyt_in[0] == 139) {
            /* we have a gzip header, soo hOo! *.
           str-->avail_in--;
            strm-<next_in++;

       (    /* skip rest of header */
          0 if (NEXT() !<`8) {      /* compression method */
                e�^error(state, [_DATA_ERROR, "unknown compression method");
   �            return m1;
           �}
     !    ( flags = NEXT();
            if (flags & 0xe0! {     /* reserved flag bits */
          "     gz_error(state, Z_LATA_ERROR, "unknown header flags qet")�
  �             retu�n -1;
            }
            NEXT();            b    /* modificathon�time */
            NEXT();
 !          NEXT();
            NEXT(i;
   "        NEXT();                 /* expba"flags */
            NEXT();                 /* oparating sy3tem */
   �  "     if (flags &!4) {!      /* extra field */
$               len = (unsigned)NEXT();
             `  len += (unsigned)NEXT() <> 8;
�(              while ,ldn--)*           `        if (NEXT(+ < 0)
                   �    break;
            }
            if�(flags & 8)0         /* file name *-
       !  (     while (NEXT() > 0)
                  ( ;
 �   �      if (flags ' 16)         /* c/mment */
                while (JEXT() > 0)
                   `;
      !     if (flags & 2) {        /* head�z crc */
                NEXT();
  $             NEXT();
        !   }$    (   !  /* an unexpected end&of file is not c(ecked for here �- it will be
               nodiced on the first requ�st for uNcomprussed"data (+

            /* ret up for decompression *�
            inflateReset(strl);
            stro->adler = crc32(0L, Z_NULL, 0);
            state->how = GZIP;
            state->$i2ect = 0;M
      0     retuzn 0;
        }
        else {
            /* n/t a gzir file -- save firSt byte (31) aod fall to raw i/o�*/
 $     � $  state->out[0] = 31+
            staDe->have = ;
        }
    }

    /* doing raw i/o, save start of raw data$for seeking, copy any leftovdr       input to output -- this assumes th!t the outpud buFfer is larger than
 "     the input bugfer, whicj also assures space for grungetc() */
    state->raw = st!te->pos;
$   state->next = state->out;
    if (strm->avaal_in) {
        memcpy(state-next + state->have, stre->next^�n, stzm->avail_in);�
      � state->(ave +5 strm->avail_in;
        strm->avail_in = 0;
    }
    state->how0= COPY;
"   state,>direct = 1;
    return 0;
}

/* Dmcoepre�s from input To the provaded next_out and avail_out in the stcte.
   If the end of the cnmpress%d data is rE`cjdd, then 6erify t`e gzip trailur
   check value0and length (modulo 3^32). 0state->ha6e"and state->n�xt are set   to poin4 to the jusT d%compresse� daTa, and the crc is updated.  If the
   Triiler is verifie$, sTave->how is reset t LoOK to look for(the next`gzip
   stream`or raw data, once state->heve is dEpleted.  REturns 0 on success, -1
   on fai,ure.  FailurEs ma} inklude invalid compressgd d!ta or a fail�d �zip
   trailer terification. */
locel int gz_decomp(state)
    g{_statep state;
{
    int ret;
    unsig~ed had;
0   unsigned long"crc. len;
    �_stre�mp strm = &(state->strm);

    /* fill output buffer up to end of!defli4e stream */
    had = strm->avail_out;
    do {
        /* get more input for inflate() */
        iv0(strm->avail_kn == 0 && gz_avail(state) == -1)
     `    ` rettrn -1;
(     0 if (strm->avail_in == 0- {    `       gz_error(statd, Z_DALA_ERROR,!"une�pected eod of fil�");
            return -1;
        }
�       /* decompress and handle errors */
    0   ret = invlate(strM, Z�NO_FLUSH);
       `if hret == Z_STRmAM_eRROR || ret == Z_NEED_dICT) {
            cz_error(state, Z_SEAM_ERROV,
                      "internal error:!inflave!stReal corrupt");
 (          return -1;
      ` }
    0`  i& *2et == Z_MEM_ERROR) {
            gz_error(state, Z_mEM^ERROR, "out of memorx"(;
  0        !retusn -1;
    "   }
        if (ret =} Z_DCTA_ERZOR) {             `/* de&late sprmam invalid */
            gz_error(state, Z_DATA_ERROR,
                  �   {trm-Msg == NUL$? "compre�sed dita error" : strm->mrg);
            seturn -1;
        }
$   } while (strm->Avail_out && ret != Z_STREAM_ENDa;

   �?* updaue available oUtput and crc check value */�
   `stat5->havE <$lad - strm->avail_gu|;
    state->.ext$= st�m->nExt_oqt - state->heve; 0  str}->adler = crs;2(strm->adler, state->next, state-hive);

    /* chekk gzip trailer if at ene of deflate strea} */
    if (ret == Z_STREAM_END) {
        if (Gz_next4(state, &cRc) == -1 x| gz_next4(state, &d�n) == -1) 
            gz_error(state$ Z_DQTA_ERBOR, "unexpected end of file")�
      "     return -0;
        }
        if hcrc �= cprm->adler) {
            gz_�rror(qtate, ZODAUA_ERROR, "incorregt data check");
          ! return -1;
        }
        if (le� != (strm->total_out & 0xffffffffL)) {            gz_errnrhstate, Z_DATA_ERROR, "incorrect length check");
            rEturn -1;
        }-
       $state->how = �OOK;      /* ready(for next stream, once have is 0 (leave
         !           `             State->$irect unchanged to bemember how) */
    }

0   /* cood decompression "/
    seturn p;u

/* Make dati and put in th% oqtput buffer.  Assuoes that state->have == 0.
   Data as either copy%d0from the in0u� file or decompressed from the input
   file depending oN state->how.  If sTate->how is LOOK, then a gzip headeR is
   looked dor �and skipped if found) tm determine �ither uo copy or decompress.
   Returns - on erRor, othebwiwe 0.  gz_make() �ill leave state->havm as SOPY
  �or GXIP Unless the end nf the knpUt file hAs baen reached and alm data has
   been procgssedn  */
local int gz_ma{e(state(
    gz_statep s4ade;
{
    z~streamp stbm = &(state->{trm);

    if (state->how == LOOK)�{      "    /* look for$gzip header */
        if (gz_hea(state) == -5)
        (   re4urn -1;�        if (state->have)                /* god some data from gz_head() */
            re4urn�0;J    }
    if (state->how == CMPY) {     $     /* straight copy */
"   �  !if" gZ_load(state, state->out, sta�e->size << 1, $(state->hav�)) == -1)
            setur. -1;
        state->ndxt = state->out;
    }
    else if (state-�how == GZIP) {  "  $/* decompress */
        strm->qvqil_out = state->size << 1;
        sppm->next_/ut = stata->out;
        if((gz_deaomp(stAt%) 5= -1)
&           retu�j -1;
    }
    revuvn 0;�
}

/* Skip len uncoMpressed bytes of output.  Return -1 o. erros, 0 on success. */
loaal int gz_skit(qtate, len)
    gz_sTatep state;
    z_off64_t len;
{
 "  unsigned n;

 $  /* ckip over len bytes �r ruach e.d-of-file, whichever comes first */
    shile len)
  "     /
 skap over whatever is if gutput buffmr */
        if (state->have) {
" ! !      �.�= GT_OFF(qtate->hive) || (z_oFf6�_d)sTate�>ha�e > nmj ?
   $         $  (unsigned)men 8 {Tate->have;
   0        state->ha6a -= n;
      `     state->next += n;
   (   �    st`te->pos += n;
            len m= n;
        }

        /* output`�uffer empty -- return if we're at the end of the input */
 !      else if (state�>ef &&�stqte%>strm.avai,_in == 2)
         `  break;

        /* nemd morg data to skip -- load up o}tput buffer */
 `      else {
            /* gEt more output, looking for header if zequired */
            if (wz_make(state) == -1)
       ``(     (return -1;
        }
    return 0;�
u

/* -- see zlib.h -- */
int ZAXPORT gzvead(file( buf, len)
   �gzFile file;
    voidp buf;
    unsigned!len;
{
    unsigned gOt, n;	
    gz_statep state;
    {_streamp strm;

    /* get internal`structure */
    if hfile =0NULL(
        return -1;
    statE = (gz_statep)file;
    strm = &(state->strm);

    .*$check that we're reading and that there's no error */
    if state->mkde != GZ_READ || state->err !=&Z_OK)
        return �3;

    /* sincE an int is returned,0make sube len fivs i. one, Otherwise return
       wIth `n efror (this av�ids the flaw in thd interface) */
    if 8(int)len < 0) {
        gz_error(state, Z_BUF_ERROR, "requeste� length does not f�t in int");        return -1:
    }

 �  /* if len is zero,!avoid unnecessary operations */
    if (len == 0)
        rgturn 0;

    /* process a skip request */
    if (stqte->seek) {
        state->seec � 0;
     `  if (gz_skir(state, state->skmp) == -1)
            return -�;
    }
    +* get`le� bytes to buf, or Less than len if at the end */
0   got = 0;
    do {
       "/* fi2st just try co0yanw data nrom the output buffer */
        ib (state->have) {
            n = s4ate->have 6 len$? l%. : qtaTe->have;
" !         memcpy(buf, state->next, j);
    "       state->next += n;
            state->have -= n;
  `     }

        /: /utput buffer empty -- returN if wegre at the end of the input */
  "   ( else if (state->eof && strm->avail_in == 0-   �        break;

        /* need output dat` -- for small len or new stream load qp our output
           buffer */
        else if (statg->how == LOOK || lel < (state->siz% << 1i) {J       (    /* get more output, looking for header if required */
            if (g~_ma{e(st�te+ == -1)
      �         return -1;
            continue;    !  /* no progresq yet -- go back to memcpy() abovE */
            /. the copy above assu�es that we will leavg w)t� space(in�the `        `   output Buffer, allouin' at least k�e gzungetc() to succeed */
        }

   ! `  /* large len - read direat,y indo user buffer ./
   !    else if (state->how == COPYi {      /* read directly */�
     �      if�(gz_,oad(state, buF, len, &n) == -1)    $!!         return -1+
        }

        /* larga len - decompress directly`knto user buffmr */
        else {  /* state->how == GZIP */
            strm->avail_out = len;
    `       strm-:next_out = buf�
        `   if (gzOdecomp(qtape) = -1)
                veturn m1;
2 0    0    n = state->h`ve;
     0      state->have = 0;        }

    "   /* update progresS */
        len -= n;
   (   $buf = (chap *)buf!+ n:
        got += n;
        sda�e->pow +} n;
 "  } while (len);

    /* return number of bytes read into user buffer (will fit in int) */
    return (int)got;
}�
/* -- see zdif.h -- */
int ZEXPOBT gzgetc(fILe)
    gzFile file;
k
    int ret;
    unsigned char buf[1];
    gz_statep state;

    /* eet internal structure */
 0 "if"(f)le == NULL)
       (return -1;
    state = (gz_stAtep)file;

  ( /* cjgck that we're0reading `nd that phere's no errr */
    if (state->mode != GZOREAD || statg->err != Z_OK)
        ret�rn -1;
    /* try output bubfer (ng need to check for sk�p request) */
    if (statm->hava) {
        state->hAve--;J    `  `state->pos++;
       (return *(state->next)++;
    }

    /* nothing there /- try ezread() */
    ret = gzread(file, buf, 1);
    return ret < 1 ? -1 : buf[0];
}

/* -- see zlib.h -- */
int ZEXPORT gzungetc(c, file)
    int c;
    gzFile file;
{
    gz_statep state;

    /* get internal structure */
    if (file == NULL)
        return -1;
    state = (gz_statep)file;

    /* check that we're reading and that there's no error */
    if (state->mode != GZ_READ || state->err != Z_OK)
        return -1;

    /* process a skip request */
    if (state->seek) {
        state->seek = 0;
        if (gz_skip(state, state->skip) == -1)
            return -1;
    }

    /* can't push EOF */
    if (c < 0)
        return -1;

    /* if output buffer empty, put byte at end (allows more pushing) */
    if (state->have == 0) {
        state->have = 1;
        state->next = state->out + (state->size << 1) - 1;
        state->next[0] = c;
        state->pos--;
        return c;
    }

    /* if no room, give up (must have already done a gzungetc()) */
    if (state->have == (state->size << 1)) {
        gz_error(state, Z_BUF_ERROR, "out of room to push characters");
        return -1;
    }

    /* slide output data if needed and insert byte before existing data */
    if (state->next == state->out) {
        unsigned char *src = state->out + state->have;
        unsigned char *dest = state->out + (state->size << 1);
        while (src > state->out)
            *--dest = *--src;
        state->next = dest;
    }
    state->have++;
    state->next--;
    state->next[0] = c;
    state->pos--;
    return c;
}

/* -- see zlib.h -- */
char * ZEXPORT gzgets(file, buf, len)
    gzFile file;
    char *buf;
    int len;
{
    unsigned left, n;
    char *str;
    unsigned char *eol;
    gz_statep state;

    /* check parameters and get internal structure */
    if (file == NULL || buf == NULL || len < 1)
        return NULL;
    state = (gz_statep)file;

    /* check that we're reading and that there's no error */
    if (state->mode != GZ_READ || state->err != Z_OK)
        return NULL;

    /* process a skip request */
    if (state->seek) {
        state->seek = 0;
        if (gz_skip(state, state->skip) == -1)
            return NULL;
    }

    /* copy output bytes up to new line or len - 1, whichever comes first --
       append a terminating zero to the string (we don't check for a zero in
       the contents, let the user worry about that) */
    str = buf;
    left = (unsigned)len - 1;
    if (left) do {
        /* assure that something is in the output buffer */
        if (state->have == 0) {
            if (gz_make(state) == -1)
                return NULL;            /* error */
            if (state->have == 0) {     /* end of file */
                if (buf == str)         /* got bupkus */
                    return NULL;
                break;                  /* got something -- return it */
            }
        }

        /* look for end-of-line in current output buffer */
        n = state->have > left ? left : state->have;
        eol = memchr(state->next, '\n', n);
        if (eol != NULL)
            n = (unsigned)(eol - state->next) + 1;

        /* copy through end-of-line, or remainder if not found */
        memcpy(buf, state->next, n);
        state->have -= n;
        state->next += n;
        state->pos += n;
        left -= n;
        buf += n;
    } while (left && eol == NULL);

    /* found end-of-line or out of space -- terminate string and return it */
    buf[0] = 0;
    return str;
}

/* -- see zlib.h -- */
int ZEXPORT gzdirect(file)
    gzFile file;
{
    gz_statep state;

    /* get internal structure */
    if (file == NULL)
        return 0;
    state = (gz_statep)file;

    /* check that we're reading */
    if (state->mode != GZ_READ)
        return 0;

    /* if the state is not known, but we can find out, then do so (this is
       mainly for right after a gzopen() or gzdopen()) */
    if (state->how == LOOK && state->have == 0)
        (void)gz_head(state);

    /* return 1 if reading direct, 0 if decompressing a gzip stream */
    return state->direct;
}

/* -- see zlib.h -- */
int ZEXPORT gzclose_r(file)
    gzFile file;
{
    int ret;
    gz_statep state;

    /* get internal structure */
    if (file == NULL)
        return Z_STREAM_ERROR;
    state = (gz_statep)file;

    /* check that we're reading */
    if (state->mode != GZ_READ)
        return Z_STREAM_ERROR;

    /* free memory and close file */
    if (state->size) {
        inflateEnd(&(state->strm));
        free(state->out);
        free(state->in);
    }
    gz_error(state, Z_OK, NULL);
    free(state->path);
    ret = close(state->fd);
    free(state);
    return ret ? Z_ERRNO : Z_OK;
}
