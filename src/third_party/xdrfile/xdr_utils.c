#include <stdlib.h>
#include <stdio.h>
#include "xdr_utils.h"
#include "xdr_seek.h"
#include "trr_header.h"

#ifndef XTC_MAGIC
#    define XTC_MAGIC 1995
#endif
#define XDR_INT_SIZE 4
const int header_size = 16;

int xdr_xtc_get_natoms(XDRFILE* handle, int* natoms){
    int magic,ok;
    int64_t pos = xdr_tell(handle); // save pos
    ok = xdrfile_read_int(&magic,1,handle);
    if(!ok) return 0;
    int ret = xdrfile_read_int(natoms,1,handle);
    if(!ok) return 0;
    xdr_seek(handle,pos,SEEK_SET); // Rewind
    return 1;
}

int xtc_at_header_start(XDRFILE* handle, int natoms, int* timestep, float* time)
{
    int       i_inp[3];
    float     f_inp[10];
    int       i;
    int64_t   off;


    if ((off = xdr_tell(handle)) < 0)
    {
        return -1;
    }
    /* read magic natoms and timestep */
    for (i = 0; i < 3; i++)
    {
        if (!xdrfile_read_int(&(i_inp[i]),1,handle))
        {
            xdr_seek(handle, off + XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* quick return */
    if (i_inp[0] != XTC_MAGIC)
    {
        if (xdr_seek(handle, off + XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        return 0;
    }
    /* read time and box */
    for (i = 0; i < 10; i++)
    {
        if (!xdrfile_read_float(&(f_inp[i]),1,handle))
        {
            xdr_seek(handle, off + XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* Make a rigourous check to see if we are in the beggining of a header
       Hopefully there are no ambiguous cases */
    /* This check makes use of the fact that the box matrix has 3 zeroes on the upper
       right triangle and that the first element must be nonzero unless the entire matrix is zero
     */
    if (i_inp[1] == natoms
        && ((f_inp[1] != 0 && f_inp[6] == 0) || (f_inp[1] == 0 && f_inp[5] == 0 && f_inp[9] == 0)))
    {
        if (xdr_seek(handle, off + XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        *time     = f_inp[0];
        *timestep = i_inp[2];
        return 1;
    }
    if (xdr_seek(handle, off + XDR_INT_SIZE, SEEK_SET))
    {
        return -1;
    }
    return 0;
}

int xtc_get_current_frame_number(XDRFILE* handle, int natoms, bool* bOK)
{
    int64_t   off;
    int       ret;
    int       step;
    float     time;
    *bOK = false;

    if ((off = xdr_tell(handle)) < 0)
    {
        return -1;
    }


    while (true)
    {
        ret = xtc_at_header_start(handle, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = true;
            if (xdr_seek(handle, off, SEEK_SET))
            {
                *bOK = false;
                return -1;
            }
            return step;
        }
        else if (ret == -1)
        {
            if (xdr_seek(handle, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (xdr_seek(handle, -2 * XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
}


int xtc_get_next_frame_number(XDRFILE* handle, int natoms)
{
    int64_t   off;
    int       step;
    float     time;
    int       ret;

    if ((off = xdr_tell(handle)) < 0)
    {
        return -1;
    }

    /* read one int just to make sure we dont read this frame but the next */
    xdrfile_read_int(&step,1,handle);
    while (true)
    {
        ret = xtc_at_header_start(handle, natoms, &step, &time);
        if (ret == 1)
        {
            if (xdr_seek(handle, off, SEEK_SET))
            {
                return -1;
            }
            return step;
        }
        else if (ret == -1)
        {
            if (xdr_seek(handle, off, SEEK_SET))
            {
                return -1;
            }
        }
    }
}

float xtc_get_current_frame_time(XDRFILE* handle, int natoms, bool* bOK)
{
    int64_t   off;
    int       step;
    float     time;
    int       ret;
    *bOK = false;

    if ((off = xdr_tell(handle)) < 0)
    {
        return -1;
    }

    while (true)
    {
        ret = xtc_at_header_start(handle, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = true;
            if (xdr_seek(handle, off, SEEK_SET))
            {
                *bOK = false;
                return -1;
            }
            return time;
        }
        else if (ret == -1)
        {
            if (xdr_seek(handle, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (xdr_seek(handle, -2 * XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
}

float xtc_get_next_frame_time(XDRFILE* handle, int natoms, bool* bOK)
{
    int64_t   off;
    float     time;
    int       step;
    int       ret;
    *bOK = false;

    if ((off = xdr_tell(handle)) < 0)
    {
        return -1;
    }
    /* read one int just to make sure we dont read this frame but the next */
    xdrfile_read_int(&step,1,handle);
    while (true)
    {
        ret = xtc_at_header_start(handle, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = true;
            if (xdr_seek(handle, off, SEEK_SET))
            {
                *bOK = false;
                return -1;
            }
            return time;
        }
        else if (ret == -1)
        {
            if (xdr_seek(handle, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
    }
}

float xdr_xtc_estimate_dt(XDRFILE* handle, int natoms, bool* bOK)
{
    float     res;
    float     tinit;
    int64_t   off;

    *bOK = false;
    if ((off = xdr_tell(handle)) < 0)
    {
        return -1;
    }

    tinit = xtc_get_current_frame_time(handle, natoms, bOK);

    if (!(*bOK))
    {
        return -1;
    }

    res = xtc_get_next_frame_time(handle, natoms, bOK);

    if (!(*bOK))
    {
        return -1;
    }

    res -= tinit;
    if (0 != xdr_seek(handle, off, SEEK_SET))
    {
        *bOK = false;
        return -1;
    }
    return res;
}

int xdr_xtc_get_last_frame_number(XDRFILE* handle, int natoms, bool* bOK)
{
    int       frame;
    int64_t off;
    *bOK = true;

    if ((off = xdr_tell(handle)) < 0)
    {
        *bOK = false;
        return -1;
    }


    if (xdr_seek(handle, -3 * XDR_INT_SIZE, SEEK_END))
    {
        *bOK = false;
        return -1;
    }

    frame = xtc_get_current_frame_number(handle, natoms, bOK);
    if (!*bOK)
    {
        return -1;
    }


    if (xdr_seek(handle, off, SEEK_SET))
    {
        *bOK = false;
        return -1;
    }

    return frame;
}

float xdr_xtc_get_last_frame_time(XDRFILE* handle, int natoms, bool* bOK)
{
    float     time;
    int64_t   off;
    *bOK = true;
    off  = xdr_tell(handle);
    if (off < 0)
    {
        *bOK = false;
        return -1;
    }

    if (xdr_seek(handle, -3 * XDR_INT_SIZE, SEEK_END) != 0)
    {
        *bOK = false;
        return -1;
    }

    time = xtc_get_current_frame_time(handle, natoms, bOK);
    if (!(*bOK))
    {
        return -1;
    }

    if (xdr_seek(handle, off, SEEK_SET) != 0)
    {
        *bOK = false;
        return -1;
    }
    return time;
}

int64_t xtc_get_next_frame_start(XDRFILE* handle, int natoms)
{
    int64_t   res;
    int       ret;
    int       step;
    float     time;
    /* read one int just to make sure we dont read this frame but the next */
    xdrfile_read_int(&step,1,handle);
    while (true)
    {
        ret = xtc_at_header_start(handle, natoms, &step, &time);
        if (ret == 1)
        {
            if ((res = xdr_tell(handle)) >= 0)
            {
                return res - XDR_INT_SIZE;
            }
            else
            {
                return res;
            }
        }
        else if (ret == -1)
        {
            return -1;
        }
    }
}


int xdr_xtc_seek_frame(int frame, XDRFILE* handle, int natoms)
{
    int64_t low = 0;
    int64_t high, pos;


    /* round to 4 bytes */
    int       fr;
    int64_t offset;
    if (xdr_seek(handle, 0, SEEK_END))
    {
        return -1;
    }

    if ((high = xdr_tell(handle)) < 0)
    {
        return -1;
    }

    /* round to 4 bytes  */
    high /= XDR_INT_SIZE;
    high *= XDR_INT_SIZE;
    offset = ((high / 2) / XDR_INT_SIZE) * XDR_INT_SIZE;

    if (xdr_seek(handle, offset, SEEK_SET))
    {
        return -1;
    }

    while (true)
    {
        fr = xtc_get_next_frame_number(handle, natoms);
        if (fr < 0)
        {
            return -1;
        }
        if (fr != frame && llabs(low - high) > header_size)
        {
            if (fr < frame)
            {
                low = offset;
            }
            else
            {
                high = offset;
            }
            /* round to 4 bytes */
            offset = (((high + low) / 2) / 4) * 4;

            if (xdr_seek(handle, offset, SEEK_SET))
            {
                return -1;
            }
        }
        else
        {
            break;
        }
    }
    if (offset <= header_size)
    {
        offset = low;
    }

    if (xdr_seek(handle, offset, SEEK_SET))
    {
        return -1;
    }

    if ((pos = xtc_get_next_frame_start(handle, natoms)) < 0)
    {
        /* we probably hit an end of file */
        return -1;
    }

    if (xdr_seek(handle, pos, SEEK_SET))
    {
        return -1;
    }

    return 0;
}

int xdr_xtc_seek_time(float time, XDRFILE* handle, int natoms, bool bSeekForwardOnly)
{
    float     t;
    float     dt;
    bool  bOK = false;
    int64_t low = 0;
    int64_t high, offset, pos;
    int       dt_sign = 0;

    if (bSeekForwardOnly)
    {
        low = xdr_tell(handle) - header_size;
    }
    if (xdr_seek(handle, 0, SEEK_END))
    {
        return -1;
    }

    if ((high = xdr_tell(handle)) < 0)
    {
        return -1;
    }
    /* round to int  */
    high /= XDR_INT_SIZE;
    high *= XDR_INT_SIZE;
    offset = (((high - low) / 2) / XDR_INT_SIZE) * XDR_INT_SIZE;

    if (xdr_seek(handle, offset, SEEK_SET))
    {
        return -1;
    }


    /*
     * No need to call xdr_xtc_estimate_dt here - since xdr_xtc_estimate_dt is called first thing in
     the loop dt = xdr_xtc_estimate_dt(fp, xdrs, natoms, &bOK);

       if (!bOK)
       {
        return -1;
       }
     */

    while (true)
    {
        dt = xdr_xtc_estimate_dt(handle, natoms, &bOK);
        if (!bOK)
        {
            return -1;
        }
        else
        {
            if (dt > 0)
            {
                if (dt_sign == -1)
                {
                    /* Found a place in the trajectory that has positive time step while
                       other has negative time step */
                    return -2;
                }
                dt_sign = 1;
            }
            else if (dt < 0)
            {
                if (dt_sign == 1)
                {
                    /* Found a place in the trajectory that has positive time step while
                       other has negative time step */
                    return -2;
                }
                dt_sign = -1;
            }
        }
        t = xtc_get_next_frame_time(handle, natoms, &bOK);
        if (!bOK)
        {
            return -1;
        }

        /* If we are before the target time and the time step is positive or 0, or we have
           after the target time and the time step is negative, or the difference between
           the current time and the target time is bigger than dt and above all the distance between
           high and low is bigger than 1 frame, then do another step of binary search. Otherwise
           stop and check if we reached the solution */
        if ((((t < time && dt_sign >= 0) || (t > time && dt_sign == -1))
             || ((t - time) >= dt && dt_sign >= 0) || ((time - t) >= -dt && dt_sign < 0))
            && (llabs(low - high) > header_size))
        {
            if (dt >= 0 && dt_sign != -1)
            {
                if (t < time)
                {
                    low = offset;
                }
                else
                {
                    high = offset;
                }
            }
            else if (dt <= 0 && dt_sign == -1)
            {
                if (t >= time)
                {
                    low = offset;
                }
                else
                {
                    high = offset;
                }
            }
            else
            {
                /* We should never reach here */
                return -1;
            }
            /* round to 4 bytes and subtract header*/
            offset = (((high + low) / 2) / XDR_INT_SIZE) * XDR_INT_SIZE;
            if (xdr_seek(handle, offset, SEEK_SET))
            {
                return -1;
            }
        }
        else
        {
            if (llabs(low - high) <= header_size)
            {
                break;
            }
            /* re-estimate dt */
            if (xdr_xtc_estimate_dt(handle, natoms, &bOK) != dt)
            {
                if (bOK)
                {
                    dt = xdr_xtc_estimate_dt(handle, natoms, &bOK);
                }
            }
            if (t >= time && t - time < dt)
            {
                break;
            }
        }
    }

    if (offset <= header_size)
    {
        offset = low;
    }

    xdr_seek(handle, offset, SEEK_SET);

    if ((pos = xtc_get_next_frame_start(handle, natoms)) < 0)
    {
        return -1;
    }

    if (xdr_seek(handle, pos, SEEK_SET))
    {
        return -1;
    }
    return 0;
}


// Check if forces and velocities are present in TRR file
int check_trr_content(XDRFILE* handle, int* natoms, int* xsz, int* vsz, int* fsz)
{
    t_trnheader sh;
    int  ret = do_trnheader(handle,1,&sh);
    if(ret != exdrOK) return ret;

    *natoms = sh.natoms;
    *xsz = sh.x_size;
    *vsz = sh.v_size;
    *fsz = sh.f_size;

    return exdrOK;
}
