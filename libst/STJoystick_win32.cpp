/* libST Joystick Implementation
 * Windows Specific Code
 * Adapted directly from libSDL, uses Windows Multimedia
 * library to access the joysticks.
 */

/*
    SDL - Simple DirectMedia Layer
    Copyright (C) 1997-2006 Sam Lantinga

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    Sam Lantinga
    slouken@libsdl.org
*/

/* Win32 MultiMedia Joystick driver, contributed by Andrei de A. Formiga */

#ifndef _ST_JOYSTICK_WIN32_H_
#define _ST_JOYSTICK_WIN32_H_

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <mmsystem.h>
#include <regstr.h>
#include <stdio.h>
#include <stdlib.h>

#define snprintf _snprintf

#include "STJoystick.h"

#define MAX_JOYSTICKS    16
#define MAX_AXES    6    /* each joystick can have up to 6 axes */
#define MAX_BUTTONS    32    /* and 32 buttons                      */
#define AXIS_MIN    -32768  /* minimum value for axis coordinate */
#define AXIS_MAX    32767   /* maximum value for axis coordinate */
/* limit axis to 256 possible positions to filter out noise */
#define JOY_AXIS_THRESHOLD      (((AXIS_MAX)-(AXIS_MIN))/256)
#define JOY_BUTTON_FLAG(n)    (1<<n)


/* array to hold joystick ID values */
static UINT    SYS_JoystickID[MAX_JOYSTICKS];
static JOYCAPS    SYS_Joystick[MAX_JOYSTICKS];
static char    *SYS_JoystickName[MAX_JOYSTICKS];

struct _transaxis {
    int offset;
    float scale;
};

/* The private structure used to keep track of a joystick */
struct joystick_hwdata
{
    /* joystick ID */
    UINT    id;

    /* values used to translate device-specific coordinates into
       SDL-standard ranges */
    _transaxis transaxis[6];
};

/* Convert a win32 Multimedia API return code to a text message */
static void SetMMerror(char *function, int code);

/* Gets joystick name from the registry. */
static char *GetJoystickName(int index, const char *szRegKey)
{
    /* added 7/24/2004 by Eckhard Stolberg */
    /*
        see if there is a joystick for the current
        index (1-16) listed in the registry
    */
    char *name = NULL;
    HKEY hKey;
    DWORD regsize;
    LONG regresult;
    unsigned char regkey[256];
    unsigned char regvalue[256];
    unsigned char regname[256];

    snprintf((char *) regkey, sizeof(regkey), "%s\\%s\\%s",
        REGSTR_PATH_JOYCONFIG,
        szRegKey,
        REGSTR_KEY_JOYCURR);
    regresult = RegOpenKeyExA(HKEY_LOCAL_MACHINE,
        (LPTSTR) &regkey, 0, KEY_READ, &hKey);
    if (regresult == ERROR_SUCCESS)
    {
        /*
            find the registry key name for the
            joystick's properties
        */
        regsize = sizeof(regname);
        snprintf((char *) regvalue, sizeof(regvalue),
            "Joystick%d%s", index+1,
            REGSTR_VAL_JOYOEMNAME);
        regresult = RegQueryValueExA(hKey,
            (char *) regvalue, 0, 0, (LPBYTE) &regname,
            (LPDWORD) &regsize);
        RegCloseKey(hKey);
        if (regresult == ERROR_SUCCESS)
        {
            /* open that registry key */
            snprintf((char *) regkey, sizeof(regkey), "%s\\%s",
                REGSTR_PATH_JOYOEM, regname);
            regresult = RegOpenKeyExA(HKEY_LOCAL_MACHINE,
                (char *) regkey, 0, KEY_READ, &hKey);
            if (regresult == ERROR_SUCCESS)
            {
                /* find the size for the OEM name text */
                regsize = sizeof(regvalue);
                regresult =
                    RegQueryValueExA(hKey,
                    REGSTR_VAL_JOYOEMNAME,
                    0, 0, NULL,
                    (LPDWORD) &regsize);
                if (regresult == ERROR_SUCCESS)
                {
                    /*
                        allocate enough memory
                        for the OEM name text ...
                    */
                    name = new char[regsize];
                    /* ... and read it from the registry */
                    regresult =
                        RegQueryValueExA(hKey,
                        REGSTR_VAL_JOYOEMNAME, 0, 0,
                        (LPBYTE) name,
                        (LPDWORD) &regsize);
                    RegCloseKey(hKey);
                }
            }
        }
    }
    return(name);
}

/* Function to scan the system for joysticks.
 * It should return 0, or -1 on an unrecoverable fatal error.
 */
int STJoystick::Sys_Initialize(void)
{
    int    i;
    int maxdevs;
    int numdevs;
    JOYINFOEX joyinfo;
    JOYCAPS    joycaps;
    MMRESULT result;

    /* Reset the joystick ID & name mapping tables */
    for ( i = 0; i < MAX_JOYSTICKS; ++i ) {
        SYS_JoystickID[i] = 0;
        SYS_JoystickName[i] = NULL;
    }

    /* Loop over all potential joystick devices */
    numdevs = 0;
    maxdevs = joyGetNumDevs();
    for ( i = JOYSTICKID1; i < maxdevs && numdevs < MAX_JOYSTICKS; ++i ) {
        
        joyinfo.dwSize = sizeof(joyinfo);
        joyinfo.dwFlags = JOY_RETURNALL;
        result = joyGetPosEx(SYS_JoystickID[i], &joyinfo);
        if ( result == JOYERR_NOERROR ) {
            result = joyGetDevCaps(i, &joycaps, sizeof(joycaps));
            if ( result == JOYERR_NOERROR ) {
                SYS_JoystickID[numdevs] = i;
                SYS_Joystick[numdevs] = joycaps;
                SYS_JoystickName[numdevs] = GetJoystickName(i, joycaps.szRegKey);
                numdevs++;
            }
        }
    }
    return(numdevs);
}

/* Function to get the device-dependent name of a joystick */
const char* STJoystick::Sys_GetName(int index)
{
    if ( SYS_JoystickName[index] != NULL ) {
        return(SYS_JoystickName[index]);
    } else {
        return(SYS_Joystick[index].szPname);
    }
}

/* Function to open a joystick for use.
   This should fill the nbuttons and naxes fields of the joystick structure.
   It returns 0, or -1 if there is an error.
 */
int STJoystick::Sys_OpenJoystick(int device_index)
{
    int index, i;
    int caps_flags[MAX_AXES-2] =
        { JOYCAPS_HASZ, JOYCAPS_HASR, JOYCAPS_HASU, JOYCAPS_HASV };
    int axis_min[MAX_AXES], axis_max[MAX_AXES];


    /* shortcut */
    index = joysticks[device_index]->index;
    axis_min[0] = SYS_Joystick[index].wXmin;
    axis_max[0] = SYS_Joystick[index].wXmax;
    axis_min[1] = SYS_Joystick[index].wYmin;
    axis_max[1] = SYS_Joystick[index].wYmax;
    axis_min[2] = SYS_Joystick[index].wZmin;
    axis_max[2] = SYS_Joystick[index].wZmax;
    axis_min[3] = SYS_Joystick[index].wRmin;
    axis_max[3] = SYS_Joystick[index].wRmax;
    axis_min[4] = SYS_Joystick[index].wUmin;
    axis_max[4] = SYS_Joystick[index].wUmax;
    axis_min[5] = SYS_Joystick[index].wVmin;
    axis_max[5] = SYS_Joystick[index].wVmax;

    /* allocate memory for system specific hardware data */
    joysticks[device_index]->hwdata = new joystick_hwdata;
    if (joysticks[device_index]->hwdata == NULL)
    {
        STJoystick::MemoryError();
        return(-1);
    }
    memset(joysticks[device_index]->hwdata, 0, sizeof(joystick_hwdata));

    /* set hardware data */
    joysticks[device_index]->hwdata->id = SYS_JoystickID[index];
    for ( i = 0; i < MAX_AXES; ++i ) {
        if ( (i<2) || (SYS_Joystick[index].wCaps & caps_flags[i-2]) ) {
            joysticks[device_index]->hwdata->transaxis[i].offset =
                AXIS_MIN - axis_min[i];
            joysticks[device_index]->hwdata->transaxis[i].scale =
                (float)(AXIS_MAX - AXIS_MIN) / (axis_max[i] - axis_min[i]);
        } else {
            joysticks[device_index]->hwdata->transaxis[i].offset = 0;
            joysticks[device_index]->hwdata->transaxis[i].scale = 1.0; /* Just in case */
        }
    }

    /* fill nbuttons, naxes, and nhats fields */
    joysticks[device_index]->nbuttons = SYS_Joystick[index].wNumButtons;
    joysticks[device_index]->naxes = SYS_Joystick[index].wNumAxes;
    if ( SYS_Joystick[index].wCaps & JOYCAPS_HASPOV ) {
        joysticks[device_index]->nhats = 1;
    } else {
        joysticks[device_index]->nhats = 0;
    }
    return(0);
}

static unsigned char TranslatePOV(DWORD value)
{
    unsigned char pos;

    pos = ST_JOYSTICK_HAT_CENTERED;
    if ( value != JOY_POVCENTERED ) {
        if ( (value > JOY_POVLEFT) || (value < JOY_POVRIGHT) ) {
            pos |= ST_JOYSTICK_HAT_UP;
        }
        if ( (value > JOY_POVFORWARD) && (value < JOY_POVBACKWARD) ) {
            pos |= ST_JOYSTICK_HAT_RIGHT;
        }
        if ( (value > JOY_POVRIGHT) && (value < JOY_POVLEFT) ) {
            pos |= ST_JOYSTICK_HAT_DOWN;
        }
        if ( value > JOY_POVBACKWARD ) {
            pos |= ST_JOYSTICK_HAT_LEFT;
        }
    }
    return(pos);
}

/* Function to update the state of a joystick - called as a device poll. */
void STJoystick::Sys_Update()
{
    MMRESULT result;
    int i;
    DWORD flags[MAX_AXES] = { JOY_RETURNX, JOY_RETURNY, JOY_RETURNZ, 
                  JOY_RETURNR, JOY_RETURNU, JOY_RETURNV };
    DWORD pos[MAX_AXES];
    struct _transaxis *transaxis;
    int value, change;
    JOYINFOEX joyinfo;

    joyinfo.dwSize = sizeof(joyinfo);
    joyinfo.dwFlags = JOY_RETURNALL|JOY_RETURNPOVCTS;
    if ( ! joystick_data->hats ) {
        joyinfo.dwFlags &= ~(JOY_RETURNPOV|JOY_RETURNPOVCTS);
    }
    result = joyGetPosEx(joystick_data->hwdata->id, &joyinfo);
    if ( result != JOYERR_NOERROR ) {
        SetMMerror("joyGetPosEx", result);
        return;
    }

    /* joystick motion events */
    pos[0] = joyinfo.dwXpos;
    pos[1] = joyinfo.dwYpos;
    pos[2] = joyinfo.dwZpos;
    pos[3] = joyinfo.dwRpos;
    pos[4] = joyinfo.dwUpos;
    pos[5] = joyinfo.dwVpos;

    transaxis = (struct _transaxis*)joystick_data->hwdata->transaxis;
    for (i = 0; i < joystick_data->naxes; i++) {
        if (joyinfo.dwFlags & flags[i]) {
            value = (int)(((float)pos[i] + transaxis[i].offset) * transaxis[i].scale);
            change = (value - joystick_data->axes[i]);
            if ( (change < -JOY_AXIS_THRESHOLD) || (change > JOY_AXIS_THRESHOLD) ) {
                joystick_data->axes[i] = value;
            }
        }
    }

    /* joystick button events */
    if ( joyinfo.dwFlags & JOY_RETURNBUTTONS ) {
        for ( i = 0; i < joystick_data->nbuttons; ++i ) {
            if ( joyinfo.dwButtons & JOY_BUTTON_FLAG(i) ) {
                if ( joystick_data->buttons[i] != ST_JOYSTICK_PRESSED ) {
                    joystick_data->buttons[i] = ST_JOYSTICK_PRESSED;
                }
            } else {
                if ( joystick_data->buttons[i] != ST_JOYSTICK_RELEASED ) {
                    joystick_data->buttons[i] = ST_JOYSTICK_RELEASED;
                }
            }
        }
    }

    /* joystick hat events */
    if ( joyinfo.dwFlags & JOY_RETURNPOV ) {
        unsigned char pos;

        pos = TranslatePOV(joyinfo.dwPOV);
        if ( pos != joystick_data->hats[0] ) {
            joystick_data->hats[0] = pos;
        }
    }
}

/* Function to close a joystick after use */
void STJoystick::Sys_Close(int device_index)
{
    if (joysticks[device_index]->hwdata != NULL) {
        /* free system specific hardware data */
        delete joysticks[device_index]->hwdata;
    }
}

/* Function to perform any system-specific joystick related cleanup */
void STJoystick::Sys_Destroy()
{
    int i;
    for (i = 0; i < MAX_JOYSTICKS; i++) {
        if ( SYS_JoystickName[i] != NULL ) {
            delete SYS_JoystickName[i];
        }
    }
}


/* WinMM Error Handler */
void SetMMerror(char *function, int code)
{
    static char *error;
    static char  errbuf[1024];

    errbuf[0] = 0;
    switch (code) 
    {
        case MMSYSERR_NODRIVER:
            error = "Joystick driver not present";
        break;

        case MMSYSERR_INVALPARAM:
        case JOYERR_PARMS:
            error = "Invalid parameter(s)";
        break;
        
        case MMSYSERR_BADDEVICEID:
            error = "Bad device ID";
        break;

        case JOYERR_UNPLUGGED:
            error = "Joystick not attached";
        break;

        case JOYERR_NOCANDO:
            error = "Can't capture joystick input";
        break;

        default:
            snprintf(errbuf, sizeof(errbuf),
                     "%s: Unknown Multimedia system error: 0x%x",
                                function, code);
        break;
    }

    if ( ! errbuf[0] ) {
        snprintf(errbuf, sizeof(errbuf), "%s: %s", function, error);
    }
    STJoystick::FatalError(errbuf);
}

#endif /* _ST_JOYSTICK_WIN32_H_ */
