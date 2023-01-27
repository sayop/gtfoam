#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2017 by Sayop Kim
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     etc/config/settings.sh
#
# Description
#     Startup file for OpenFOAM
#     Sourced from OpenFOAM-<VERSION>/etc/bashrc
#
#------------------------------------------------------------------------------

# prefix to PATH
_foamAddPath()
{
    while [ $# -ge 1 ]
    do
        export PATH=$1:$PATH
        shift
    done
}

# prefix to LD_LIBRARY_PATH
_foamAddLib()
{
    while [ $# -ge 1 ]
    do
        export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
        shift
    done
}

# Location of GTFOAM source files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export GT_SRC=$GT_PROJECT_DIR/src

# Location of GTFOAM application files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export GT_APP=$GT_PROJECT_DIR/applications

# Location of GTFOAM solver files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export GT_SOL=$GT_PROJECT_DIR/solvers


# base executables/libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
export GT_APPBIN=$GT_PROJECT_DIR/platforms/$WM_OPTIONS/bin
export GT_LIBBIN=$GT_PROJECT_DIR/platforms/$WM_OPTIONS/lib

# Add GTFOAM app binary files into $PATH
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_foamAddPath $GT_APPBIN

# Add GTFOAM libs into $LD_LIBRARY_PATH
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_foamAddLib $GT_LIBBIN


# ----------------------------------------------------------------- end-of-file
