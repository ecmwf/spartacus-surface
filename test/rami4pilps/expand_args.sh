#!/bin/bash
#
# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

INPUT=$(echo $1 | awk -F_ '{print $1}')
BAND=$(echo $INPUT | awk -F- '{print $1}')
SURF=$(echo $INPUT | awk -F- '{print $2}')
FRAC=$(echo $INPUT | awk -F- '{print $3}')
VREGS=$(echo $INPUT | awk -F- '{print $4}')
STREAMS=$(echo $INPUT | awk -F- '{print $5}')

if [ "$BAND" = vis ]; then
    SSA=0.1301
    if [ "$SURF" = med ]; then
	ALBEDO=0.1217
    elif [ "$SURF" = snw ]; then
	ALBEDO=0.9640
    else
	echo "Surface \"$SURF\" not understood"
	return 1
    fi
elif [ "$BAND" = nir ]; then
    SSA=0.8058
    if [ "$SURF" = med ]; then
	ALBEDO=0.2142
    elif [ "$SURF" = snw ]; then
	ALBEDO=0.5568
    else
	echo "Surface \"$SURF\" not understood"
	return 1
    fi
else
    echo "Band \"$BAND\" not understood"
fi

if [ "$VREGS" ]; then
    VREG="n_vegetation_region_forest=$VREGS"
fi

if [ "$STREAMS" ]; then
    STR="n_stream_sw_forest=$STREAMS"
else
    STR=""
fi

echo vegetation_fraction=$FRAC ground_sw_albedo=$ALBEDO vegetation_sw_ssa=$SSA $VREG $STR
